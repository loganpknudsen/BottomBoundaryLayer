### Load in Packages
using Pkg
Pkg.activate(".")
# Pkg.instantiate()
using Oceananigans
using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.AbstractOperations: @at, Average
using Oceananigans.Grids: Center, Face
using Oceananigans.Fields: FunctionField
using Oceananigans.Units
using Oceananigans.OutputWriters: Checkpointer
using Random
using Printf
using ArgParse
# using CUDA: has_cuda_gpu
using CUDA 
using Oceanostics

# Path file is saved under
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "path"
        help = "pathname to save data under"
        default = ""
    end
    return parse_args(s)
end

args=parse_commandline()

@info("command line args:")
for (arg,val) in args
  @info("   $arg => $val")
end

path_name = args["path"]


# grid specifications
arch = has_cuda_gpu() ? GPU() : CPU()
@info("Arch => $arch")

Lx = 2000meters
Lz = 200meters
Nx = 500
Nz = 100 # Note to self, maintain 2 to 1 resolution ration
# arch;
grid = RectilinearGrid(arch; topology = (Periodic, Flat, Bounded),
                       size = (Nx, Nz),
                       x = (0, Lx),
                       z = (0,Lz))


# tilted domain parameters
θ = 5*10^(-3) # degrees 10^(-2) is previous value for 110 meter layer
ĝ = [sind(θ), 0, cosd(θ)] # gravity vector

# realistic mid latitude for now
buoyancy = Buoyancy(model = BuoyancyTracer(), gravity_unit_vector = -ĝ)
coriolis = ConstantCartesianCoriolis(f = 1e-4, rotation_axis = ĝ)

# parameters for simulation
const V∞ = 0.05 # m s⁻¹ interior velocity
const f = 1e-4  # coriolis parameter
const N² = 1e-5 # interior stratification
const S∞ = (N²*θ^2)/(f^2) # sloep burger number
const γ = (1+S∞)^(-1) # 0 PV parameter
const hu = ceil((f*V∞)/(γ*N²*θ)) # Height of Boundary Layer
const fˢ=(f^2+θ^2*N²)^(0.5) # modified oscillation
const uₒ = 0 # Initial u shear perturbation
const vₒ = γ*(N²*θ)/(f)*0.5 # Initial v shear perturbation
const bₒ = vₒ*((θ*N²)/(f))*0.1 # initial stratification perturbation
const a1 = (f*vₒ+bₒ*θ)/(fˢ) # a1-h1 are constants for the following oscillations, calculate here for efficiency
const b1 = (f^2*vₒ+f*bₒ*θ)/(fˢ)^2
const c1 = (f*uₒ)/(fˢ)
const d1 = ((fˢ^2-f^2)*vₒ-f*bₒ*θ)/(fˢ)^2
const e1 = N²*θ*(f*vₒ+bₒ*θ)/(fˢ)^2
const h1 = (N²*θ*uₒ)/(fˢ)

# array of paramerers for background function
p =(; N², θ, f, V∞, hu, γ, uₒ, vₒ, bₒ, fˢ, a1, b1, c1, d1, e1, h1)

# background flow with geostrophic and ageostrophic shear 

# heaviside function for boundary layer
@inline heaviside(x,z) = ifelse(z < 0, zero(z), one(z))
# @inline heavisideh(x,z) = ifelse(hu < z, zero(z), one(z))
# @inline heavisideu(x,z,t,p) = ifelse(p.hu < z, zero(x), u_pert(x,z,t,p))
# @inline heavisidev(x,z,t,p) = ifelse(p.hu < z, zero(x), v_pert(x,z,t,p))
# @inline heavisideb(x,z,t,p) = ifelse(x < 0, zero(x), b_pert(x,z,t,p))


# oscillation functions for background
@inline sn_fn(x,z,t,p) = sin(p.fˢ*t)
@inline cs_fn(x,z,t,p) = cos(p.fˢ*t)

u_pert(x,z,t,p) = p.uₒ*cs_fn(x,z,t,p) +p.a1*sn_fn(x,z,t,p) # shear
v_pert(x,z,t,p) = p.b1*cs_fn(x,z,t,p) - p.c1*sn_fn(x,z,t,p)+p.d1
b_pert(x,z,t,p) = p.e1*cs_fn(x,z,t,p) - p.h1*sn_fn(x,z,t,p)+p.bₒ-p.e1

u_adjustment(x, z, t, p) = u_pert(x,z,t,p)*(p.hu-z)*heaviside(x,p.hu-z)
v_adjustment(x, z, t, p) = p.V∞-p.γ*(p.θ * p.N²)/(p.f)*(p.hu-z)*heaviside(x,p.hu-z) + v_pert(x,z,t,p)*(p.hu-z)*heaviside(x,p.hu-z)
constant_stratification(x, z, t, p) = p.N²*x*p.θ + p.N²*z + p.N²*p.γ*(p.hu-z)*heaviside(x,p.hu-z) + b_pert(x,z,t,p)*(p.hu-z)*heaviside(x,p.hu-z)

U_field = BackgroundField(u_adjustment, parameters=p)
V_field = BackgroundField(v_adjustment, parameters=p)
B_field = BackgroundField(constant_stratification, parameters=p)

# Boundary condition set up
# Free Slip Boundary Conditions

b_bc_top= GradientBoundaryCondition(N²)
# b_bc_bottom= GradientBoundaryCondition(-1*N²*(1-γ))
buoyancy_grad = FieldBoundaryConditions(top=b_bc_top) # ,bottom=b_bc_bottom

# diffusitivity and viscosity values for closure
const ν = 1e-4
closure = ScalarDiffusivity(ν=ν, κ=1e-4)

start_time = time_ns()

# model set up 
model = NonhydrostaticModel(; grid, buoyancy, coriolis, closure,
                            timestepper = :RungeKutta3,
                            advection = WENO(),
                            tracers = :b,
                            boundary_conditions = (; b=buoyancy_grad),
                            background_fields = (; u=U_field, v=V_field, b=B_field))

ns = 10^(-4) # standard deviation for noise

# initial conditions to start instability
ui(x, z) = ns*Random.randn()
vi(x, z) = ns*Random.randn()
wi(x, z) = ns*Random.randn()
# bₒ(x,y,z) = 0.005*Random.randn()

# set simulation and decide run time
set!(model, u=ui, v=vi, w=wi)

simulation = Simulation(model, Δt = 1seconds, stop_time = 10.01*(2*pi)/f) # stop_iteration=10

# time step wizard
wizard = TimeStepWizard(cfl=0.95, max_change=1.1seconds, max_Δt=100.0seconds, min_Δt=0.01seconds) 
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5)) 

simulation.output_writers[:checkpointer] = Checkpointer(model; schedule=TimeInterval((5*(2*pi)/f)seconds), prefix="model_checkpoint")

# progress message to check time step size
progress_message(sim) =
        @printf("i: %04d, t: %s, Δt: %s, wall time: %s, Intertial Period %s\n",
        sim.model.clock.iteration, prettytime(sim.model.clock.time),
        prettytime(sim.Δt), prettytime((time_ns() - start_time) * 1e-9),sim.model.clock.time*f/2π)

simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(1000) ) # TimeInterval(0.5*(2*pi)/f) 

# diagnostic calculations, it is saved in 2 files with one saving the flow field and the other tke diagnostics
# @inline bottom_mask(x,z)= heaviside(hu-z) # mask functions to contrict to boundary layer
# full_mask(x,z) = bottom_mask(x,z)
# umask = Oceananigans.Fields.FunctionField{Center,Center,Center}(full_mask,model.grid)
# vmask = Oceananigans.Fields.FunctionField{Center,Face,Center}(full_mask,model.grid)
# calculate the pertubation in velocity
# ua, va, wa = model.velocities # change back to ua, va, wa
# um = Field(@at (Center, Center, Center) Average(ua, dims=1))
# vm = Field(@at (Center, Center, Center) Average(va, dims=1))
# wm = Field(@at (Center, Center, Center) Average(wa, dims=1))
# u = ua - um
# v = va - vm
# w = wa - wm
# # background fields
# ub = model.background_fields.velocities.u
# vb = model.background_fields.velocities.v
# ba = model.tracers.b
# bm = Field(@at (Center, Center, Center) Average(ba, dims=1))
# b = ba - bm
# B∞ = model.background_fields.tracers.b
# # pr = model.pressures.pHY′
# # wapr = wa*pr
# # wmpm = Field(@at (Center, Center, Center) Average(wapr, dims=1))
# # wp = wa*pr - wmpm

# U = u + ub
# V = v + vb
# B = b + B∞

# dbdz = Field(@at (Center, Center, Center) ∂z(b)) #stratification pertubation calculation
# dBdz = Field(@at (Center, Center, Center) ∂z(b+B∞)) # stratification total calculation
# PV = ErtelPotentialVorticity(model, add_background=true) # potential vorticity calculation
# # KE = KineticEnergy(model) # total kinetic energy calculation
# E = KineticEnergyDissipationRate(model; U = u, V = v, W = w) # kinetic energy dissaption calcualtion
# k = 0.5*(u^2 + v^2 + w^2) # TKE calculation
# # waka = wa*ka
# # wmkm = Field(@at (Center, Center, Center) Average(waka, dims=1))
# # wk = wa*ka - wmkm
# # uh = u - θ*w
# # wh = w + θ*u
# # uz = Field(@at (Center, Center, Center) ∂z(u)) 
# # vz = Field(@at (Center, Center, Center) ∂z(v)) 
# # wz = Field(@at (Center, Center, Center) ∂z(w))
# ### AGSP calculation
# # upx = Field(@at (Center, Center, Center) Average(u, dims=1))
# # vpx = Field(@at (Center, Center, Center) Average(v, dims=1))
# # wpx = Field(@at (Center, Center, Center) Average(w, dims=1))
# dudz = Field(@at (Center, Center, Center) ∂z(u))
# dvdz = Field(@at (Center, Center, Center) ∂z(v))
# dwdz = Field(@at (Center, Center, Center) ∂z(w))
# AGSPu = -1*(u*w)*(dudz) # AGSP contribution 
# AGSPv = -1*(v*w)*(dvdz)
# AGSPw = -1*(w*w)*(dwdz) 
# AGSP = AGSPu + AGSPv + AGSPw
# # zC = znodes(grid, Center())
# # @inline builder(z,h) = permutedims(heaviside(h.*ones(100,)-z).*ones(100,1,500),(3,2,1))
# # const hv = builder(zC,hu)
# ### wave shear production calculation
# # uheaviside = FunctionField(location(u), heavisideh,u.grid)
# # vheaviside = FunctionField(location(v), heavisideh,v.grid)
# # xF, zC = xnodes(grid, Face()), znodes(grid, Center())

# # XF = [xF[i] for i in 1:Nx, j in 1:Nz]
# # ZC = [zC[j] for i in 1:Nx, j in 1:Nz]
# # hvu(model) = heavisideu(XF,ZC,simulation.model.clock.time,p)
# WSPu = (u*w)*(u_pert(0,0,simulation.model.clock.time,p))#*hvu#*umask # AGSP contribution 
# WSPv = (v*w)*(v_pert(0,0,simulation.model.clock.time,p))#*vmask 
# WSP = WSPu + WSPv# + AGSPw

# GSP = -1*(v*w)*(γ*(θ * N²)/(f))#*vheaviside#*vmask # geostrophic shear production
# BFLUX = (w+u*θ)*b # flux from buoyancy
# dpudx = Field(@at (Center, Center, Center) ∂z(θ*pr*u))
# dpvdy = Field(@at (Center, Center, Center) ∂y(pr*v))
# dpwdz = Field(@at (Center, Center, Center) ∂z(wp))
# dkwdz = Field(@at (Center, Center, Center) ∂z(wk))
# PWORK= -1*dpwdz # work due to pressure
# KTRANS = -1*dkwdz
# dk2dz2 = Field(@at (Center, Center, Center) ∂z(∂z(k)))
# KDISS = ν*dk2dz2
# Ri = RichardsonNumber(model, add_background=true)
# Ro = RossbyNumber(model)

# output writers
# output = (; u, U, v, V, w, b, B) # , PV , dbdz, dBdz, ε , Ri, Ro
# output2 = (; GSP, BFLUX, WSP) #AGSP, k, E,
# output2 = (; k, E, AGSP, WSP, GSP, BFLUX)
# output = (; u, v, w)
# # output3 = (; k)

# simulation.output_writers[:fields] = NetCDFOutputWriter(model, output;
#                                                           schedule = TimeInterval(0.05*(2*pi)/f),
#                                                           filename = path_name*"BBL_w_O_flow.nc",
#                                                           overwrite_existing = true)

# simulation.output_writers[:diagnostics] = NetCDFOutputWriter(model, output2;
#                                                           schedule = TimeInterval(0.005*(2*pi)/f),
#                                                           filename = path_name*"BBL_w_O_TKE.nc",
#                                                           overwrite_existing = true)


# With initial conditions set and an output writer at the ready, we run the simulation
simulation.stop_time = 15.01*((2π)/f)seconds

run!(simulation,pickup=true)
