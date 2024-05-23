using Oceananigans
using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.Grids: Center, Face
using Oceananigans.Units
using Random
using Printf
using ArgParse
using CUDA: has_cuda_gpu
using Oceanostics

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

# made grid correct shape, need to modify z boundaries to make sure they are no slip

# grid specifications
arch = has_cuda_gpu() ? GPU() : CPU()
@info("Arch => $arch")

Lx = 2000meters
Lz = 200meters
Nx = 500
Nz = 100 # Two to one ratio similar to previous would be Nx = 500 Nz = 100

grid = RectilinearGrid(arch; topology = (Periodic, Flat, Bounded),
                       size = (Nx, Nz),
                       x = (0, Lx),
                       z = (0,Lz))


# tilted domain parameters
θ = 10^(-2) # degrees 
# ĝ = [θ, 0, 1] # gravity vector small angle
ĝ = [sind(θ), 0, cosd(θ)] # gravity vector

# realustic mid latitude for now
buoyancy = Buoyancy(model = BuoyancyTracer(), gravity_unit_vector = -ĝ)
coriolis = ConstantCartesianCoriolis(f = 1e-4, rotation_axis = ĝ)

# parameters
const V∞ = 0.1 # m s⁻¹
const f = coriolis.fz
const N² = 1e-5 # interior stratification
#ϕ = 0
const S∞ = (N²*θ^2)/(f^2)
const γ = (1+S∞)^(-1) #(θ^2+1)*(1+S∞*(θ^2+1))^(-1)
const hu = ceil((f*V∞)/(γ*N²*θ)) # set to negative
const fˢ=(f^2+θ^2*N²)^(0.5)
const uₒ = 0#γ*(N²*θ)/(f)*cos(ϕ)
const vₒ = γ*(N²*θ)/(f)*0.5#*sin(ϕ)
const bₒ = vₒ*((θ*N²)/(f))*0.1 # initial stratification
const a1 = (f*vₒ+bₒ*θ)/(fˢ)
const b1 = (f^2*vₒ+f*bₒ*θ)/(fˢ)^2
const c1 = (f*uₒ)/(fˢ)
const d1 = ((fˢ^2-f^2)*vₒ-f*bₒ*θ)/(fˢ)^2
const e1 = N²*θ*(f*vₒ+bₒ*θ)/(fˢ)^2
const h1 = (N²*θ*uₒ)/(fˢ)


p =(; N², θ, f, V∞, hu, γ, uₒ, vₒ, bₒ, fˢ, a1, b1, c1, d1, e1, h1)

# background flow with geostrophic and ageostrophic shear 

# @inline interval(x,a,b) = ifelse(a<=x<=b, one(x), zero(x))

@inline heaviside(x) = ifelse(x < 0, zero(x), one(x))

@inline sn_fn(x,z,t,p) = sin(p.fˢ*t)
@inline cs_fn(x,z,t,p) = cos(p.fˢ*t)

u_pert(x,z,t,p) = p.uₒ*cs_fn(x,z,t,p) +p.a1*sn_fn(x,z,t,p)
v_pert(x,z,t,p) = p.b1*cs_fn(x,z,t,p) - p.c1*sn_fn(x,z,t,p)+p.d1
b_pert(x,z,t,p) = p.e1*cs_fn(x,z,t,p) - p.h1*sn_fn(x,z,t,p)+p.bₒ-p.e1

u_adjustment(x, z, t, p) =  u_pert(x,z,t,p)*(p.hu-z)*heaviside(p.hu-z)
v_adjustment(x, z, t, p) = p.V∞-p.γ*(p.θ * p.N²)/(p.f)*(p.hu-z)*heaviside(p.hu-z) + v_pert(x,z,t,p)*(p.hu-z)*heaviside(p.hu-z)
constant_stratification(x, z, t, p) = p.N²*x*p.θ + p.N²*z + p.N²*p.γ*(p.hu-z)*heaviside(p.hu-z) + b_pert(x,z,t,p)*(p.hu-z)*heaviside(p.hu-z)

U_field = BackgroundField(u_adjustment, parameters=p)
V_field = BackgroundField(v_adjustment, parameters=p)
B_field = BackgroundField(constant_stratification, parameters=p)

# Boundary condition set up
# Free Slip Boundary Conditions

b_bc_top= GradientBoundaryCondition(N²)
# b_bc_bottom= GradientBoundaryCondition(N²*(1-γ))
buoyancy_grad = FieldBoundaryConditions(top=b_bc_top) # ,bottom=b_bc_bottom

# boundary_conditions=(;b=buoyancy_grad),
closure = ScalarDiffusivity(ν=1e-4, κ=1e-4)

start_time = time_ns()

model = NonhydrostaticModel(; grid, buoyancy, coriolis, closure,
                            timestepper = :RungeKutta3,
                            advection = WENO(),
                            tracers = :b,
                            boundary_conditions = (; b=buoyancy_grad),
                            background_fields = (; u=U_field, v=V_field, b=B_field))

ns = 10^(-4) # standard deviation for noise

u₀(x, z) = ns*Random.randn()
v₀(x, z) = ns*Random.randn()
w₀(x, z) = ns*Random.randn()
# bₒ(x,y,z) = 0.005*Random.randn()

set!(model, u=u₀, v=v₀, w=w₀)

simulation = Simulation(model, Δt = 1, stop_time = 50*(2*pi)/f)


wizard = TimeStepWizard(cfl=0.7, max_change=1.1, max_Δt=10.0, min_Δt=0.01) 
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5)) 

progress_message(sim) =
        @printf("i: %04d, t: %s, Δt: %s, wall time: %s\n",
        sim.model.clock.iteration, prettytime(sim.model.clock.time),
        prettytime(sim.Δt), prettytime((time_ns() - start_time) * 1e-9))

simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(10000) ) # TimeInterval(0.5*(2*pi)/f) 

# and add an output writer that saves the vertical velocity field every two iterations:

u, v, w = model.velocities
ub = model.background_fields.velocities.u
vb = model.background_fields.velocities.v
b = model.tracers.b
B∞ = model.background_fields.tracers.b

U = u + ub
V = v + vb #+ V∞
B = b + B∞
# dBdz = Field(@at (Center, Center, Center) ∂z(b+B∞))

PV = ErtelPotentialVorticity(model, add_background=true)
KE = KineticEnergy(model)
E = KineticEnergyDissipationRate(model)
AGSPu = u*w*u_pert(x,z,t,p)
AGSPv = v*w*v_pert(x,z,t,p)
AGSP = AGSPu + AGSPv
GSP = - v*w*γ*(θ * N²)/(f)
BFLUX = w*b
k = 0.5*(u^2+v^2+w^2) # pertubation kinetic energy
# Ri = RichardsonNumber(model, add_background=true)
# Ro = RossbyNumber(model)


output = (; u, U, v, V, w, b, B, PV, KE, E, AGSP, GSP, BFLUX, k) # , ε , Ri, Ro

# u,v,w = model.velocities

# output = (;u,v,w,model.tracers.b,U=(model.background_fields.velocities.u+0*u),V=(model.background_fields.velocities.v+0*v),B=(model.background_fields.tracers.b+0*model.tracers.b))

# ε = Field(KineticEnergyDissipationRate(model))
# dBdz = Field(@at (Center, Center, Center) ∂z(model.tracers.b+model.background_fields.tracers.b))
# u_m_flux = u*w
# v_m_flux = v*w

# output = merge(output, (; E=ε, N2=dBdz, UM=u_m_flux, VM=v_m_flux,))

simulation.output_writers[:fields] = NetCDFOutputWriter(model, output;
                                                          schedule = TimeInterval(0.05*(2*pi)/f),
                                                          filename = path_name*"BBL_w_O_50_base_test.nc",
                                                          overwrite_existing = true)

# With initial conditions set and an output writer at the ready, we run the simulation

run!(simulation)
