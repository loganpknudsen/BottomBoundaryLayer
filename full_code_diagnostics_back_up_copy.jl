### Load in Packages
using Pkg
Pkg.activate(".")
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
using CUDA 
using Oceanostics

# Path file is saved under
# function parse_commandline()
#     s = ArgParseSettings()
#     @add_arg_table s begin
#         "path"
#         help = "pathname to save data under"
#         default = ""
#     end
#     return parse_args(s)
# end

# args=parse_commandline()

# @info("command line args:")
# for (arg,val) in args
#   @info("   $arg => $val")
# end

path_name = "/glade/derecho/scratch/knudsenl/data/new_data/paper_data/" #args["path"]


# grid specifications
arch = has_cuda_gpu() ? GPU() : CPU()
@info("Arch => $arch")

Lx = 6000meters
Lz = 1100meters
Nx = 3072 #1024 # 512 originally
Nz = 1408 #256 #128 # # 128 originally Note to self, maintain 2 to 1 resolution ration

grid = RectilinearGrid(arch; topology = (Periodic, Flat, Bounded),
                       size = (Nx, Nz),
                       x = (0, Lx),
                       z = (0,Lz))


# tilted domain parameters
const θ = 0.29 # degrees 10^(-2) is previous value for 110 meter layer
const f = 1e-4
ĝ = [sind(θ), 0, cosd(θ)] # gravity vector

# realistic mid latitude for now
buoyancy = Buoyancy(model = BuoyancyTracer(), gravity_unit_vector = -ĝ)
coriolis = ConstantCartesianCoriolis(f = f, rotation_axis = ĝ)

# parameters for simulation
const V∞ = 0.01 # m s⁻¹ interior velocity
const f = 1e-4 # coriolis parameter
const N² = 2e-7 # interior stratification
const S∞ = (N²*tand(θ)^2)/(f^2) # slope burger number
const fˢ = cosd(θ)*(f^2+tand(θ)^2*N²)^(0.5) # modified oscillation
const δ = 0.5
const γ = ((cosd(θ)*(1+S∞*(1-δ)))^(-1)+(3-S∞)*((3*cosd(θ)*(1+S∞)-4*δ*cosd(θ)*S∞))^(-1))/2 # 0 PV parameter
const hu = (f*V∞)/(γ*N²*tand(θ)) # Height of Boundary Layer
const uₒ = 0 # Initial u shear perturbation
const vₒ = δ*γ*(N²*tand(θ))/(f) # Initial v shear perturbation
const bₒ = 0 #-γ*((N²*θ)/(f))^2*δ#vₒ*((θ*N²)/(f))*0.1 # initial stratification perturbation
# a1-h1 are constants for the following oscillations, calculate here for efficiency
const a1 = (f*cosd(θ)*vₒ)/(fˢ) 
const b1 = (f^2*cosd(θ)^2*vₒ)/(fˢ)^2
# const c1 = (f*uₒ)/(fˢ)
# const d1 = ((fˢ^2-f^2)*vₒ-f*bₒ*θ)/(fˢ)^2
const e1 = N²*sind(θ)*f*cosd(θ)*vₒ/(fˢ)^2
# const h1 = (N²*θ*uₒ)/(fˢ)

# array of paramerers for background function
p =(; N², θ, f, V∞, hu, γ, uₒ, vₒ, bₒ, fˢ, a1, b1, e1)

# heaviside function for boundary layer
heaviside(x,z) = 0.5*(1+tanh(10000*z))

# oscillation functions for background
@inline sn_fn(x,z,t,p) = sin(p.fˢ*t)
@inline cs_fn(x,z,t,p) = cos(p.fˢ*t)

u_pert(x,z,t,p) = p.a1*sn_fn(x,z,t,p) # shear
v_pert(x,z,t,p) = p.vₒ+p.b1*(cs_fn(x,z,t,p)-1)
b_pert(x,z,t,p) = p.e1*(cs_fn(x,z,t,p) - 1)

u_adjustment(x, z, t, p) = u_pert(x,z,t,p)*(p.hu-z)*heaviside(x,p.hu-z)
v_adjustment(x, z, t, p) = p.V∞ - p.γ*(tand(p.θ) * p.N²)/(p.f)*(p.hu-z)*heaviside(x,p.hu-z) + v_pert(x,z,t,p)*(p.hu-z)*heaviside(x,p.hu-z)
constant_stratification(x, z, t, p) = p.N²*x*sind(p.θ) + p.N²*z*cosd(p.θ) + p.N²*p.γ*(p.hu-z)*heaviside(x,p.hu-z) + b_pert(x,z,t,p)*(p.hu-z)*heaviside(x,p.hu-z)

U_field = BackgroundField(u_adjustment, parameters=p)
V_field = BackgroundField(v_adjustment, parameters=p)
B_field = BackgroundField(constant_stratification, parameters=p)

b_bc_top= GradientBoundaryCondition(-1*N²*cosd(θ))

buoyancy_grad = FieldBoundaryConditions(top = b_bc_top) 

# diffusitivity and viscosity values for closure
const ν1 = 1e-6
closure = ScalarDiffusivity(ν=ν1, κ=ν1)

start_time = time_ns()

# model set up 
model = NonhydrostaticModel(; grid, buoyancy, coriolis, closure,
                            timestepper = :RungeKutta3,
                            advection =  Centered(order=2), # WENO(), #
                            tracers = :b,
                            boundary_conditions = (; b=buoyancy_grad),
                            background_fields = (; u=U_field, v=V_field, b=B_field))

ns = 10^(-10) # standard deviation for noise

# initial conditions to start instability
ui(x, z) = ns*Random.randn() #*heaviside(x,hu-1-z)
vi(x, z) = ns*Random.randn() #*heaviside(x,hu-1-z)
wi(x, z) = ns*Random.randn() #*heaviside(x,hu-1-z)
# bp(x,z) = ns*Random.randn()

# set simulation and decide run time
set!(model, u=ui, v=vi) #, w=wi)

simulation = Simulation(model, Δt = 1seconds, stop_time = 200.1*((2*pi)/fˢ)seconds) # stop_iteration=10

# time step wizard
wizard = TimeStepWizard(cfl=0.75, max_change=1.1seconds, max_Δt=100.0seconds, min_Δt=0.01seconds) 
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5)) 

# simulation.output_writers[:checkpointer] = Checkpointer(model; schedule=TimeInterval((5*(2*pi)/f)seconds), prefix="model_checkpoint")

# progress message to check time step size
progress_message(sim) =
        @printf("i: %04d, t: %s, Δt: %s, wall time: %s, Intertial Period %s\n",
        sim.model.clock.iteration, prettytime(sim.model.clock.time),
        prettytime(sim.Δt), prettytime((time_ns() - start_time) * 1e-9),sim.model.clock.time*f/2π)

simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(1000)) # TimeInterval(0.5*(2*pi)/f) 

# diagnostic calculations, it is saved in 2 files with one saving the flow field and the other tke diagnostics
# calculate the pertubation in velocities

ua, va, w = model.velocities
um = Field(Average(ua, dims=(1))) #averaging
vm = Field(Average(va, dims=(1)))
# wm = Field(Average(wa, dims=(1)))
u = Field(ua - um) # calculating the Pertubations
v = Field(va - vm)
# w = Field(wa - wm)
ub = model.background_fields.velocities.u
vb = model.background_fields.velocities.v
B = model.background_fields.tracers.b

# buoyancy pertubation calculation
ba = model.tracers.b
bm = Field(Average(ba, dims=(1)))
b = Field(ba - bm)
# bt = Field(@at (Center, Center, Center) B+ba)

# Ri = RichardsonNumber(model, ut, vt, wa, bt)
# Ro = RossbyNumber(model, ut, vt, wa, coriolis)
PV = ErtelPotentialVorticity(model, ub+ua, vb+va, w, B, coriolis) # potential vorticity calculation
eps = KineticEnergyDissipationRate(model; U = um, V = vm, W = 0)
E = Field(Average(eps)) # kinetic energy dissaption calcualtion
k_c = Oceanostics.TurbulentKineticEnergy(model, u, v, w)
k = Field(Average(k_c)) # TKE calculation

### AGSP calculation
AGSP_c =Oceanostics.ZShearProductionRate(model, u, v, w, um, vm, 0)
AGSP = Field(Average(AGSP_c))

### wave shear production calculation
@inline sn_fn(x,z,t,p) = sin(p.fˢ*t)
@inline cs_fn(x,z,t,p) = cos(p.fˢ*t)

upert(x,z,t,p) = p.a1*sn_fn(x,z,t,p)*(p.hu-z)*heaviside(x,p.hu-z)# shear
vpert(x,z,t,p) = (p.vₒ+p.b1*(cs_fn(x,z,t,p)-1))*(p.hu-z)*heaviside(x,p.hu-z)
bpert(x,z,t,p) = p.e1*(cs_fn(x,z,t,p) - 1)*(p.hu-z)*heaviside(x,p.hu-z)

UPERT = Oceananigans.Fields.FunctionField{Center, Center, Center}(upert, grid, clock= model.clock, parameters = p)
VPERT = Oceananigans.Fields.FunctionField{Center, Center, Center}(vpert, grid, clock= model.clock, parameters = p)

WSP_c = Oceanostics.ZShearProductionRate(model, u, v, w, UPERT, VPERT, 0)
WSP = Field(Average(WSP_c))

@inline tnd_fn(x,z,t,p) = tand(p.θ)

### geostrophic shear production calcualtion
gshear(x,z,t,p) = p.V∞-((p.γ * tnd_fn(x,z,t,p) * p.N²)/(p.f))*(p.hu-z)*heaviside(x,p.hu-z)
GSHEAR = Oceananigans.Fields.FunctionField{Center, Center, Center}(gshear, grid, clock= model.clock, parameters = p)
GSP_c = Oceanostics.ZShearProductionRate(model, u, v, w, 0, GSHEAR, 0)
GSP = Field(Average(GSP_c))

### Buoyancy Flux Calcuation
BFLUX_c = Oceanostics.BuoyancyProductionTerm(model; velocities=(u=u, v=v, w=w), tracers=(b=b,))
BFLUX =  Field(Average(BFLUX_c))

# output writers
output = (; u, ua, ub, v, va, vb, w, b, ba, B, PV) # pertubation fields and PV
output2 = (; k, E, GSP, WSP, AGSP, BFLUX) # TKE Diagnostic Calculations 

simulation.output_writers[:fields] = NetCDFOutputWriter(model, output;
                                                          schedule = TimeInterval(0.05*(2*pi)/fˢ),
                                                          filename = path_name*"flow_fields_height_"*string(hu)*"_theta_"*string(θ)*"_stratification_"*string(N²)*"_interior_velocity_"*string(V∞)*"_visc_"*string(ν1)*"_Sinf_"*string(S∞)*"_gamma_"*string(γ)*"_f_"*string(f)*".nc",
                                                          overwrite_existing = true)

simulation.output_writers[:diagnostics] = NetCDFOutputWriter(model, output2;
                                                          schedule = TimeInterval(0.005*(2*pi)/fˢ),
                                                          filename = path_name*"TKE_terms_height_"*string(hu)*"_theta_"*string(θ)*"_stratification_"*string(N²)*"_interior_velocity_"*string(V∞)*"_visc_"*string(ν1)*"_Sinf_"*string(S∞)*"_gamma_"*string(γ)*"_f_"*string(f)*".nc",
                                                          overwrite_existing = true)

# With initial conditions set and an output writer at the ready, we run the simulation
# simulation.stop_time = 15*((2π)/f)seconds

run!(simulation) # , pickup=true
