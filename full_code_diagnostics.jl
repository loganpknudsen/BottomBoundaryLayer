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

path_name = "/glade/derecho/scratch/knudsenl/data" #args["path"]


# grid specifications
arch = has_cuda_gpu() ? GPU() : CPU()
@info("Arch => $arch")

Lx = 2000meters
Lz = 200meters
Nx = 512
Nz = 128 # Note to self, maintain 2 to 1 resolution ration

grid = RectilinearGrid(arch; topology = (Periodic, Flat, Bounded),
                       size = (Nx, Nz),
                       x = (0, Lx),
                       z = (0,Lz))


# tilted domain parameters
θ = 1*10^(-3) # degrees 10^(-2) is previous value for 110 meter layer
ĝ = [sind(θ), 0, cosd(θ)] # gravity vector

# realistic mid latitude for now
buoyancy = Buoyancy(model = BuoyancyTracer(), gravity_unit_vector = -ĝ)
coriolis = ConstantCartesianCoriolis(f = 1e-4, rotation_axis = ĝ)

# parameters for simulation
const V∞ = 0.05 # m s⁻¹ interior velocity
const f = 1e-4 # coriolis parameter
const N² = 5e-5 # interior stratification
const S∞ = (N²*θ^2)/(f^2) # slope burger number
const γ = (1+S∞)^(-1) # 0 PV parameter
const hu = ceil((f*V∞)/(γ*N²*θ)) # Height of Boundary Layer
const fˢ=(f^2+θ^2*N²)^(0.5) # modified oscillation
const uₒ = 0 # Initial u shear perturbation
const vₒ = γ*(N²*θ)/(f)*0.5 # Initial v shear perturbation
const bₒ = vₒ*((θ*N²)/(f))*0.1 # initial stratification perturbation
# a1-h1 are constants for the following oscillations, calculate here for efficiency
const a1 = (f*vₒ+bₒ*θ)/(fˢ) 
const b1 = (f^2*vₒ+f*bₒ*θ)/(fˢ)^2
const c1 = (f*uₒ)/(fˢ)
const d1 = ((fˢ^2-f^2)*vₒ-f*bₒ*θ)/(fˢ)^2
const e1 = N²*θ*(f*vₒ+bₒ*θ)/(fˢ)^2
const h1 = (N²*θ*uₒ)/(fˢ)

# array of paramerers for background function
p =(; N², θ, f, V∞, hu, γ, uₒ, vₒ, bₒ, fˢ, a1, b1, c1, d1, e1, h1)

# heaviside function for boundary layer
@inline heaviside(x,z) = ifelse(z < 0, zero(z), one(z))

# oscillation functions for background
@inline sn_fn(x,z,t,p) = sin(p.fˢ*t)
@inline cs_fn(x,z,t,p) = cos(p.fˢ*t)

u_pert(x,z,t,p) = p.uₒ*cs_fn(x,z,t,p) +p.a1*sn_fn(x,z,t,p) # shear
v_pert(x,z,t,p) = p.b1*cs_fn(x,z,t,p) - p.c1*sn_fn(x,z,t,p)+p.d1
b_pert(x,z,t,p) = p.e1*cs_fn(x,z,t,p) - p.h1*sn_fn(x,z,t,p)+p.bₒ-p.e1

u_adjustment(x, z, t, p) = u_pert(x,z,t,p)*(p.hu-z)*heaviside(x,p.hu-z)
v_adjustment(x, z, t, p) = p.V∞ - p.γ*(p.θ * p.N²)/(p.f)*(p.hu-z)*heaviside(x,p.hu-z) + v_pert(x,z,t,p)*(p.hu-z)*heaviside(x,p.hu-z)
constant_stratification(x, z, t, p) = p.N²*x*p.θ + p.N²*z + p.N²*p.γ*(p.hu-z)*heaviside(x,p.hu-z) + b_pert(x,z,t,p)*(p.hu-z)*heaviside(x,p.hu-z)

U_field = BackgroundField(u_adjustment, parameters=p)
V_field = BackgroundField(v_adjustment, parameters=p)
B_field = BackgroundField(constant_stratification, parameters=p)

# Boundary condition set up
# Free Slip Boundary Conditions

b_bc_top= GradientBoundaryCondition(0)
b_bc_bottom= GradientBoundaryCondition(-1*N²)
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
# bp(x,z) = ns*Random.randn()

# set simulation and decide run time
set!(model, u=ui, v=vi, w=wi)

simulation = Simulation(model, Δt = 1seconds, stop_time = 60.01*((2*pi)/f)seconds) # stop_iteration=10

# time step wizard
wizard = TimeStepWizard(cfl=0.95, max_change=1.1seconds, max_Δt=100.0seconds, min_Δt=0.01seconds) 
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
ua, va, wa = model.velocities # change back to ua, va, wa
um = Field(@at (Face, Center, Center) Average(ua, dims=1)) #averaging
vm = Field(@at (Center, Face, Center) Average(va, dims=1))
wm = Field(@at (Center, Center, Face) Average(wa, dims=1))
u = Field(@at (Center, Center, Center) ua - um) # calculating the Pertubations
v = Field(@at (Center, Center, Center) va - vm)
w = Field(@at (Center, Center, Center) wa - wm)
ub = model.background_fields.velocities.u
vb = model.background_fields.velocities.u
B = model.background_fields.tracers.b

# buoyancy pertubation calculation
ba = model.tracers.b
bm = Field(@at (Center, Center, Center) Average(ba, dims=1))
b = Field(@at (Center, Center, Center) ba - bm)

PV = ErtelPotentialVorticity(model, add_background=true) # potential vorticity calculation
E = KineticEnergyDissipationRate(model; U = um, V = vm, W = wm) # kinetic energy dissaption calcualtion
k = Oceanostics.TurbulentKineticEnergy(model, u, v, w) # TKE calculation

### AGSP calculation
AGSP = Oceanostics.ZShearProductionRate(model, u, v, w, um, vm, wm)

### wave shear production calculation
upert(x,z,t,p) = (p.uₒ*cs_fn(x,z,t,p) + p.a1*sn_fn(x,z,t,p))*(p.hu-z)*heaviside(x,p.hu-z)
vpert(x,z,t,p) = (p.b1*cs_fn(x,z,t,p) - p.c1*sn_fn(x,z,t,p)+p.d1)*(p.hu-z)*heaviside(x,p.hu-z)

UPERT = Oceananigans.Fields.FunctionField{Center, Center, Center}(upert, grid, clock= model.clock, parameters = p)
VPERT = Oceananigans.Fields.FunctionField{Center, Center, Center}(vpert, grid, clock= model.clock, parameters = p)

WSP = Oceanostics.ZShearProductionRate(model, u, v, w, UPERT, VPERT, 0)

### geostrophic shear production calcualtion
gshear(x,z,t,p) = p.V∞-((p.γ * p.θ * p.N²)/(p.f))*(p.hu-z)*heaviside(x,p.hu-z)
GSHEAR = Oceananigans.Fields.FunctionField{Center, Center, Center}(gshear, grid, clock= model.clock, parameters = p)
GSP = Oceanostics.ZShearProductionRate(model, u, v, w, 0, GSHEAR, 0)

### Buoyancy Flux Calcuation
BFLUX =  Oceanostics.BuoyancyProductionTerm(model; velocities=(u=u, v=v, w=w), tracers=(b=b,))

# output writers
output = (; u, ua, ub, v, va, vb, w, wa, b, ba, B, PV) # pertubation fields and PV
output2 = (; k, E, GSP, WSP, AGSP, BFLUX) # TKE Diagnostic Calculations

simulation.output_writers[:fields] = NetCDFOutputWriter(model, output;
                                                          schedule = TimeInterval(0.05*(2*pi)/f),
                                                          filename = path_name*"flow_fields_height_"*string(hu)*"_theta_"*string(θ)*"_stratification_"*string(N²)*"_interior_velocity_"*string(V∞)*".nc",
                                                          overwrite_existing = true)

simulation.output_writers[:diagnostics] = NetCDFOutputWriter(model, output2;
                                                          schedule = TimeInterval(0.005*(2*pi)/f),
                                                          filename = path_name*"TKE_terms_height_"*string(hu)*"_theta_"*string(θ)*"_stratification_"*string(N²)*"_interior_velocity_"*string(V∞)*".nc",
                                                          overwrite_existing = true)

# With initial conditions set and an output writer at the ready, we run the simulation
# simulation.stop_time = 15*((2π)/f)seconds

run!(simulation) # , pickup=true
