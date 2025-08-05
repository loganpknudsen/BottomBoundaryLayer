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
using Parameters
using ArgParse
using CUDA 
using Oceanostics

### Load in Parameters

function parse_commandline()
        s = ArgParseSettings()
        @add_arg_table s begin
            "paramset"
                help = "sets which parameters to use"
                default = "f1e4theta029N21e5delta05Vinf005gammau"
            "--fruntime", "-T"
                help = "how many inertial periods sim runs"
                arg_type = Float64
                default = 30.0
            "--path"
                help = "pathname to save data under"
                default = "/glade/derecho/scratch/knudsenl/data/new_data/"
            "--Sinf"
                help = "Sinf"
                default = 1.0
            "--gamma"
                help = "PV parameter"
                default = 0.5
            "--suffix"
                help = "parameter set name"
                default = "f1e4theta029N21e5delta05Vinf005gammau"
        end
        return parse_args(s)
end

args=parse_commandline()

@info("command line args:")
for (arg,val) in args
  @info("   $arg => $val")
end

### Path Saved To

path_name = "/glade/derecho/scratch/knudsenl/data/new_data/test/" # args["path"]
print(args["suffix"])
setname = args["suffix"]

### Load in Parameters

@info "Loading parameters..."
include("parameters_stblty.jl")
params = getproperty(DenseParams(), Symbol(setname))

print(params)

### grid specifications
arch = has_cuda_gpu() ? GPU() : CPU()
@info("Arch => $arch")

Lx = 2000meters
Lz = params.Lz
Nx = 1024 # 512 originally
Nz = params.Nz # # 128 originally Note to self, maintain 2 to 1 resolution ration

grid = RectilinearGrid(arch; topology = (Periodic, Flat, Bounded),
                       size = (Nx, Nz),
                       x = (0, Lx),
                       z = (0,Lz))


### tilted domain parameters
const θ = 0.01 # params.θ 
const f = 1e-4 # params.f
ĝ = [θ, 0, 1] # gravity vector

### realistic mid latitude for now
buoyancy = Buoyancy(model = BuoyancyTracer(), gravity_unit_vector = -ĝ)
coriolis = ConstantCartesianCoriolis(f = f, rotation_axis = ĝ)

### parameters for simulation
const V∞ = params.V∞ # m s⁻¹ interior velocity
const N² = params.N² # interior stratification
const S∞ = params.S # slope burger number
const fˢ = f*(1+S∞^2)^(0.5) # modified oscillation
const δ = params.δ # geostrophic scaling factor
const γ = params.γ  # stratification parameter
const Λ = N²*γ*θ/f
const H = V∞/Λ # Height of Boundary Layer
const uₒ = δ*Λ  # Initial shear perturbation
const ϕ = params.ϕ

# array of paramerers for background function

p =(; N², θ, f, V∞, H, γ, uₒ, fˢ, Λ, ϕ)

# heaviside function for boundary layer

heaviside(x,z) = 0.5*(1+tanh(10000*z))

# oscillation functions for background

@inline sn_fn(x,z,t,p) = sin(p.fˢ*t+p.ϕ)
@inline cs_fn(x,z,t,p) = cos(p.fˢ*t+p.ϕ)

u_pert(x,z,t,p) = p.uₒ*cs_fn(x,z,t,p) 
v_pert(x,z,t,p) = -f*p.uₒ/(p.fˢ)*sn_fn(x,z,t,p)
b_pert(x,z,t,p) = -p.N²*p.θ*p.uₒ/(p.fˢ)*sn_fn(x,z,t,p)

### Total Background Velocity and Buoyancy

u_adjustment(x, z, t, p) = u_pert(x,z,t,p)*(p.H-z)*heaviside(x,p.H-z)
v_adjustment(x, z, t, p) = p.V∞ - p.Λ*(p.H-z)*heaviside(x,p.H-z) + v_pert(x,z,t,p)*(p.H-z)*heaviside(x,p.H-z)
constant_stratification(x, z, t, p) = p.N²*(x*p.θ + z) + p.N²*p.γ*(p.H-z)*heaviside(x,p.H-z) + b_pert(x,z,t,p)*(p.H-z)*heaviside(x,p.H-z)

U_field = BackgroundField(u_adjustment, parameters=p)
V_field = BackgroundField(v_adjustment, parameters=p)
B_field = BackgroundField(constant_stratification, parameters=p)

### Boundary Conditions for Buoyancy

b_bc_top= GradientBoundaryCondition(-1*N²)
b_bc_bottom= ValueBoundaryCondition(0) 

buoyancy_grad = FieldBoundaryConditions(top = b_bc_top, bottom=b_bc_bottom) 

### diffusitivity and viscosity values for closure

const ν1 = 1e-5
closure = ScalarDiffusivity(ν=ν1, κ=ν1)

start_time = time_ns()

### model set up 

model = NonhydrostaticModel(; grid, buoyancy, coriolis, closure,
                            timestepper = :RungeKutta3,
                            advection =  Centered(order=2), # Advection 
                            tracers = :b,
                            boundary_conditions = (; b=buoyancy_grad),
                            background_fields = (; u=U_field, v=V_field, b=B_field))

### initial conditions to start instability

ns = V∞*10^(-6) # standard deviation for noise

ui(x, z) = ns*Random.randn()
vi(x, z) = ns*Random.randn()
wi(x, z) = ns*Random.randn()

### set simulation and decide run time

set!(model, u=ui, v=vi, w=wi)

simulation = Simulation(model, Δt = 1seconds, stop_time = params.T*((2*pi)/fˢ)seconds) 

### time step wizard

wizard = TimeStepWizard(cfl=0.95, max_change=1.1seconds, max_Δt=100.0seconds, min_Δt=0.01seconds) 
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5)) 

### progress message to check time step size

progress_message(sim) =
        @printf("i: %04d, t: %s, Δt: %s, wall time: %s, Intertial Period %s\n",
        sim.model.clock.iteration, prettytime(sim.model.clock.time),
        prettytime(sim.Δt), prettytime((time_ns() - start_time) * 1e-9),sim.model.clock.time*fˢ/2π)

simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(2π/fˢ) ) 

### diagnostic calculations, it is saved in 2 files with one saving the flow field and the other tke diagnostics
### calculate the pertubation in velocities

ua, va, w = model.velocities # unaveraged velocities
um = Field(Average(ua, dims=(1))) # averging functions
vm = Field(Average(va, dims=(1)))
u = Field(ua - um) # calculating the Pertubations
v = Field(va - vm)

### background velocity and buoyancy

ub = model.background_fields.velocities.u
vb = model.background_fields.velocities.v
B = model.background_fields.tracers.b

### Buoyancy Pertubation calculation

ba = model.tracers.b
bm = Field(Average(ba, dims=(1)))
b = Field(ba - bm)

### Total Richardson Number calculation
Ri = RichardsonNumber(model, ub+ua, vb+va, w, B+ba)

### Potential Vorticity calculation

PV = ErtelPotentialVorticity(model, ub+ua, vb+va, w, B+ba, coriolis) 

### Dissaption calcuation

eps = KineticEnergyDissipationRate(model; U = um, V = vm, W = 0)
E = Field(Average(eps)) # kinetic energy dissaption calcualtion

### TKE caluclation

k_c = Oceanostics.TurbulentKineticEnergy(model, u, v, w)
k = Field(Average(k_c)) # TKE calculation

### AGSP calculation

AGSP_c =Oceanostics.ZShearProductionRate(model, u, v, w, um, vm, 0)
AGSP = Field(Average(AGSP_c))

### WSP calculation

@inline sn_fn(x,z,t,p) = sin(p.fˢ*t+p.ϕ)
@inline cs_fn(x,z,t,p) = cos(p.fˢ*t+p.ϕ)

upert(x,z,t,p) =  p.uₒ*cs_fn(x,z,t,p) *(p.H-z)*heaviside(x,p.H-z)# shear
vpert(x,z,t,p) = -f*p.uₒ/(p.fˢ)*sn_fn(x,z,t,p)*(p.H-z)*heaviside(x,p.H-z)

UPERT = Oceananigans.Fields.FunctionField{Center, Center, Center}(upert, grid, clock= model.clock, parameters = p)
VPERT = Oceananigans.Fields.FunctionField{Center, Center, Center}(vpert, grid, clock= model.clock, parameters = p)

WSP_c = Oceanostics.ZShearProductionRate(model, u, v, w, UPERT, VPERT, 0)
WSP = Field(Average(WSP_c))

### Tangent calculation

@inline tnd_fn(x,z,t,p) = tand(p.θ)

### GSP calcualtion

gshear(x,z,t,p) = p.V∞-p.Λ*(p.H-z)*heaviside(x,p.H-z)
GSHEAR = Oceananigans.Fields.FunctionField{Center, Center, Center}(gshear, grid, clock= model.clock, parameters = p)
GSP_c = Oceanostics.ZShearProductionRate(model, u, v, w, 0, GSHEAR, 0)
GSP = Field(Average(GSP_c))

### BP calcuation

BFLUX_c = Oceanostics.BuoyancyProductionTerm(model; velocities=(u=u, v=v, w=w), tracers=(b=b,))
BFLUX =  Field(Average(BFLUX_c))

### Output Writers array

output = (; u, ua, ub, v, va, vb, w, b, ba, B, PV) # pertubation fields and PV
output2 = (; k, E, GSP, WSP, AGSP, BFLUX) # TKE Diagnostic Calculations 

### Output Writers

simulation.output_writers[:fields] = NetCDFOutputWriter(model, output;
                                                          schedule = TimeInterval(0.05*(2*pi)/fˢ),
                                                          filename = path_name*"flow_fields_height_"*string(hu)*"_interior_velocity_"*string(V∞)*"_visc_"*string(ν1)*"_Sinf_"*string(S∞)*"_gamma_"*string(γ)*"_theta_"*string(θ)*"_f_"*string(f)*"_N2_"*string(N²)*".nc",
                                                          overwrite_existing = true)

simulation.output_writers[:diagnostics] = NetCDFOutputWriter(model, output2;
                                                          schedule = TimeInterval(0.005*(2*pi)/fˢ),
                                                          filename = path_name*"TKE_terms_height_"*string(hu)*"_interior_velocity_"*string(V∞)*"_visc_"*string(ν1)*"_Sinf_"*string(S∞)*"_gamma_"*string(γ)*"_theta_"*string(θ)*"_f_"*string(f)*"_N2_"*string(N²)*".nc",
                                                          overwrite_existing = true)

### Run Simulation

run!(simulation) 
