### Load in Packages
using Pkg
Pkg.activate(".")
Pkg.instantiate()
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

path_name = "/glade/derecho/scratch/knudsenl/data/" #args["path"]


# grid specifications
arch = has_cuda_gpu() ? GPU() : CPU()
@info("Arch => $arch")

Lx = 2000meters
Lz = 200meters
Nx = 512 # 512 originally
Nz = 128 # # 128 originally Note to self, maintain 2 to 1 resolution ration

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
const f = 1e-4 # coriolis parameter
const N² = 1e-5 # interior stratification
const S∞ = (N²*θ^2)/(f^2) # slope burger number
const γ = (1+S∞)^(-1) # 0 PV parameter
const hu = ceil((f*V∞)/(γ*N²*θ)) # Height of Boundary Layer
const fˢ=(f^2+θ^2*N²)^(0.5) # modified oscillation
const δ = 0.5
const vₒ = γ*(N²*θ)/(f)*δ # Initial v shear perturbation
# a1-h1 are constants for the following oscillations, calculate here for efficiency
const a1 = (f*vₒ)/(fˢ) 
const b1 = (f^2*vₒ)/(fˢ)^2
const d1 = ((fˢ^2-f^2)*vₒ)/(fˢ)^2
const e1 = N²*θ*(f*vₒ)/(fˢ)^2

# array of paramerers for background function
p =(; N², θ, f, V∞, hu, γ, vₒ, fˢ, a1, b1, d1, e1)

# heaviside function for boundary layer
@inline heaviside(x,z) = ifelse(z < 0, zero(z), one(z))

# oscillation functions for background
@inline sn_fn(x,z,t,p) = sin(p.fˢ*t)
@inline cs_fn(x,z,t,p) = cos(p.fˢ*t)

u_pert(x,z,t,p) = p.a1*sn_fn(x,z,t,p) # shear
v_pert(x,z,t,p) = p.b1*cs_fn(x,z,t,p) + p.d1
b_pert(x,z,t,p) = p.e1*(cs_fn(x,z,t,p) - 1)

u_adjustment(x, z, t, p) = u_pert(x,z,t,p)*(p.hu-z)*heaviside(x,p.hu-z)
v_adjustment(x, z, t, p) = p.V∞ - p.γ*(p.θ * p.N²)/(p.f)*(p.hu-z)*heaviside(x,p.hu-z) + v_pert(x,z,t,p)*(p.hu-z)*heaviside(x,p.hu-z)
constant_stratification(x, z, t, p) = p.N²*x*p.θ + p.N²*z + p.N²*p.γ*(p.hu-z)*heaviside(x,p.hu-z) + b_pert(x,z,t,p)*(p.hu-z)*heaviside(x,p.hu-z)

U_field = BackgroundField(u_adjustment, parameters=p)
V_field = BackgroundField(v_adjustment, parameters=p)
B_field = BackgroundField(constant_stratification, parameters=p)

# Boundary condition set up
# Free Slip Boundary Conditions

b_bc_top= FluxBoundaryCondition(-1*N²)
# b_bc_bottom= GradientBoundaryCondition(N²)
buoyancy_grad = FieldBoundaryConditions(top=b_bc_top) # ,bottom=b_bc_bottom

# diffusitivity and viscosity values for closure
const ν1 = 1e-4
closure = ScalarDiffusivity(ν=ν1, κ=ν1)

start_time = time_ns()

# model set up 
model = NonhydrostaticModel(; grid, buoyancy, coriolis, closure,
                            timestepper = :RungeKutta3,
                            advection = WENO(),
                            tracers = :b,
                            boundary_conditions = (; b=buoyancy_grad),
                            background_fields = (; u=U_field, v=V_field, b=B_field))

simulation = Simulation(model, Δt = 1seconds, stop_time = 1.0*((2*pi)/f)seconds, verbose=false)

wizard = TimeStepWizard(cfl=0.95, max_change=1.1seconds, max_Δt=100.0seconds, min_Δt=0.01seconds) 
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10)) 

function grow_instability!(simulation, energy)
    # Initialize
    simulation.model.clock.iteration = 0
    t₀ = simulation.model.clock.time = 0
    compute!(energy)
    energy₀ = CUDA.@allowscalar energy[1,1,1]

    # Grow
    run!(simulation)

    # Analyze
    compute!(energy)
    energy₁ =  CUDA.@allowscalar energy[1,1,1]
    Δτ = simulation.model.clock.time - t₀

    # ½(u² + v²) ~ exp(2 σ Δτ)
    σ = growth_rate = log(energy₁ / energy₀) / 2Δτ

    return growth_rate
end

"""
    rescale!(model, energy; target_kinetic_energy = 1e-3)

Rescales all model fields so that `energy = target_kinetic_energy`.
"""
function rescale!(model, energy; target_kinetic_energy = 1e-6)
    compute!(energy)
    rescale_factor = √(target_kinetic_energy / energy)

    for f in merge(model.velocities, model.tracers)
        f .*= rescale_factor
    end

    return nothing
end

using Printf

@inline convergence(σ) = length(σ) > 1 ? abs((σ[end] - σ[end-1]) / σ[end]) : 9.1e18 # pretty big (not Inf tho)

"""
    estimate_growth_rate(simulation, energy, ω; convergence_criterion=1e-3)

Estimates the growth rate iteratively until the relative change
in the estimated growth rate ``σ`` falls below `convergence_criterion`.

Returns ``σ``.
"""
function estimate_growth_rate(simulation, energy, convergence_criterion=1e-3)
    σ = Vector()
    power_method_data = Vector()
    push!(power_method_data, (deepcopy(σ)))

    while convergence(σ) > convergence_criterion
        compute!(energy)
        
        es =  CUDA.@allowscalar energy[1,1,1]
        @info @printf("About to start power method iteration %d; kinetic energy: %.2e\n", length(σ)+1,es) # , energy
        push!(σ, grow_instability!(simulation, energy))
        compute!(energy)
        es =  CUDA.@allowscalar energy[1,1,1]
        @info @printf("Power method iteration %d, kinetic energy: %.2e, σⁿ: %.2e, relative Δσ: %.2e\n",
                       length(σ), es, σ[end], convergence(σ))

        rescale!(simulation.model, energy)
        push!(power_method_data, (deepcopy(σ)))
    end

    return σ, power_method_data
end

ua, va, wa = model.velocities # change back to ua, va, wa
um = Field(@at (Face, Center, Center) Average(ua, dims=1)) #averaging
vm = Field(@at (Center, Face, Center) Average(va, dims=1))
wm = Field(@at (Center, Center, Face) Average(wa, dims=1))
u = Field(@at (Center, Center, Center) ua - um) # calculating the Pertubations
v = Field(@at (Center, Center, Center) va - vm)
w = Field(@at (Center, Center, Center) wa - wm)

mean_perturbation_kinetic_energy = Field(Average(Oceanostics.TurbulentKineticEnergy(model, u, v, w))) # TKE calculation

ns = 10^(-4) # standard deviation for noise

# initial conditions to start instability
ui(x, z) = ns*Random.randn()
vi(x, z) = ns*Random.randn()
wi(x, z) = ns*Random.randn()

# set simulation and decide run time
set!(model, u=ui, v=vi, w=wi)

rescale!(simulation.model, mean_perturbation_kinetic_energy, target_kinetic_energy=1e-6)
growth_rates, power_method_data = estimate_growth_rate(simulation, mean_perturbation_kinetic_energy)
@info "Power iterations converged! Estimated growth rate: $(growth_rates[end])"
