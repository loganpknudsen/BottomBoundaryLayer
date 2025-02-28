### Load in Packages
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Oceananigans
using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.AbstractOperations: @at, Average
using Oceananigans.Grids: Center, Face
using Oceananigans.Fields: FunctionField, interior
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

simulation = Simulation(model, Δt = 10seconds, stop_time = 1.01*((2*pi)/f)seconds, verbose=false)

function grow_instability!(simulation, energy)
    # Initialize
    simulation.model.clock.iteration = 0
    t₀ = simulation.model.clock.time = 0
    compute!(energy)
    energy₀ = energy[1, 1, 1]

    # Grow
    run!(simulation)

    # Analyze
    compute!(energy)
    energy₁ = energy[1, 1, 1]
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

convergence(σ) = length(σ) > 1 ? abs((σ[end] - σ[end-1]) / σ[end]) : 9.1e18 # pretty big (not Inf tho)

"""
    estimate_growth_rate(simulation, energy, ω; convergence_criterion=1e-3)

Estimates the growth rate iteratively until the relative change
in the estimated growth rate ``σ`` falls below `convergence_criterion`.

Returns ``σ``.
"""
function estimate_growth_rate(simulation, energy, ω, b; convergence_criterion=1e-3)
    σ = []
    power_method_data = []
    compute!(ω)
    push!(power_method_data, (ω=collect(interior(ω, :, 1, :)), b=collect(interior(b, :, 1, :)), σ=deepcopy(σ)))

    while convergence(σ) > convergence_criterion
        compute!(energy)

        @info @sprintf("About to start power method iteration %d; kinetic energy: %.2e", length(σ)+1, energy)
        push!(σ, grow_instability!(simulation, energy))
        compute!(energy)

        @info @sprintf("Power method iteration %d, kinetic energy: %.2e, σⁿ: %.2e, relative Δσ: %.2e",
                       length(σ), energy, σ[end], convergence(σ))

        compute!(ω)
        rescale!(simulation.model, energy)
        push!(power_method_data, (ω=collect(interior(ω, :, 1, :)), b=collect(interior(b, :, 1, :)), σ=deepcopy(σ)))
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

# buoyancy pertubation calculation
ba = model.tracers.b
bm = Field(@at (Center, Center, Center) Average(ba, dims=1))
b = Field(@at (Center, Center, Center) ba - bm)


mean_perturbation_kinetic_energy = Field(Average(Oceanostics.TurbulentKineticEnergy(model, u, v, w))) # TKE calculation

perturbation_vorticity = ErtelPotentialVorticity(model, u, v, w, b, coriolis) # potential vorticity calculation

xω, yω, zω = nodes(perturbation_vorticity)
xb, yb, zb = nodes(b)

ns = 10^(-6) # standard deviation for noise

# initial conditions to start instability
ui(x, z) = ns*Random.randn()
vi(x, z) = ns*Random.randn()
wi(x, z) = ns*Random.randn()
# bp(x,z) = ns*Random.randn()

# set simulation and decide run time
set!(model, u=ui, v=vi, w=wi)

rescale!(simulation.model, mean_perturbation_kinetic_energy, target_kinetic_energy=1e-6)
growth_rates, power_method_data = estimate_growth_rate(simulation, mean_perturbation_kinetic_energy, perturbation_vorticity, b)
@info "Power iterations converged! Estimated growth rate: $(growth_rates[end])"

n = Observable(1)

fig = Figure(size=(800, 600))

kwargs = (xlabel="x", ylabel="z", limits = ((xω[1], xω[end]), (zω[1], zω[end])), aspect=1,)

ω_title(t) = t === nothing ? @sprintf("vorticity") : @sprintf("vorticity at t = %.2f", t)
b_title(t) = t === nothing ? @sprintf("buoyancy")  : @sprintf("buoyancy at t = %.2f", t)

ax_ω = Axis(fig[2, 1]; title = ω_title(nothing), kwargs...)

ax_b = Axis(fig[2, 3]; title = b_title(nothing), kwargs...)

ωₙ = @lift power_method_data[$n].ω
bₙ = @lift power_method_data[$n].b

σₙ = @lift [(i-1, i==1 ? NaN : growth_rates[i-1]) for i in 1:$n]

ω_lims = @lift (-maximum(abs, power_method_data[$n].ω) - 1e-16, maximum(abs, power_method_data[$n].ω) + 1e-16)
b_lims = @lift (-maximum(abs, power_method_data[$n].b) - 1e-16, maximum(abs, power_method_data[$n].b) + 1e-16)

hm_ω = heatmap!(ax_ω, xω, zω, ωₙ; colorrange = ω_lims, colormap = :balance)
Colorbar(fig[2, 2], hm_ω)

hm_b = heatmap!(ax_b, xb, zb, bₙ; colorrange = b_lims, colormap = :balance)
Colorbar(fig[2, 4], hm_b)

eigentitle(σ, t) = length(σ) > 0 ? @sprintf("Iteration #%i; growth rate %.2e", length(σ), σ[end]) : @sprintf("Initial perturbation fields")
σ_title = @lift eigentitle(power_method_data[$n].σ, nothing)

ax_σ = Axis(fig[1, :];
            xlabel = "Power iteration",
            ylabel = "Growth rate",
            title = σ_title,
            xticks = 1:length(power_method_data)-1,
            limits = ((0.5, length(power_method_data)-0.5), (-0.25, 0.25)))

scatter!(ax_σ, σₙ; color = :blue)

frames = 1:length(power_method_data)

record(fig, "powermethod.mp4", frames, framerate=1) do i
       n[] = i
end
