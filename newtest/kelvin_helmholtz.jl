using Oceananigans

grid = RectilinearGrid(size=(64, 64), x=(-5, 5), z=(-5, 5),
                       topology=(Periodic, Flat, Bounded))
shear_flow(x, z, t) = tanh(z)

stratification(x, z, t, p) = p.h * p.Ri * tanh(z / p.h)

U = BackgroundField(shear_flow)

B = BackgroundField(stratification, parameters=(Ri=0.1, h=1/4))

model = NonhydrostaticModel(timestepper = :RungeKutta3,
                              advection = UpwindBiasedFifthOrder(),
                                   grid = grid,
                               coriolis = nothing,
                      background_fields = (u=U, b=B),
                                closure = ScalarDiffusivity(ν=2e-4, κ=2e-4),
                               buoyancy = BuoyancyTracer(),
                                tracers = :b)
simulation = Simulation(model, Δt=0.1, stop_iteration=150, verbose=false)

"""
    grow_instability!(simulation, energy)

Grow an instability by running `simulation`.

Estimates the growth rate ``σ`` of the instability
using the fractional change in volume-mean kinetic energy,
over the course of the `simulation`

``
energy(t₀ + Δτ) / energy(t₀) ≈ exp(2 σ Δτ)
``

where ``t₀`` is the starting time of the simulation and ``t₀ + Δτ``
the ending time of the simulation. We thus find that the growth rate
is measured by

``
σ = log(energy(t₀ + Δτ) / energy(t₀)) / (2 Δτ) .
``
"""
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
    rescale_factor = √(target_kinetic_energy / energy[1, 1, 1])

    for f in merge(model.velocities, model.tracers)
        f .*= rescale_factor
    end

    return nothing
end

using Printf
"""
    convergence(σ)

Check if the growth rate has converged. If the array `σ` has at least 2 elements then returns the
relative difference between ``σ[end]`` and ``σ[end-1]``.
"""
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

        @info @sprintf("About to start power method iteration %d; kinetic energy: %.2e", length(σ)+1, energy[1, 1, 1])
        push!(σ, grow_instability!(simulation, energy))
        compute!(energy)

        @info @sprintf("Power method iteration %d, kinetic energy: %.2e, σⁿ: %.2e, relative Δσ: %.2e",
                       length(σ), energy[1, 1, 1], σ[end], convergence(σ))

        compute!(ω)
        rescale!(simulation.model, energy)
        push!(power_method_data, (ω=collect(interior(ω, :, 1, :)), b=collect(interior(b, :, 1, :)), σ=deepcopy(σ)))
    end

    return σ, power_method_data
end
u, v, w = model.velocities
b = model.tracers.b

perturbation_vorticity = Field(∂z(u) - ∂x(w))

xω, yω, zω = nodes(perturbation_vorticity)
xb, yb, zb = nodes(b)
using Random, Statistics

mean_perturbation_kinetic_energy = Field(Average(1/2 * (u^2 + w^2)))
noise(x, z) = randn()
set!(model, u=noise, w=noise, b=noise)
rescale!(simulation.model, mean_perturbation_kinetic_energy, target_kinetic_energy=1e-6)
growth_rates, power_method_data = estimate_growth_rate(simulation, mean_perturbation_kinetic_energy, perturbation_vorticity, b)

@info "Power iterations converged! Estimated growth rate: $(growth_rates[end])"

# Reset the clock
model.clock.iteration = 0
model.clock.time = 0

estimated_growth_rate = growth_rates[end]

simulation.stop_time = 5 / estimated_growth_rate
simulation.stop_iteration = 9.1e18 # pretty big (not Inf tho)

# Rescale the eigenmode
initial_eigenmode_energy = 5e-5
rescale!(simulation.model, mean_perturbation_kinetic_energy, target_kinetic_energy=initial_eigenmode_energy)

total_vorticity = Field(∂z(u) + ∂z(model.background_fields.velocities.u) - ∂x(w))

total_b = Field(b + model.background_fields.tracers.b)

simulation.output_writers[:vorticity] =
    JLD2OutputWriter(model, (ω=perturbation_vorticity, Ω=total_vorticity, b=b, B=total_b, KE=mean_perturbation_kinetic_energy),
                     schedule = TimeInterval(0.10 / estimated_growth_rate),
                     filename = "kelvin_helmholtz_instability.jld2",
                     overwrite_existing = true)

@info "*** Running a simulation of Kelvin-Helmholtz instability..."
run!(simulation)