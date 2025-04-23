module TKEBudgetTerms
using DocStringExtensions

export TurbulentKineticEnergy, KineticEnergy
export KineticEnergyTendency
export AdvectionTerm
export KineticEnergyStressTerm
export KineticEnergyForcingTerm
export IsotropicKineticEnergyDissipationRate, KineticEnergyDissipationRate
export PressureRedistributionTerm
export BuoyancyProductionTerm
export XShearProductionRate, YShearProductionRate, ZShearProductionRate

using Oceananigans: NonhydrostaticModel, HydrostaticFreeSurfaceModel, fields
using Oceananigans.Operators
using Oceananigans.AbstractOperations
using Oceananigans.AbstractOperations: KernelFunctionOperation
using Oceananigans.Grids: Center, Face
using Oceananigans.Fields: ZeroField
using Oceananigans.Models.NonhydrostaticModels: u_velocity_tendency, v_velocity_tendency, w_velocity_tendency
using Oceananigans.Advection: div_𝐯u, div_𝐯v, div_𝐯w
using Oceananigans.TurbulenceClosures: viscous_flux_ux, viscous_flux_uy, viscous_flux_uz,
                                       viscous_flux_vx, viscous_flux_vy, viscous_flux_vz,
                                       viscous_flux_wx, viscous_flux_wy, viscous_flux_wz,
                                       ∂ⱼ_τ₁ⱼ, ∂ⱼ_τ₂ⱼ, ∂ⱼ_τ₃ⱼ
using Oceananigans.TurbulenceClosures: immersed_∂ⱼ_τ₁ⱼ, immersed_∂ⱼ_τ₂ⱼ, immersed_∂ⱼ_τ₃ⱼ
using Oceananigans.BuoyancyFormulations: x_dot_g_bᶠᶜᶜ, y_dot_g_bᶜᶠᶜ, z_dot_g_bᶜᶜᶠ

using Oceanostics: _νᶜᶜᶜ
using Oceanostics: validate_location, validate_dissipative_closure, perturbation_fields


# Some useful operators
@inline ψ²(i, j, k, grid, ψ) = @inbounds ψ[i, j, k]^2

@inline ψ′²(i, j, k, grid, ψ, ψ̄) = @inbounds (ψ[i, j, k] - ψ̄[i, j, k])^2
@inline ψ′²(i, j, k, grid, ψ, ψ̄::Number) = @inbounds (ψ[i, j, k] - ψ̄)^2

@inline fψ²(i, j, k, grid, f, ψ) = f(i, j, k, grid, ψ)^2

@inline fψ_plus_gφ²(i, j, k, grid, f, ψ, g, φ) = (f(i, j, k, grid, ψ) + g(i, j, k, grid, φ))^2

@inline ψf(i, j, k, grid, ψ, f, args...) = @inbounds ψ[i, j, k] * f(i, j, k, grid, args...)

@inline function uᵢ∂ⱼ_τᵢⱼᶜᶜᶜ(i, j, k, grid, closure,
                                            diffusivity_fields,
                                            clock,
                                            model_fields,
                                            buoyancy,
                                            velocities)

    u∂ⱼ_τ₁ⱼ = ℑxᶜᵃᵃ(i, j, k, grid, ψf, velocities.u, ∂ⱼ_τ₁ⱼ, closure, diffusivity_fields, clock, model_fields, buoyancy)
    v∂ⱼ_τ₂ⱼ = ℑyᵃᶜᵃ(i, j, k, grid, ψf, velocities.v, ∂ⱼ_τ₂ⱼ, closure, diffusivity_fields, clock, model_fields, buoyancy)
    w∂ⱼ_τ₃ⱼ = ℑzᵃᵃᶜ(i, j, k, grid, ψf, velocities.w, ∂ⱼ_τ₃ⱼ, closure, diffusivity_fields, clock, model_fields, buoyancy)

    return u∂ⱼ_τ₁ⱼ+ v∂ⱼ_τ₂ⱼ + w∂ⱼ_τ₃ⱼ
end

"""
    $(SIGNATURES)

Return a `KernelFunctionOperation` that computes the diffusive term of the KE prognostic equation:

```
DIFF = uᵢ∂ⱼτᵢⱼ
```

where `uᵢ` are the velocity components and `τᵢⱼ` is the diffusive flux of `i` momentum in the
`j`-th direction.

"""
function KineticEnergyStress(model; velocities=model.velocities ,location = (Center, Center, Center))
    validate_location(location, "KineticEnergyStressTerm")
    model_fields = fields(model)

    if model isa HydrostaticFreeSurfaceModel
        model_fields = (; model_fields..., w=ZeroField())
    end
    dependencies = (model.closure,
                    model.diffusivity_fields,
                    model.clock,
                    fields(model),
                    model.buoyancy,
                    velocities)
    return KernelFunctionOperation{Center, Center, Center}(uᵢ∂ⱼ_τᵢⱼᶜᶜᶜ, model.grid, dependencies...)
end