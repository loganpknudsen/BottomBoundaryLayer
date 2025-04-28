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
@inline ψ′(i, j, k, grid, ψ, ψ̄) = @inbounds (ψ[i, j, k] - ψ̄[i, j, k])
@inline ψ′²(i, j, k, grid, ψ, ψ̄::Number) = @inbounds (ψ[i, j, k] - ψ̄)^2

@inline fψ²(i, j, k, grid, f, ψ) = f(i, j, k, grid, ψ)^2

@inline fψ_plus_gφ²(i, j, k, grid, f, ψ, g, φ) = (f(i, j, k, grid, ψ) + g(i, j, k, grid, φ))^2

@inline ψf(i, j, k, grid, ψ, f) = @inbounds ψ[i, j, k] * f[i, j, k]

@inline function TRNS(i, j, k, grid, closure,
                                            diffusivity_fields,
                                            clock,
                                            model_fields,
                                            buoyancy,
                                            mean_velocities)
    k = (ℑxᶜᵃᵃ(i, j, k, grid, ψ′², model.velocities.u, mean_velocities.U) + ℑyᵃᶜᵃ(i, j, k, grid, ψ′², model.velocities.v, mean_velocities.V) + ℑzᵃᵃᶜ(i, j, k, grid, ψ′², model.velocities.w, mean_velocities.W)) / 2
    wpert = ℑzᵃᵃᶜ(i, j, k, grid, ψ′, model.velocities.w, mean_velocities.W)
    wk = ψf(i,j,k,grid, wpert,k)
    dwkdz = ∂zᶜᶜᶜ(i, j, k, grid, wk)
    return dwkdz
end


function KineticEnergyTransport(model; U=ZeroField(),V=ZeroField(),W=ZeroField(),location = (Center, Center, Center))
    validate_location(location, "KineticEnergyTransport")
    model_fields = fields(model)

    dependencies = (model.closure,
                    model.diffusivity_fields,
                    model.clock,
                    model_fields,
                    model.buoyancy,
                    mean_velocities=(U,V,W))
    return KernelFunctionOperation{Center, Center, Center}(TRNS, model.grid, dependencies...)
end