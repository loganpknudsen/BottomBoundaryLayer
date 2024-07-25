# # Tilted bottom boundary layer example
#
# This example simulates a two-dimensional oceanic bottom boundary layer
# in a domain that's tilted with respect to gravity. We simulate the perturbation
# away from a constant along-slope (y-direction) velocity constant density stratification.
# This perturbation develops into a turbulent bottom boundary layer due to momentum
# loss at the bottom boundary modeled with a quadratic drag law.
# 
# This example illustrates
#
#   * changing the direction of gravitational acceleration in the buoyancy model;
#   * changing the axis of rotation for Coriolis forces.
#
# ## Install dependencies
#
# First let's make sure we have all required packages installed.

# ```julia
# using Pkg
# pkg"add Oceananigans, NCDatasets, CairoMakie"
# ```
#
# ## The domain
#
# We create a grid with finer resolution near the bottom,

using Oceananigans
using Oceananigans.Units

Lx = 200meters
Lz = 100meters
Nx = 64
Nz = 64

## Creates a grid with near-constant spacing `refinement * Lz / Nz`
## near the bottom:
refinement = 1.8 # controls spacing near surface (higher means finer spaced)
stretching = 10  # controls rate of stretching at bottom 

## "Warped" height coordinate
h(k) = (Nz + 1 - k) / Nz

## Linear near-surface generator
ζ(k) = 1 + (h(k) - 1) / refinement

## Bottom-intensified stretching function 
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))

## Generating function
z_faces(k) = - Lz * (ζ(k) * Σ(k) - 1)

grid = RectilinearGrid(topology = (Periodic, Flat, Bounded),
                       size = (Nx, Nz),
                       x = (0, Lx),
                       z = z_faces)

# Let's make sure the grid spacing is both finer and near-uniform at the bottom,

# ## Tilting the domain
#
# We use a domain that's tilted with respect to gravity by

θ = 3 # degrees

# so that ``x`` is the along-slope direction, ``z`` is the across-slope direction that
# is perpendicular to the bottom, and the unit vector anti-aligned with gravity is

ĝ = [sind(θ), 0, cosd(θ)]

# Changing the vertical direction impacts both the `gravity_unit_vector`
# for `Buoyancy` as well as the `rotation_axis` for Coriolis forces,

buoyancy = Buoyancy(model = BuoyancyTracer(), gravity_unit_vector = -ĝ)
coriolis = ConstantCartesianCoriolis(f = 1e-4, rotation_axis = ĝ)

# where above we used a constant Coriolis parameter ``f = 10^{-4} \, \rm{s}^{-1}``.
# The tilting also affects the kind of density stratified flows we can model.
# In particular, a constant density stratification in the tilted
# coordinate system

@inline constant_stratification(x, z, t, p) = p.N² * (x * p.ĝ[1] + z * p.ĝ[3])

# is _not_ periodic in ``x``. Thus we cannot explicitly model a constant stratification
# on an ``x``-periodic grid such as the one used here. Instead, we simulate periodic
# _perturbations_ away from the constant density stratification by imposing
# a constant stratification as a `BackgroundField`,

N² = 1e-5 # s⁻² # background vertical buoyancy gradient
B∞_field = BackgroundField(constant_stratification, parameters=(; ĝ, N² = N²))

# We choose to impose a bottom boundary condition of zero *total* diffusive buoyancy
# flux across the seafloor,
# ```math
# ∂_z B = ∂_z b + N^{2} \cos{\theta} = 0.
# ```
# This shows that to impose a no-flux boundary condition on the total buoyancy field ``B``, we must apply a boundary condition to the perturbation buoyancy ``b``,
# ```math
# ∂_z b = - N^{2} \cos{\theta}.
#```

∂z_b_bottom = - N² * cosd(θ)
negative_background_diffusive_flux = GradientBoundaryCondition(∂z_b_bottom)
b_bcs = FieldBoundaryConditions(bottom = negative_background_diffusive_flux)

# ## Bottom drag
#
# We impose bottom drag that follows Monin--Obukhov theory.
# We include the background flow in the drag calculation,
# which is the only effect the background flow enters the problem,

V∞ = 0.1 # m s⁻¹
z₀ = 0.1 # m (roughness length)
κ = 0.4  # von Karman constant

z₁ = first(znodes(grid, Center())) # Closest grid center to the bottom
cᴰ = (κ / log(z₁ / z₀))^2 # Drag coefficient

@inline drag_u(x, t, u, v, p) = - p.cᴰ * √(u^2 + (v + p.V∞)^2) * u
@inline drag_v(x, t, u, v, p) = - p.cᴰ * √(u^2 + (v + p.V∞)^2) * (v + p.V∞)

drag_bc_u = FluxBoundaryCondition(drag_u, field_dependencies=(:u, :v), parameters=(; cᴰ, V∞))
drag_bc_v = FluxBoundaryCondition(drag_v, field_dependencies=(:u, :v), parameters=(; cᴰ, V∞))

u_bcs = FieldBoundaryConditions(bottom = drag_bc_u)
v_bcs = FieldBoundaryConditions(bottom = drag_bc_v)

# ## Create the `NonhydrostaticModel`
#
# We are now ready to create the model. We create a `NonhydrostaticModel` with an
# `UpwindBiasedFifthOrder` advection scheme, a `RungeKutta3` timestepper,
# and a constant viscosity and diffusivity. Here we use a smallish value
# of ``10^{-4} \, \rm{m}^2\, \rm{s}^{-1}``.

ν = 1e-4
κ = 1e-4
closure = ScalarDiffusivity(ν=ν, κ=κ)

model = NonhydrostaticModel(; grid, buoyancy, coriolis, closure,
                            timestepper = :RungeKutta3,
                            advection = UpwindBiasedFifthOrder(),
                            tracers = :b,
                            boundary_conditions = (u=u_bcs, v=v_bcs, b=b_bcs),
                            background_fields = (; b=B∞_field))

# Let's introduce a bit of random noise at the bottom of the domain to speed up the onset of
# turbulence:

noise(x, z) = 1e-3 * randn() * exp(-(10z)^2 / grid.Lz^2)
set!(model, u=noise, w=noise)

# ## Create and run a simulation
#
# We are now ready to create the simulation. We begin by setting the initial time step
# conservatively, based on the smallest grid size of our domain and either an advective
# or diffusive time scaling, depending on which is shorter.

Δt₀ = 0.5 * minimum([minimum_zspacing(grid) / V∞, minimum_zspacing(grid)^2/κ])
simulation = Simulation(model, Δt = Δt₀, stop_time = 0.5day)

# We use a `TimeStepWizard` to adapt our time-step,

wizard = TimeStepWizard(max_change=1.1, cfl=0.7)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(4))

# and also we add another callback to print a progress message,

using Printf

start_time = time_ns() # so we can print the total elapsed wall time

progress_message(sim) =
    @printf("Iteration: %04d, time: %s, Δt: %s, max|w|: %.1e m s⁻¹, wall time: %s\n",
            iteration(sim), prettytime(time(sim)),
            prettytime(sim.Δt), maximum(abs, sim.model.velocities.w),
            prettytime((time_ns() - start_time) * 1e-9))

simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(200))

# ## Add outputs to the simulation
#
# We add outputs to our model using the `NetCDFOutputWriter`,

u, v, w = model.velocities
b = model.tracers.b
B∞ = model.background_fields.tracers.b

B = b + B∞
V = v + V∞
ωy = ∂z(u) - ∂x(w)

outputs = (; u, V, w, B, ωy)

simulation.output_writers[:fields] = NetCDFOutputWriter(model, outputs;
                                                        filename = joinpath(@__DIR__, "tilted_bottom_boundary_layer.nc"),
                                                        schedule = TimeInterval(20minutes),
                                                        overwrite_existing = true)

# Now we just run it!

run!(simulation)