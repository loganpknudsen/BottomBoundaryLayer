# # Internal tide by a seamount
#
# In this example, we show how internal tide is generated from a barotropic tidal flow
# sloshing back and forth over a sea mount.
#
# ## Install dependencies
#
# First let's make sure we have all required packages installed.

# ```julia
# using Pkg
# pkg"add Oceananigans, CairoMakie"
# ```

using Oceananigans
using Oceananigans.Units

# ## Grid

# We create an `ImmersedBoundaryGrid` wrapped around an underlying two-dimensional `RectilinearGrid`
# that is periodic in ``x`` and bounded in ``z``.

Nx, Nz = 250, 125

H = 2kilometers

underlying_grid = RectilinearGrid(size = (Nx, Nz),
                                  x = (-1000kilometers, 1000kilometers),
                                  z = (-H, 0),
                                  halo = (4, 4),
                                  topology = (Periodic, Flat, Bounded))

# Now we can create the non-trivial bathymetry. We use `GridFittedBottom` that gets as input either
# *(i)* a two-dimensional function whose arguments are the grid's native horizontal coordinates and
# it returns the ``z`` of the bottom, or *(ii)* a two-dimensional array with the values of ``z`` at
# the bottom cell centers.
#
# In this example we'd like to have a Gaussian hill at the center of the domain.
#
# ```math
# h(x) = -H + h_0 \exp(-x^2 / 2σ^2)
# ```

h₀ = 250meters
width = 20kilometers
hill(x) = h₀ * exp(-x^2 / 2width^2)
bottom(x) = - H + hill(x)

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom))

# Let's see how the domain with the bathymetry is.

x = xnodes(grid, Center())
bottom_boundary = interior(grid.immersed_boundary.bottom_height, :, 1, 1)
top_boundary = 0*x

using CairoMakie

fig = Figure(size = (700, 200))
ax = Axis(fig[1, 1],
          xlabel="x [km]",
          ylabel="z [m]",
          limits=((-grid.Lx/2e3, grid.Lx/2e3), (-grid.Lz, 0)))

band!(ax, x/1e3, bottom_boundary, top_boundary, color = :mediumblue)

fig

# Now we want to add a barotropic tide forcing. For example, to add the lunar semi-diurnal ``M_2`` tide 
# we need to add forcing in the ``u``-momentum equation of the form:
# ```math
# F_0 \sin(\omega_2 t)
# ```
# where ``\omega_2 = 2π / T_2``, with ``T_2 = 12.421 \,\mathrm{hours}`` the period of the ``M_2`` tide.

# The excursion parameter is a nondimensional number that expresses the ratio of the flow movement
# due to the tide compared to the size of the width of the hill.
#
# ```math
# \epsilon = \frac{U_{\mathrm{tidal}} / \omega_2}{\sigma}
# ```
# 
# We prescribe the excursion parameter which, in turn, implies a tidal velocity ``U_{\mathrm{tidal}}``
# which then allows us to determing the tidal forcing amplitude ``F_0``. For the last step, we
# use Fourier decomposition on the inviscid, linearized momentum equations to determine the
# flow response for a given tidal forcing. Doing so we get that for the sinusoidal forcing above,
# the tidal velocity and tidal forcing amplitudes are related via:
#
# ```math
# U_{\mathrm{tidal}} = \frac{\omega_2}{\omega_2^2 - f^2} F_0
# ```
#
# Now we have the way to find the value of the tidal forcing amplitude that would correspond to a
# given excursion parameter. The Coriolis frequency is needed, so we start by constructing a Coriolis on an ``f``-plane at the
# mid-latitudes.

coriolis = FPlane(latitude = -45)

# Now we have everything we require to construct the tidal forcing given a value of the
# excursion parameter.

T₂ = 12.421hours
ω₂ = 2π / T₂ # radians/sec

ϵ = 0.1 # excursion parameter

U_tidal = ϵ * ω₂ * width

tidal_forcing_amplitude = U_tidal * (ω₂^2 - coriolis.f^2) / ω₂

@inline tidal_forcing(x, z, t, p) = p.tidal_forcing_amplitude * sin(p.ω₂ * t)

u_forcing = Forcing(tidal_forcing, parameters=(; tidal_forcing_amplitude, ω₂))

# ## Model

# We built a `HydrostaticFreeSurfaceModel`:

model = HydrostaticFreeSurfaceModel(; grid, coriolis,
                                      buoyancy = BuoyancyTracer(),
                                      tracers = :b,
                                      momentum_advection = WENO(),
                                      tracer_advection = WENO(),
                                      forcing = (; u = u_forcing))

# We initialize the model with the tidal flow and a linear stratification.

uᵢ(x, z) = U_tidal

Nᵢ² = 1e-4  # [s⁻²] initial buoyancy frequency / stratification
bᵢ(x, z) = Nᵢ² * z

set!(model, u=uᵢ, b=bᵢ)

# Now let's built a `Simulation`.

Δt = 5minutes
stop_time = 4days

simulation = Simulation(model; Δt, stop_time)

# We add a callback to print a message about how the simulation is going,

using Printf

wall_clock = Ref(time_ns())

function progress(sim)
    elapsed = 1e-9 * (time_ns() - wall_clock[])

    msg = @sprintf("iteration: %d, time: %s, wall time: %s, max|w|: %6.3e, m s⁻¹\n",
                   iteration(sim), prettytime(sim), prettytime(elapsed),
                   maximum(abs, sim.model.velocities.w))

    wall_clock[] = time_ns()

    @info msg

    return nothing
end

add_callback!(simulation, progress, name=:progress, IterationInterval(200))
nothing #hide

# ## Diagnostics/Output

# Add some diagnostics. Instead of ``u`` we save the deviation of ``u`` from its instantaneous
# domain average, ``u' = u - (L_x H)^{-1} \int u \, \mathrm{d}x \mathrm{d}z``. We also save
# the stratification ``N^2 = \partial_z b``.

b = model.tracers.b
u, v, w = model.velocities

U = Field(Average(u))

u′ = u - U

N² = ∂z(b)

filename = "internal_tide"
save_fields_interval = 30minutes

simulation.output_writers[:fields] = JLD2OutputWriter(model, (; u, u′, w, b, N²);
                                                      filename,
                                                      schedule = TimeInterval(save_fields_interval),
                                                      overwrite_existing = true)

# We are ready -- let's run!

run!(simulation)