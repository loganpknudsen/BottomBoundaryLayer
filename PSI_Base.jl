## Note for self: I will be building off of this program to recreate the simulations in Thomas and Taylor 2014
#
# # Internal wave example
#
# ## The physical domain
#

using Oceananigans
using Random
using LaTeXStrings

# made grid correct shape, need to modify z boundaries to make sure they are no slip
grid = RectilinearGrid(size=(1, 1024, 200), x=(0,1),y=(0,3000),z=(-200,0), topology=(Periodic, Periodic, Bounded))

# ## Internal wave parameters
#
# Inertia-gravity waves propagate in fluids that are both _(i)_ rotating, and
# _(ii)_ density-stratified. We use Oceananigans' Coriolis abstraction
# to implement a background rotation rate:

# realustuc mid latitude for now
coriolis = FPlane(rotation_rate=7.292115e-5, latitude=45)

## initial buoyancy frequency, horizontal buoyancy, scaling factor, inital phase of inertial oscillation, corilois
ps = (Nₒ = 81.7*coriolis.f, S = 7.7*coriolis.f, γ = 0.6, ϕ = 0, f = coriolis.f)

# background flow with geostrophic and ageostrophic shear 
U_func(x, y, z, t, ps) = ((ps.S^2*z)/ps.f)*(1+ps.γ*cos(ps.f*t-ps.ϕ))
V_func(x, y, z, t, ps) = -1*((ps.S^2*z*ps.γ)/ps.f)*(sin(ps.f*t-ps.ϕ))
W_func(x,y,z,t) = 0
B_func(x, y, z, t, ps) = (ps.Nₒ^2-ps.γ*(ps.S^4/ps.f^2)*(cos(ps.ϕ)-cos(ps.f*t-ps.ϕ)))*z - ps.S^2*y #multiply by z since we integrate N^2 w.r.t z
U = BackgroundField(U_func, parameters=ps)
V = BackgroundField(V_func, parameters=ps)
W = BackgroundField(W_func)
B = BackgroundField(B_func, parameters=ps)

# Boundary condition set up
no_slip_field_bcs = FieldBoundaryConditions(FluxBoundaryCondition(0.0));
buoyancy_grad = FieldBoundaryConditions(top=GradientBoundaryCondition(ps.Nₒ),bottom=GradientBoundaryCondition(ps.Nₒ))

model = NonhydrostaticModel(; grid,
                            boundary_conditions=(u=no_slip_field_bcs, v=no_slip_field_bcs, w=no_slip_field_bcs, b=buoyancy_grad),
                            coriolis,
                            advection = CenteredFourthOrder(),
                            timestepper = :RungeKutta3,
                            closure = ScalarDiffusivity(ν=1e-3,κ=1e-3), # removed molecular diffusiviy 
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            background_fields = (; u=U, v=V, w=W, b=B)) # `background_fields` is a `NamedTuple`

u₀(x, y, z) = 0.25*Random.randn()
v₀(x, y, z) = 0.25*Random.randn()
w₀(x, y, z) = 0.25*Random.randn()
b₀(x, y, z) = 0

set!(model, u=u₀, v=v₀, w=w₀, b=b₀)

simulation = Simulation(model, Δt = 1, stop_time = 10)

wizard = TimeStepWizard(cfl=1, max_change=1.1, max_Δt=10.0, min_Δt=0.0001) # dec cfl 0.9 -> .5
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5)) # dec. int 500 -> 5 (max)

# and add an output writer that saves the vertical velocity field every two iterations:

filename = "/Users/loganknudsen/Documents/UMD_Research/BottomBoundaryLayer/internal_wave.jld2"
simulation.output_writers[:velocities] = JLD2OutputWriter(model, model.velocities; filename,
                                                          schedule = TimeInterval(1),
                                                          overwrite_existing = true)

# With initial conditions set and an output writer at the ready, we run the simulation

run!(simulation)

# ## Animating a propagating packet
#
# To animate a the propagating wavepacket we just simulated, we load CairoMakie
# and make a Figure and an Axis for the animation,
