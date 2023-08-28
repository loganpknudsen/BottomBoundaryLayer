using Oceananigans
using Random
using Printf
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "path"
        help = "pathname to save data under"
        default = ""
    end
    return parse_args(s)
end

args=parse_commandline()

@info("command line args:")
for (arg,val) in args
  @info("   $arg => $val")
end

path_name = args["path"]

# made grid correct shape, need to modify z boundaries to make sure they are no slip
grid = RectilinearGrid(size=(1, 1024, 200), x=(0,1),y=(0,3000),z=(-200,0), topology=(Periodic, Periodic, Bounded))

# realustuc mid latitude for now
coriolis = FPlane(rotation_rate=7.292115e-5, latitude=45)

## initial buoyancy frequency, horizontal buoyancy, scaling factor, inital phase of inertial oscillation, corilois
ps = (Nₒ = 81.7*coriolis.f, S = 7.7*coriolis.f, γ = 0.6, ϕ = 0, f = coriolis.f)

# background flow with geostrophic and ageostrophic shear 
U_func(x, y, z, t, ps) = ((ps.S^2*z)/ps.f)*(1+ps.γ*cos(ps.f*t-ps.ϕ))
V_func(x, y, z, t, ps) = -1*((ps.S^2*z*ps.γ)/ps.f)*(sin(ps.f*t-ps.ϕ))
B_func(x, y, z, t, ps) = (ps.Nₒ^2-ps.γ*(ps.S^4/ps.f^2)*(cos(ps.ϕ)-cos(ps.f*t-ps.ϕ)))*z - ps.S^2*y #multiply by z since we integrate N^2 w.r.t z
U = BackgroundField(U_func, parameters=ps)
V = BackgroundField(V_func, parameters=ps)
B = BackgroundField(B_func, parameters=ps)

# Boundary condition set up
no_slip_field_bcs = FieldBoundaryConditions(FluxBoundaryCondition(0.0));
buoyancy_grad = FieldBoundaryConditions(top=GradientBoundaryCondition(ps.Nₒ^2),bottom=GradientBoundaryCondition(ps.Nₒ^2))

start_time = time_ns()

model = NonhydrostaticModel(; grid,
                            boundary_conditions=(u=no_slip_field_bcs, v=no_slip_field_bcs, w=no_slip_field_bcs, b=buoyancy_grad),
                            coriolis,
                            advection = CenteredFourthOrder(),
                            timestepper = :RungeKutta3,
                            closure = ScalarDiffusivity(ν=1e-3,κ=1e-3), # removed molecular diffusiviy 
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            background_fields = ( u=U, v=V, b=B)) # `background_fields` is a `NamedTuple`

u₀(x, y, z) = 0.01*Random.randn()
v₀(x, y, z) = 0.01*Random.randn()
w₀(x, y, z) = 0.01*Random.randn()
bₒ(x,y,z) = 0.005*Random.randn()

set!(model, u=u₀, v=v₀, w=w₀, b=bₒ)

simulation = Simulation(model, Δt = 1, stop_time = 100000)


wizard = TimeStepWizard(cfl=0.5, max_change=1.1, max_Δt=10.0, min_Δt=0.0001) 
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5)) 

progress_message(sim) =
        @printf("i: %04d, t: %s, Δt: %s, wall time: %s\n",
        sim.model.clock.iteration, prettytime(sim.model.clock.time),
        prettytime(sim.Δt), prettytime((time_ns() - start_time) * 1e-9))

simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(1000.0))

# and add an output writer that saves the vertical velocity field every two iterations:

u,v,w = model.velocities

output = (;u,v,w,model.tracers.b,U=(model.background_fields.velocities.u+0*u),V=(model.background_fields.velocities.v+0*v),B=(model.background_fields.tracers.b+0*model.tracers.b))

simulation.output_writers[:fields] = NetCDFOutputWriter(model, output;
                                                          schedule = TimeInterval(500.0),
                                                          filename = path_name*"PSI_100k.nc",
                                                          overwrite_existing = true)

# With initial conditions set and an output writer at the ready, we run the simulation
 
run!(simulation)
