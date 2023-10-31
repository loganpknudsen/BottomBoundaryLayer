using Oceananigans
using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.Grids: Center, Face
using Random
using Printf
using ArgParse
using CUDA: has_cuda_gpu
using Oceanostics

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

# grid specifications
arch = has_cuda_gpu() ? GPU() : CPU()

grid = RectilinearGrid(arch; size=(1024, 200), y=(0,3000),z=(-200,0), topology=(Flat, Periodic, Bounded))

# realustuc mid latitude for now
coriolis = FPlane(rotation_rate=7.292115e-5, latitude=45)

## initial buoyancy frequency, horizontal buoyancy, scaling factor, inital phase of inertial oscillation, corilois
ps = (Nₒ = 81.7*coriolis.f, S = 7.7*coriolis.f, γ =0.6, ϕ = 0, f = coriolis.f)

# background flow with geostrophic and ageostrophic shear 
@inline cs_func(t,ps) = cos(ps.f*t-ps.ϕ)
@inline sn_func(t,ps) = sin(ps.f*t-ps.ϕ)
@inline phs_dff(t,ps) = cos(ps.ϕ)-cs_func(t,ps)

U_func(x, y, z, t, ps) = (ps.S^2/ps.f)*(1+ps.γ*cs_func(t,ps))*z # current run is set on gamma=0.6
V_func(x, y, z, t, ps) = -1*((ps.S^2*ps.γ)/ps.f)*sn_func(t,ps)*z # change to 0.6 on next run if current on collapses
B_func(x, y, z, t, ps) = (ps.Nₒ^2-ps.γ*(ps.S^4/ps.f^2)*phs_dff(t,ps))*z - ps.S^2*y #multiply by z since we integrate N^2 w.r.t z
U = BackgroundField(U_func, parameters=ps)
V = BackgroundField(V_func, parameters=ps)
B = BackgroundField(B_func, parameters=ps)

# Boundary condition set up

b_bc = GradientBoundaryCondition(ps.Nₒ^2)
buoyancy_grad = FieldBoundaryConditions(top=b_bc,bottom=b_bc)
# boundary_conditions=(;b=buoyancy_grad),

U = (ps.S^2*ps.γ*200)/(coriolis.f)
eddy_visc = (U*200)/(1.3*10^6)
diffus = eddy_visc

start_time = time_ns()

model = NonhydrostaticModel(; grid,
                            boundary_conditions=(;b=buoyancy_grad),
                            coriolis,
                            advection = CenteredFourthOrder(),
                            timestepper = :RungeKutta3,
                            closure = ScalarDiffusivity(ν=eddy_visc,κ=diffus), # removed molecular diffusiviy 
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            background_fields = (; u=U, v=V, b=B)) # `background_fields` is a `NamedTuple`

ns = 10^(-4) # standard deviation for noise

u₀(x, y, z) = ns*Random.randn()
v₀(x, y, z) = ns*Random.randn()
w₀(x, y, z) = ns*Random.randn()
# bₒ(x,y,z) = 0.005*Random.randn()

set!(model, u=u₀, v=v₀, w=w₀)

simulation = Simulation(model, Δt = 1, stop_time = 20*(2*pi)/ps.f)


wizard = TimeStepWizard(cfl=0.5, max_change=1.1, max_Δt=10.0, min_Δt=0.001) 
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5)) 

progress_message(sim) =
        @printf("i: %04d, t: %s, Δt: %s, wall time: %s\n",
        sim.model.clock.iteration, prettytime(sim.model.clock.time),
        prettytime(sim.Δt), prettytime((time_ns() - start_time) * 1e-9))

simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(1*(2*pi)/ps.f))

# and add an output writer that saves the vertical velocity field every two iterations:

u,v,w = model.velocities

output = (;u,v,w,model.tracers.b,U=(model.background_fields.velocities.u+0*u),V=(model.background_fields.velocities.v+0*v),B=(model.background_fields.tracers.b+0*model.tracers.b))

ε = Field(KineticEnergyDissipationRate(model))
dBdz = Field(@at (Center, Center, Center) ∂z(model.tracers.b+model.background_fields.tracers.b))

output = merge(output, (; E=ε, N2=dBdz,))

simulation.output_writers[:fields] = NetCDFOutputWriter(model, output;
                                                          schedule = TimeInterval(0.05*(2*pi)/ps.f),
                                                          filename = path_name*"psi_base_test_ocng_w_b.nc",
                                                          overwrite_existing = true)

# With initial conditions set and an output writer at the ready, we run the simulation
 
run!(simulation)
