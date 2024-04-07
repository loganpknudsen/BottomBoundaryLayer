using Oceananigans
using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.Grids: Center, Face
using Oceananigans.Units
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

Lx = 2000meters
Lz = 200meters
Nx = 500
Nz = 100

grid = RectilinearGrid(arch; topology = (Periodic, Flat, Bounded),
                       size = (Nx, Nz),
                       x = (0, Lx),
                    #    z = z_faces
                       z = (0,Lz))

# grid = RectilinearGrid(arch; size=(1024, 200), y=(0,3000),z=(-200,0), topology=(Flat, Periodic, Bounded))

# tilted domain parameters
θ = 10^(-2) # degrees 
# ĝ = [θ, 0, 1] # gravity vector small angle
ĝ = [sind(θ), 0, cosd(θ)] # gravity vector

# realustic mid latitude for now
buoyancy = Buoyancy(model = BuoyancyTracer(), gravity_unit_vector = -ĝ)
coriolis = ConstantCartesianCoriolis(f = 1e-4, rotation_axis = ĝ)

# parameters
V∞ = 0.1 # m s⁻¹
N² = 1e-5 # interior stratification
f = coriolis.fz
ϕ = 0
S∞ = (N²*θ^2)/(f^2)
γ = (1+S∞)^(-1)
hu = (f*V∞)/(γ*N²*θ) # set to negative
uₒ = 0 #γ*(N²*θ)/(f)*cos(ϕ)
vₒ = γ*(θ*N²)/(f) #γ*(N²*θ)/(f)*sin(ϕ)
bₒ = 0 # initial stratification
fˢ=(f^2+θ^2*N²)^(0.5)

p =(;N²,θ,f,V∞,hu,γ,uₒ,vₒ,bₒ,fˢ)

# background flow with geostrophic and ageostrophic shear 

@inline interval(x,a,b) = ifelse(a<=x<=b, one(x), zero(x))

sn_fn(x,z,t,p) = sin(p.fˢ*t)
cs_fn(x,z,t,p) = cos(p.fˢ*t)

u_pert(x,z,t,p) = p.uₒ*cs_fn(x,z,t,p) -(p.f*p.vₒ+p.bₒ*p.θ)/(p.fˢ)*sn_fn(x,z,t,p)
v_pert(x,z,t,p) = (p.f^2*p.vₒ+p.f*p.bₒ*p.θ)/(p.fˢ)^2*cs_fn(x,z,t,p) -(p.f*p.uₒ)/(p.fˢ)*sn_fn(x,z,t,p)+((p.fˢ^2- p.f^2)*p.vₒ-p.f*p.bₒ*p.θ)/(p.fˢ)^2
b_pert(x,z,t,p) = p.N²*p.θ*(p.f*p.vₒ+p.bₒ*p.θ)/(p.fˢ)^2*cs_fn(x,z,t,p) -(p.N²*p.θ*p.uₒ)/(p.fˢ)*sn_fn(x,z,t,p)+((1- p.N²*p.θ)*p.bₒ-p.f*p.vₒ*p.θ*p.N²)/(p.fˢ)^2

u_adjustment(x, z, t, p) =  u_pert(x,z,t,p)*(p.hu-z)*interval(z,0,p.hu)
v_adjustment(x, z, t, p) = -p.γ*(p.θ * p.N²)/(p.f)*(p.hu-z)*interval(z,0,p.hu)+p.V∞a + v_pert(x,z,t,p)*(p.hu-z)*interval(z,0,p.hu)
constant_stratification(x, z, t, p) = p.N²*x*p.θ + p.N²*z - p.N²*p.γ*(p.hu-z)*interval(z,0,p.hu)+ b_pert(x,z,t,p)*(p.hu-z)*interval(z,0,p.hu)

U_field = BackgroundField(u_adjustment, parameters=p)
V_field = BackgroundField(v_adjustment, parameters=p)
B_field = BackgroundField(constant_stratification, parameters=p)

# Boundary condition set up

z₀ = 0.1 # m (roughness length)
κ = 0.4 # von Karman constant
z₁ = znodes(grid, Center())[1] # Closest grid center to the bottom
cᴰ = (κ / log(z₁ / z₀))^2 # Drag coefficient

@inline drag_u(x, t, u, v, p) = - p.cᴰ * √(u^2 + (v -p.γ*(p.θ * p.N²)/(p.f)*(p.hu)+ p.V∞)^2) * u
@inline drag_v(x, t, u, v, p) = - p.cᴰ * √(u^2 + (v -p.γ*(p.θ * p.N²)/(p.f)*(p.hu)+ p.V∞)^2) * (v-p.γ*(p.θ * p.N²)/(p.f)*(p.hu) + p.V∞)

ps = (;cᴰ, V∞, hu, θ, f, N², γ)

drag_bc_u = FluxBoundaryCondition(drag_u, field_dependencies=(:u, :v), parameters=ps)
drag_bc_v = FluxBoundaryCondition(drag_v, field_dependencies=(:u, :v), parameters=ps)

u_bcs = FieldBoundaryConditions(bottom = drag_bc_u)
v_bcs = FieldBoundaryConditions(bottom = drag_bc_v)

# boundary_conditions=(;b=buoyancy_grad),
closure = ScalarDiffusivity(ν=1e-4, κ=1e-4)

start_time = time_ns()

model = NonhydrostaticModel(; grid, buoyancy, coriolis, closure,
                            timestepper = :RungeKutta3,
                            advection = CenteredFourthOrder(),
                            tracers = :b,
                            boundary_conditions = (; u=u_bcs, v=v_bcs),
                            background_fields = (; u=U_field, v=V_field, b=B_field))

ns = 10^(-4) # standard deviation for noise

u₀(x, z) = ns*Random.randn()
v₀(x, z) = ns*Random.randn()
w₀(x, z) = ns*Random.randn()
# bₒ(x,y,z) = 0.005*Random.randn()

set!(model, u=u₀, v=v₀, w=w₀)

simulation = Simulation(model, Δt = 1, stop_time = 20*(2*pi)/f)


wizard = TimeStepWizard(cfl=0.7, max_change=1.1, max_Δt=10.0, min_Δt=0.001) 
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5)) 

progress_message(sim) =
        @printf("i: %04d, t: %s, Δt: %s, wall time: %s\n",
        sim.model.clock.iteration, prettytime(sim.model.clock.time),
        prettytime(sim.Δt), prettytime((time_ns() - start_time) * 1e-9))

simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(1*(2*pi)/f))

# and add an output writer that saves the vertical velocity field every two iterations:

u, v, w = model.velocities
ub = model.background_fields.velocities.u
vb = model.background_fields.velocities.v
b = model.tracers.b
B∞ = model.background_fields.tracers.b

U = u + ub
V = v + vb #+ V∞
B = b + B∞
# dBdz = Field(@at (Center, Center, Center) ∂z(b+B∞))


output = (; u, U, v, V, w, b, B)

# u,v,w = model.velocities

# output = (;u,v,w,model.tracers.b,U=(model.background_fields.velocities.u+0*u),V=(model.background_fields.velocities.v+0*v),B=(model.background_fields.tracers.b+0*model.tracers.b))

# ε = Field(KineticEnergyDissipationRate(model))
# dBdz = Field(@at (Center, Center, Center) ∂z(model.tracers.b+model.background_fields.tracers.b))
# u_m_flux = u*w
# v_m_flux = v*w

# output = merge(output, (; E=ε, N2=dBdz, UM=u_m_flux, VM=v_m_flux,))

simulation.output_writers[:fields] = NetCDFOutputWriter(model, output;
                                                          schedule = TimeInterval(0.1*(2*pi)/f),
                                                          filename = path_name*"BLL_w_O_test.nc",
                                                          overwrite_existing = true)

# With initial conditions set and an output writer at the ready, we run the simulation
 
run!(simulation)
