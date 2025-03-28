### Load in Packages
using Pkg
Pkg.activate(".")
using Oceananigans
using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.AbstractOperations: @at, Average
using Oceananigans.Grids: Center, Face
using Oceananigans.Fields: FunctionField
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

path_name = "/Users/loganknudsen/Documents/GitHub/BottomBoundaryLayer" #args["path"]


# grid specifications
arch = has_cuda_gpu() ? GPU() : CPU()
@info("Arch => $arch")

Lx = 200meters
Lz = 200meters
Nx = 56 # 512 originally
Nz = 56 # # 128 originally Note to self, maintain 2 to 1 resolution ration

grid = RectilinearGrid(arch; topology = (Periodic, Flat, Bounded),
                       size = (Nx, Nz),
                       x = (0, Lx),
                       z = (0,Lz))


# tilted domain parameters
const θ = 1e-1 # degrees 10^(-2) is previous value for 110 meter layer
ĝ = [sind(θ), 0, cosd(θ)] # gravity vector

# realistic mid latitude for now
buoyancy = Buoyancy(model = BuoyancyTracer(), gravity_unit_vector = -ĝ)
coriolis = ConstantCartesianCoriolis(f = 1e-4, rotation_axis = ĝ)

# parameters for simulation
const V∞ = 0.05 # m s⁻¹ interior velocity
const f = 1e-4 # coriolis parameter
const N² = 1e-6 # interior stratification
const S∞ = (N²*θ^2)/(f^2) # slope burger number
const fˢ = (f^2+θ^2*N²)^(0.5) # modified oscillation
const δ = 0.5
const γ = ((1+S∞*(1-δ))^(-1)+(3-S∞)*(3*(1+S∞)-S∞*(1-δ))^(-1))/2 # 0 PV parameter
const hu = (f*V∞)/(γ*N²*θ) # Height of Boundary Layer

# array of paramerers for background function
p =(; N², θ, f, V∞, hu, γ)

# heaviside function for boundary layer
heaviside(x,z) = 0.5*(1+tanh(10000*z))

u_adjustment(x, z, t, p) = 0
v_adjustment(x, z, t, p) = p.V∞ - p.γ*(p.θ * p.N²)/(p.f)*(p.hu-z)*heaviside(x,p.hu-z) 
constant_stratification(x, z, t, p) = p.N²*x*p.θ + p.N²*z + p.N²*p.γ*(p.hu-z)*heaviside(x,p.hu-z)

U_field = BackgroundField(u_adjustment, parameters=p)
V_field = BackgroundField(v_adjustment, parameters=p)
B_field = BackgroundField(constant_stratification, parameters=p)

b_bc_top= GradientBoundaryCondition(-1*N²)
# b_bc_bottom= GradientBoundaryCondition(N²)
buoyancy_grad = FieldBoundaryConditions(top=b_bc_top) # ,bottom=b_bc_bottom

# diffusitivity and viscosity values for closure
const ν1 = 5e-5
closure = ScalarDiffusivity(ν=ν1, κ=ν1)

start_time = time_ns()

# model set up 
model = NonhydrostaticModel(; grid, buoyancy, coriolis, closure,
                            timestepper = :RungeKutta3,
                            advection = WENO(),
                            tracers = :b,
                            boundary_conditions = (; b=buoyancy_grad),
                            background_fields = (; u=U_field, v=V_field, b=B_field))

ns = 10^(-4) # standard deviation for noise

# initial conditions to start instability
ui(x, z) = ns*Random.randn()
vi(x, z) = ns*Random.randn()
wi(x, z) = ns*Random.randn()
# bp(x,z) = ns*Random.randn()

# set simulation and decide run time
set!(model, u=ui, v=vi, w=wi)

simulation = Simulation(model, Δt = 1seconds, stop_time = 0.05*((2*pi)/fˢ)seconds) # stop_iteration=10

# time step wizard
wizard = TimeStepWizard(cfl=1, max_change=1.1seconds, max_Δt=100.0seconds, min_Δt=0.01seconds) 
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5)) 

# simulation.output_writers[:checkpointer] = Checkpointer(model; schedule=TimeInterval((5*(2*pi)/f)seconds), prefix="model_checkpoint")

# progress message to check time step size
progress_message(sim) =
        @printf("i: %04d, t: %s, Δt: %s, wall time: %s, Intertial Period %s\n",
        sim.model.clock.iteration, prettytime(sim.model.clock.time),
        prettytime(sim.Δt), prettytime((time_ns() - start_time) * 1e-9),sim.model.clock.time*fˢ/2π)

simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(1000)) # TimeInterval(0.5*(2*pi)/f) 

# diagnostic calculations, it is saved in 2 files with one saving the flow field and the other tke diagnostics
# calculate the pertubation in velocities
ua, va, wa = model.velocities # change back to ua, va, wa
um = Field(@at (Face, Center, Center) Average(ua, dims=1)) #averaging
vm = Field(@at (Center, Face, Center) Average(va, dims=1))
wm = Field(@at (Center, Center, Face) Average(wa, dims=1))
u = Field(@at (Center, Center, Center) ua - um) # calculating the Pertubations
v = Field(@at (Center, Center, Center) va - vm)
w = Field(@at (Center, Center, Center) wa - wm)
ub = model.background_fields.velocities.u
vb = model.background_fields.velocities.v
B = model.background_fields.tracers.b

ut = Field(@at (Center, Center, Center) ub+ua)
vt = Field(@at (Center, Center, Center) vb+va)

# buoyancy pertubation calculation
ba = model.tracers.b
bm = Field(@at (Center, Center, Center) Average(ba, dims=1))
b = Field(@at (Center, Center, Center) ba - bm)
bt = Field(@at (Center, Center, Center) B+ba)

PV = ErtelPotentialVorticity(model, ub, vb, 0, B, coriolis) # potential vorticity calculation


# output writers
output = (; u, ua, ub, v, va, vb, w, wa, b, ba, B, PV) # pertubation fields and PV

simulation.output_writers[:fields] = NetCDFOutputWriter(model, output;
                                                          schedule = TimeInterval(0.05*(2*pi)/fˢ),
                                                          filename = path_name*"flow_fields_height_"*string(hu)*"_theta_"*string(θ)*"_stratification_"*string(N²)*"_interior_velocity_"*string(V∞)*"_visc_"*string(ν1)*".nc",
                                                          overwrite_existing = true)


# With initial conditions set and an output writer at the ready, we run the simulation
# simulation.stop_time = 15*((2π)/f)seconds

run!(simulation) # , pickup=true

using NCDatasets, CairoMakie

xω, yω, zω = nodes(ωy)

ds = NCDataset(simulation.output_writers[:fields].filepath, "r")

fig = Figure(size = (800, 600))

axis_kwargs = (xlabel = "Across-slope distance (m)",
               ylabel = "Slope-normal\ndistance (m)",
               limits = ((0, Lx), (0, Lz)),
               )

ax_ω = Axis(fig[2, 1]; title = "Ertel PV", axis_kwargs...)

n = Observable(1)

ωy = @lift ds["PV"][:, :, $n]
hm_ω = heatmap!(ax_ω, xω, zω, ωy, colorrange = (-0.015, +0.015), colormap = :balance)
Colorbar(fig[2, 2], hm_ω; label = "s⁻¹")
ct_b = contour!(ax_ω, xb, zb, B, levels=-1e-3:0.5e-4:1e-3, color=:black)


times = collect(ds["time"])
title = @lift "t = " * string(prettytime(times[$n]))
fig[1, :] = Label(fig, title, fontsize=20, tellwidth=false)

fig