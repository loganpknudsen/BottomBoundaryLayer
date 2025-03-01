### Load in Packages
using Pkg
Pkg.activate(".")
Pkg.instantiate()
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

path_name = "/glade/derecho/scratch/knudsenl/data/" #args["path"]


# grid specifications
arch = has_cuda_gpu() ? GPU() : CPU()
@info("Arch => $arch")

Lx = 2000meters
Lz = 200meters
Nx = 512 # 512 originally
Nz = 128 # # 128 originally Note to self, maintain 2 to 1 resolution ration

grid = RectilinearGrid(arch; topology = (Periodic, Flat, Bounded),
                       size = (Nx, Nz),
                       x = (0, Lx),
                       z = (0,Lz))


# tilted domain parameters
θ = 5*10^(-3) # degrees 10^(-2) is previous value for 110 meter layer
ĝ = [sind(θ), 0, cosd(θ)] # gravity vector

# realistic mid latitude for now
buoyancy = Buoyancy(model = BuoyancyTracer(), gravity_unit_vector = -ĝ)
coriolis = ConstantCartesianCoriolis(f = 1e-4, rotation_axis = ĝ)

# parameters for simulation
const V∞ = 0.05 # m s⁻¹ interior velocity
const f = 1e-4 # coriolis parameter
const N² = 1e-5 # interior stratification
const S∞ = (N²*θ^2)/(f^2) # slope burger number
const γ = (1+S∞)^(-1) # 0 PV parameter
const hu = ceil((f*V∞)/(γ*N²*θ)) # Height of Boundary Layer
const fˢ=(f^2+θ^2*N²)^(0.5) # modified oscillation
const δ = 0.5
const vₒ = γ*(N²*θ)/(f)*δ # Initial v shear perturbation
# a1-h1 are constants for the following oscillations, calculate here for efficiency
const a1 = (f*vₒ)/(fˢ) 
const b1 = (f^2*vₒ)/(fˢ)^2
const d1 = ((fˢ^2-f^2)*vₒ)/(fˢ)^2
const e1 = N²*θ*(f*vₒ)/(fˢ)^2

# array of paramerers for background function
p =(; N², θ, f, V∞, hu, γ, vₒ, fˢ, a1, b1, d1, e1)

# heaviside function for boundary layer
@inline heaviside(x,z) = ifelse(z < 0, zero(z), one(z))

# oscillation functions for background
@inline sn_fn(x,z,t,p) = sin(p.fˢ*t)
@inline cs_fn(x,z,t,p) = cos(p.fˢ*t)

u_pert(x,z,t,p) = p.a1*sn_fn(x,z,t,p) # shear
v_pert(x,z,t,p) = p.b1*cs_fn(x,z,t,p) + p.d1
b_pert(x,z,t,p) = p.e1*(cs_fn(x,z,t,p) - 1)

u_adjustment(x, z, t, p) = u_pert(x,z,t,p)*(p.hu-z)*heaviside(x,p.hu-z)
v_adjustment(x, z, t, p) = p.V∞ - p.γ*(p.θ * p.N²)/(p.f)*(p.hu-z)*heaviside(x,p.hu-z) + v_pert(x,z,t,p)*(p.hu-z)*heaviside(x,p.hu-z)
constant_stratification(x, z, t, p) = p.N²*x*p.θ + p.N²*z + p.N²*p.γ*(p.hu-z)*heaviside(x,p.hu-z) + b_pert(x,z,t,p)*(p.hu-z)*heaviside(x,p.hu-z)

U_field = BackgroundField(u_adjustment, parameters=p)
V_field = BackgroundField(v_adjustment, parameters=p)
B_field = BackgroundField(constant_stratification, parameters=p)

# Boundary condition set up
# Free Slip Boundary Conditions

b_bc_top= FluxBoundaryCondition(-1*N²)
# b_bc_bottom= GradientBoundaryCondition(N²)
buoyancy_grad = FieldBoundaryConditions(top=b_bc_top) # ,bottom=b_bc_bottom

# diffusitivity and viscosity values for closure
const ν1 = 1e-6
closure = ScalarDiffusivity(ν=ν1, κ=ν1)

start_time = time_ns()

# model set up 
model = NonhydrostaticModel(; grid, buoyancy, coriolis, closure,
                            timestepper = :RungeKutta3,
                            advection = WENO(),
                            tracers = :b,
                            boundary_conditions = (; b=buoyancy_grad),
                            background_fields = (; u=U_field, v=V_field, b=B