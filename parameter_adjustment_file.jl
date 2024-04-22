using Oceananigans
using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.Grids: Center, Face
using Oceananigans.Units
using Random
using Printf
using ArgParse

# grid = RectilinearGrid(arch; size=(1024, 200), y=(0,3000),z=(-200,0), topology=(Flat, Periodic, Bounded))

# tilted domain parameters
θ = 10^(-2) # degrees 
# ĝ = [θ, 0, 1] # gravity vector small angle
ĝ = [sind(θ), 0, cosd(θ)] # gravity vector

# realustic mid latitude for now
coriolis = ConstantCartesianCoriolis(f = 1e-4, rotation_axis = ĝ)

# parameters
V∞ = 0.05 # m s⁻¹
f = coriolis.fz
N² = 5e-6 # interior stratification
ϕ = 0
S∞ = (N²*θ^2)/(f^2)
γ = (1+S∞)^(-1)#(θ^2+1)*(1+S∞*(θ^2+1))^(-1)
hu = (f*V∞)/(γ*N²*θ) # set to negative
print(hu)
print((f^2)/(N²*(1-γ)^2))
