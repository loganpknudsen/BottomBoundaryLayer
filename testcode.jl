using Pkg
Pkg.activate(".")
Pkg.add("CUDA")
Pkg.add("Oceananigans")
Pkg.add("Oceanostics")
Pkg.add("Random")
Pkg.add("Printf")
Pkg.add("ArgParse")
# using Oceananigans
# using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
# using Oceananigans.AbstractOperations: @at, Average
# using Oceananigans.Grids: Center, Face
# using Oceananigans.Units
# using Random
# using Printf
# using ArgParse
# using CUDA: has_cuda_gpu
# versioninfo()
using CUDA 
# using Oceanostics

