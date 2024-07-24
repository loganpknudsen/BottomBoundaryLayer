# using Pkg
# Pkg.activate(".")
# Pkg.add("Oceananigans")
# Pkg.add("CUDA")
# Pkg.add("Random")
# Pkg.add("Printf")
# Pkg.add("ArgParse")
# Pkg.add("Oceanostics")
# Pkg.instantiate()

using Oceananigans
using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.AbstractOperations: @at, Average
using Oceananigans.Grids: Center, Face
using Oceananigans.Units
using Random
using Printf
using ArgParse
using CUDA: has_cuda_gpu
using CUDA 
using Oceanostics