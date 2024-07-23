using Pkg
Pkg.rm("CUDA_Runtime_jll")
Pkg.add("CUDA_Runtime_jll")
Pkg.instantiate()
using CUDA

