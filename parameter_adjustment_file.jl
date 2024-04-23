using Oceananigans

# tilted domain parameters
θ = 10^(-2) # degrees 
# ĝ = [θ, 0, 1] # gravity vector small angle
ĝ = [sind(θ), 0, cosd(θ)] # gravity vector

# realustic mid latitude for now
coriolis = ConstantCartesianCoriolis(f = 1e-4, rotation_axis = ĝ)

# parameters
V∞ = 0.5 # m s⁻¹
f = coriolis.fz
N² = 1e-4 # interior stratification
ϕ = 0
S∞ = (N²*θ^2)/(f^2)
γ = (1+S∞)^(-1)#(θ^2+1)*(1+S∞*(θ^2+1))^(-1)
hu = (f*V∞)/(γ*N²*θ) # set to negative
fˢ=(f^2+θ^2*N²)^(0.5)
uₒ = 0#γ*(N²*θ)/(f)*cos(ϕ)
vₒ = γ*(N²*θ)/(f)*0.1#*sin(ϕ)
bₒ = vₒ*((θ*N²)/(f)) # initial stratification
q = vₒ*(θ*N²)
println("$(hu)")
println("$((f^2)/(N²*(1-γ)^2))")
println("$(bₒ)")
println("$(q)")
