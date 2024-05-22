using Oceananigans

# tilted domain parameters
θ = 10^(-1) # degrees 
# ĝ = [θ, 0, 1] # gravity vector small angle
ĝ = [sind(θ), 0, cosd(θ)] # gravity vector

# realustic mid latitude for now
coriolis = ConstantCartesianCoriolis(f = 1e-4, rotation_axis = ĝ)

# parameters
V∞ = 0.1 # m s⁻¹
f = coriolis.fz
N² =1e-5 #1e-4 # interior stratification
ϕ = 0
S∞ = (N²*θ^2)/(f^2)
γ = (1+S∞)^(-1)#(θ^2+1)*(1+S∞*(θ^2+1))^(-1)
hu = ceil(Int,((f*V∞)/(γ*N²*θ))) # set to negative
fˢ=(f^2+θ^2*N²)^(0.5)
uₒ = 0#γ*(N²*θ)/(f)*cos(ϕ)
vₒ = γ*(N²*θ)/(f)*0.1#*sin(ϕ)
bₒ = vₒ*((θ*N²)/(f)) # initial stratification
q = vₒ*(θ*N²)
Ri = (N²*(1-γ))/((N²^2*γ^2*θ^2)/f^2)
RiPSI = (4)/(3*γ)-(4*f*vₒ)/(3*N²*θ*γ^2)+(5*f^2*bₒ)/(3*N²^2*θ^2*γ^2)-(bₒ-f*vₒ)/(3*N²*θ*γ)
RiPSI2 =  4/(3*γ)+(f^2*bₒ)/(N²^2*θ^2*γ^2)+(f*vₒ)/(N²*θ*γ^2)*(f^2/(3*fˢ^2)-4/3)+(f^2*bₒ)/(3*fˢ^2*N²*θ*γ^2)
Term1 = (4)/(3*γ)
Term2 = (4*f*vₒ)/(3*N²*θ*γ^2)
Term3 = (5*f^2*bₒ)/(3*N²^2*θ^2*γ^2)
Term4 = (bₒ-f*vₒ)/(3*N²*θ*γ)
ν = (V∞*hu)/(1.3*10^6)
println("$(hu)")
# println("$((f^2)/(N²*(1-γ)^2))")
# println("$(bₒ)")
# println("v_o $(vₒ)")
# println("$(q)")
# println("$(γ)")
println("$(Ri)")
v = 1+S∞
println("$(v)")
q = (1+(f*vₒ)/(N²*θ)+bₒ/(N²*θ))/((4*f*vₒ)/(N²*θ)-(f^2*bₒ)/(N²^2*θ^2))
println("$(RiPSI)")
println("RiPsi 2 $(RiPSI2)")
println("q $(q)")
println("T1 $(Term1)")
println("T2 $(Term2)")
println("T3 $(Term3)")
println("T4 $(Term4)")
# println("$(ν)")
