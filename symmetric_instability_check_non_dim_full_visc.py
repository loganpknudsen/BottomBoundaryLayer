# Libraries being imported
import numpy as np
import dedalus.public as d3
import xarray as xr

# Parameters
N_list = [(10**(-5))**(0.5)] # stratification
f = 10**(-4) # coriolis
theta = np.arctan((0.4225*f**(2)/(N_list[0]**(2)))**(0.5)) # angle of slope for problem
gm = 0.9 # stratification modification parameter
H = 1 # non-dimensionalized height
H_1 = f*0.05/(gm*N_list[0]**2*np.sin(theta)) # dimensional height
visc= 10**(-4) # viscosity/diffusitivity for this problem

# Basis
coord = d3.Coordinate('z')
dist = d3.Distributor(coord, dtype=np.complex128)
basis = d3.ChebyshevT(coord, 128, bounds=(0, H),dealias=3/2)

# Fields
u = dist.Field(name="u",bases=basis) # u-velocity
uz = dist.Field(name="uz",bases=basis) # w-velocity
uzz = dist.Field(name="uzz",bases=basis) # w-velocity
v = dist.Field(name="v",bases=basis) # v-velocity
vz = dist.Field(name="vz",bases=basis) # w-velocity
vzz = dist.Field(name="vzz",bases=basis) # w-velocity
w = dist.Field(name="w",bases=basis) # w-velocity
wz = dist.Field(name="wz",bases=basis) # w-velocity
wzz = dist.Field(name="wzz",bases=basis) # w-velocity
b = dist.Field(name="b",bases=basis) # buoyancy
bz = dist.Field(name="bz",bases=basis) # w-velocity
bzz = dist.Field(name="bzz",bases=basis) # w-velocity
p = dist.Field(name="p",bases=basis) # pressure

# Tau Fields for Boundary Conditions
tau_1 = dist.Field(name="tau_1")
tau_2 = dist.Field(name="tau_2")
tau_3 = dist.Field(name="tau_3")
tau_4 = dist.Field(name="tau_4")
tau_5 = dist.Field(name="tau_5")
tau_6 = dist.Field(name="tau_6")
tau_7 = dist.Field(name="tau_7")
tau_8 = dist.Field(name="tau_8")
tau_p = dist.Field(name="tau_p") # pressure gauge

# Eigenvalues for problem
omega = dist.Field(name="omega") # Frequency/Growth Rate, what we are solving for
k = dist.Field(name="k") # horizontal wavenumber, which is set in loop

# Parameters
Ri = dist.Field(name="Ri") # Richardson Number
Ek = dist.Field(name="Ek") # Ekman Number
n = dist.Field(name="n") # non-hydrostatic parameter
alpha = dist.Field(name="alpha") # slope parameter

# Build Lift basis and vertical derivatives for problem
lift_basis = basis.derivative_basis(1)
lift = lambda A: d3.Lift(A,lift_basis,-1)
dz = lambda A: d3.Differentiate(A, coord)
# uz = dz(u)+lift(tau_1)
# uzz = dz(uz)
# vz = dz(v)+lift(tau_3) 
# vzz = dz(vz)
# wz = dz(w)+lift(tau_5)
# wzz = dz(wz)
# pz = dz(p)
# bz = dz(b)+lift(tau_7)
# bzz = dz(bz)

# Substituion for Frequency and Wavenumber (for x-derivative only, constant in y-direction assumed)
dx = lambda A: 1j*k*A
dt = lambda A: -1j*omega*A

# Problem
problem = d3.EVP([u,uz,uzz,v,vz,vzz,w,wz,wzz,b,bz,bzz,p,tau_1,tau_2,tau_3,tau_4,tau_5,tau_6,tau_7,tau_8,tau_p], eigenvalue=omega, namespace=locals())

# Add Equations to Eigenvalue Problem
problem.add_equation("dt(u)-v*np.cos(theta)+Ri*dx(p)-alpha*b*np.cos(theta)-Ek*uzz+lift(tau_2)= 0")
problem.add_equation("dt(v)+w+u*np.cos(theta)-np.sin(theta)*n*w-Ek*vzz+lift(tau_4)= 0")
problem.add_equation("n**2*dt(w)+n*np.sin(theta)*v+Ri*dz(p)-Ri*b*np.cos(theta)-n**(2)*Ek*wzz+lift(tau_6)= 0")
problem.add_equation("dx(u)+wz+lift(tau_p)=0")
problem.add_equation("dt(b)+ Ri**(-1)*(1+alpha)*u*np.cos(theta)+(1-Ri**(-1)*n*np.tan(theta))*w*np.cos(theta)-Ek*bzz+lift(tau_8)= 0") 
problem.add_equation("uz- dz(u)-lift(tau_1)=0")
problem.add_equation("uzz- dz(uz)=0")
problem.add_equation("vz- dz(v)-lift(tau_3)=0")
problem.add_equation("vzz- dz(vz)=0")
problem.add_equation("wz- dz(w)-lift(tau_5)=0")
problem.add_equation("wzz- dz(wz)=0")
problem.add_equation("bz- dz(b)-lift(tau_7)=0")
problem.add_equation("bzz- dz(bz)=0")
# Setting Boundary Values
problem.add_equation("u(z=0)=0")
problem.add_equation("u(z="+str(H)+")=0")
problem.add_equation("v(z=0)=0")
problem.add_equation("v(z="+str(H)+")=0")
problem.add_equation("w(z=0)=0")
problem.add_equation("w(z="+str(H)+")=0")
problem.add_equation("integ(p)=0")
problem.add_equation("b(z=0)=0")
problem.add_equation("b(z="+str(H)+")=0")

# Solver
solver = problem.build_solver()
evals = [] # list to save output
k_list = np.arange(0,21,5) # horizontal wavenumber values solver is run for
Ni = N_list[0] # short cut so stratification frequency does not need to be indexed everytime
Gsheari = (np.sin(theta)*(Ni)**2*(gm))/(f) # Geostrophic Shear
Ri['g'] = Ni**2*(1-gm)/(Gsheari**2) # Richardson number calculation
Ek['g'] = visc/(f*H_1**2) # Ekman number calculation
n['g'] = f/Gsheari # non-hydrostatic parameter calculation
alpha['g'] = (Ni**2*(1-gm)*np.tan(theta))/(f*Gsheari) # slope parameter calculation
for ki in k_list:
    k['g'] = ki # add current wavenumber to problem
    solver.solve_dense(solver.subproblems[0], rebuild_matrices=True) # run solver
    # access and store maximum growth rate
    omg = solver.eigenvalues
    omg[np.isnan(omg)] = 0.
    omg[np.isinf(omg)] = 0.
    idx = np.sort(omg.imag)
    sorted_evals = idx[-1:]
    evals.append(sorted_evals)

# Saving the data to an external .nc file
evals = np.array(evals)
gr_data = xr.Dataset(data_vars={"growth_rate":(["k"],evals[:,0])},coords={"k":k_list})

gr_data.to_netcdf("SI_non_dim_visc.nc")

