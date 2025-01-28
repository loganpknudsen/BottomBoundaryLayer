import numpy as np
# import matplotlib.pyplot as plt
import dedalus.public as d3
import xarray as xr
import logging
logger = logging.getLogger(__name__)


# Parameters
N_list = [(1e-5)**(0.5)] #np.linspace((1e-7)**(0.5),(8e-4)**(0.5),51) #np.array([(1e-5)**(0.5)]) # stratification
delta_list = [0] #np.linspace(0, 1, 26)
f = 10**(-4) #*(1+N_list[0]**2*theta**2/10**(-8))**(0.5)
gm = 0.9
S2 = 0.4225
theta = np.arctan((S2)**(0.5)*10**(-4)/(N_list[0]))
H = 1
H_1 = f*0.05/(gm*N_list[0]**2*np.sin(theta))
lmbd = N_list[0]**2*np.sin(theta)*gm/f

# Basis
coord = d3.Coordinate('z')
dist = d3.Distributor(coord, dtype=np.complex128)
basis = d3.Chebyshev(coord, 64, bounds=(0, H))

# Fields
u = dist.Field(name="u",bases=basis)
v = dist.Field(name="v",bases=basis)
w = dist.Field(name="w",bases=basis)
b = dist.Field(name="b",bases=basis)
p = dist.Field(name="p",bases=basis)
tau_1 = dist.Field(name="tau_1")
tau_2 = dist.Field(name="tau_2")
tau_3 = dist.Field(name="tau_3")
tau_4 = dist.Field(name="tau_4")
omega = dist.Field(name="omega")

# Substitutions
# z = dist.local_grid(basis)
delta = dist.Field()
gamma = dist.Field()
Gshear = dist.Field()
Ri = dist.Field()
n = dist.Field()
alpha = dist.Field()
k = dist.Field()
lift_basis = basis.derivative_basis(1)
lift = lambda A: d3.Lift(A,lift_basis,-1)
dz = lambda A: d3.Differentiate(A, coord)
wz = dz(w) + lift(tau_1) + lift(tau_2)
pz = dz(p) + lift(tau_3) + lift(tau_4)
dx = lambda A: 1j*k*A
dt = lambda A: -1j*omega*A


# Problem
problem = d3.EVP([u,v,w,b,p,tau_1,tau_2,tau_3,tau_4], eigenvalue=omega, namespace=locals()) # 

problem.add_equation("dt(u)-v+Ri*dx(p)= 0")
problem.add_equation("dt(v)+w+(1-n**(-1)*np.tan(theta))*u= 0") 
problem.add_equation("pz-b= 0") # -1j*omega*n**2*w+n*theta*v+Ri*pz-Ri*b= 0
problem.add_equation("dx(u)+wz=0")
problem.add_equation("dt(b) -np.tan(theta)*n**(-1)*u+w= 0") 
problem.add_equation("w(z=0)=0")
problem.add_equation("w(z="+str(H)+")=0")
problem.add_equation("p(z=0)=0")
problem.add_equation("p(z="+str(H)+")=0")


# # Solver
# solver = problem.build_solver()
solver = problem.build_solver() # entry_cutoff=0


evals = []
gammas = []
k_list = np.arange(0,41,1)
time = [0] # np.arange(0,2*np.pi,0.1)
k_l = []
for ti in time:
    gammas5 = []
    evals5 = []
    for Ni in N_list:
        gammas2 = []
        eval2 = []
        for deltai in delta_list:
            delta['g'] = deltai
            gamma_list = [gm]#np.linspace(gamma_lower_limit(Ni,deltai),gamma_upper_limit(Ni,deltai),21)
            eval3 = []
            gammas3 = []
            for gammai in gamma_list:
                eval4 = []
                gammas4 = []
                gamma['g'] = gammai
                Gsheari = -1*(theta*(Ni)**2*(gammai))/(f)
                Rii = Ni**2*(1-gammai)/(Gsheari**2)
                Ri['g'] = Rii
                ni = f/Gsheari
                n['g'] = ni
                alphai = (Ni**2*(1-gammai)*theta)/(f*Gsheari)
                alpha['g'] = alphai
                Gshear['g'] = Gsheari
                for ki in k_list:
                    k['g'] = ki
                    k_l.append(ki)
                    solver.solve_dense(solver.subproblems[0],rebuild_matrices=True)
                    # solver.solve_dense(solver.subproblems[0], rebuild_matrices=True)
                    omg = solver.eigenvalues
                    omg[np.isnan(omg)] = 0.
                    omg[np.isinf(omg)] = 0.
                    # omg[omg>0] = 0 #  TO AVOID SPURIOUS VALUES....
                    idx = np.sort(omg.imag)
                    sorted_evals = 1*idx[-1:]
                    eval4.append(sorted_evals)
                    gammas4.append([gammai])
                eval3.append(eval4)
                gammas3.append(gammas4)
            eval2.append(eval3)
            gammas2.append(gammas3)
        evals5.append(eval2)
        gammas5.append(gammas2)
    evals.append(evals5)
    gammas.append(gammas5)
    
evals = np.array(evals)
gammas = np.array(gammas)
print(np.shape(evals))
print(np.shape(gammas))
g_index= np.linspace(0,len(gamma_list)+1,len(gamma_list))
gr_data = xr.Dataset(data_vars={"growth_rate":(["t","N","delta","gamma_index","k"],evals[:,:,:,:,:,0]),"gamma":(["t","N","delta","gamma_index","k"],gammas[:,:,:,:,:,0])},coords={"t":time,"N":N_list,"delta":delta_list,"gamma_index":g_index,"k":k_l})
gr_data.to_netcdf("SI_non_dim.nc")

