import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as d3
# import logging
# logger = logging.getLogger(__name__)


# Parameters
N_list = [(1e-5)**(0.5)]#np.linspace((1e-7)**(0.5),(8e-4)**(0.5),21) #np.array([(1e-5)**(0.5)]) # stratification
delta_list = [0] # np.linspace(0, 1, 11)
theta = 0.4225*(1e-4)/((1e-5)**(0.5))
f = 10**(-4)

phi = np.pi/2
uz = np.sin(phi)
vz = np.cos(phi)-1
bz = np.cos(phi)-1

# Basis
coord = d3.Coordinate('z')
dist = d3.Distributor(coord, dtype=np.complex128)
basis = d3.Legendre(coord, 56, bounds=(0, 1),dealias=3/2)

# Fields
u = dist.Field(bases=basis)
v = dist.Field(bases=basis)
w = dist.Field(bases=basis)
b = dist.Field(bases=basis)
p = dist.Field(bases=basis)
tau_1 = dist.Field(name="tau_1")
tau_2 = dist.Field(name="tau_2")
tau_3 = dist.Field(name="tau_3")
tau_4 = dist.Field(name="tau_4")
omega = dist.Field()

# Substitutions
z = dist.local_grid(basis)
delta = dist.Field()
N = dist.Field()
# S = dist.Field()
gamma = dist.Field()
Gshear = dist.Field()
n = dist.Field()
Ri = dist.Field()
alpha = dist.Field()
k = dist.Field()
one_z = dist.Field(bases=basis)
one_z['g'] = 1-z
lift_basis = basis.derivative_basis(1)
lift = lambda A: d3.Lift(A,lift_basis,-1)
dz = lambda A: d3.Differentiate(A, coord)
wz = dz(w)+lift(tau_1)+lift(tau_2)
pz = dz(p)+lift(tau_3)+lift(tau_4)

# Problem
problem = d3.EVP([u,v,w,b,p,tau_1,tau_2,tau_3,tau_4], eigenvalue=omega, namespace=locals())

problem.add_equation("-1j*omega*u-gamma**(0.5)*delta*uz*w+1j*k*gamma**(0.5)*delta*uz*one_z*u+1j*k*Ri*p-alpha*b= 0")
problem.add_equation("-1j*omega*v-(1+delta+delta*gamma*vz)*w+1j*k*gamma**(0.5)*delta*uz*one_z*v+u-n*theta*w= 0")
problem.add_equation("-1j*omega*n**2*w+1j*k*gamma**(0.5)*delta*n*uz*one_z*w+n*theta*v+Ri*pz-Ri*b= 0")
problem.add_equation("1j*k*u+wz=0")
problem.add_equation("-1j*omega*b + u*theta/(n*(1-gamma)) + (1-gamma)**(-1)*(1-gamma-gamma*delta*bz*theta)*w + 1j*k*gamma**(0.5)*delta*uz*one_z*b = 0")

# problem.add_equation("-1j*omega*u-gamma**(0.5)*delta*uz*w+1j*k*gamma**(0.5)*delta*uz*one_z*u+1j*k*Ri*p-alpha*b = 0")
# problem.add_equation("-1j*omega*v-(1+delta+delta*gamma*vz)*w+1j*k*gamma**(0.5)*delta*uz*one_z*v+u-n*theta*w= 0")
# problem.add_equation("-1j*omega*n**2*w+1j*k*gamma**(0.5)*delta*n*uz*one_z*w+n*theta*v+Ri*pz-Ri*b = 0")
# problem.add_equation("1j*k*u+wz=0")
# problem.add_equation("-1j*omega*b + u*alpha/Ri + (1-gamma-gamma*delta*theta*bz)*w + 1j*k*gamma**(0.5)*delta*uz*one_z*b = 0")
problem.add_equation("w(z=0)=0")
problem.add_equation("w(z=1)=0")
problem.add_equation("p(z=0)=0")
problem.add_equation("p(z=1)=0")

def gamma_lower_limit(N,delta):
    S2 = N**2*theta**2/f**2
    return (3-S2)/(3*(1+S2)-delta*S2)

def gamma_upper_limit(N,delta):
    S2 = N**2*theta**2/f**2
    return 1.00001#2*(1+S2)**(-1)

# Solver
solver = problem.build_solver()
evals = []
gammas = []
k_list = np.arange(1,100,1)
for Ni in N_list:
    gammas2 = []
    N['g'] = Ni
    eval2 = []
    for deltai in delta_list:
        delta['g'] = deltai
        gamma_list = [gamma_upper_limit(Ni,deltai)]#np.linspace(gamma_lower_limit(Ni,deltai),gamma_upper_limit(Ni,deltai),5)
        gammas2.append(gamma_list)
        eval3 = []
        for gammai in gamma_list:
            eval4 = []
            gamma['g'] = gammai
            Gsheari = (theta*(Ni)**2*gammai)/f
            Gshear['g'] = Gsheari
            n['g'] = f/(Gsheari)
            Ri['g'] = Ni**2*(1-gammai)/((Gsheari)**2)
            alpha['g'] = (Ni**2*theta*(1-gammai))/(f*Gsheari)
            Hi = 1 # (f*0.05)/(gammai*Ni**2*theta)
            for ki in k_list:
                k['g'] = ki*np.pi*(Gsheari*Hi)/f
                solver.solve_dense(solver.subproblems[0], rebuild_matrices=True)
                sorted_evals = np.sort(solver.eigenvalues.imag)
                eval4.append(sorted_evals[:1])
            eval3.append(eval4)
        eval2.append(eval3)
    evals.append(eval2)
    gammas.append(gammas2)
evals = np.array(evals)

import matplotlib.colors as colors

# fig, ax = plt.subplots()
# cs = ax.contourf(delta_list,N_list*theta/f,np.log10(evals[:,:,0,0]*(-1)),levels=20,cmap="bwr")
# cbar = fig.colorbar(cs,ticks=[-2,-1,0,1,2,3,4,5,6,7,8,9])
# ax.set_xlabel(r"$\delta$")
# ax.set_ylabel(r"$S_\infty$")
# # ax.title(r"")
# cbar.ax.set_ylabel("log of growth rate")
# plt.show()
fig, ax = plt.subplots()
ax.set_xlabel("Wavenumber")
ax.set_ylabel("Growth Rate")
# ax.title("Growth Rate for delta = 0.5 S=0.158 and var. k")
ax.plot(k_list,evals[0,0,0,:]*(-1))
plt.show()
