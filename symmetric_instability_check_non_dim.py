import numpy as np
import dedalus.public as d3
import matplotlib.pyplot as plt
import xarray as xr
import logging
logger = logging.getLogger(__name__)


# Parameters
N_list = [(1e-5)**(0.5)] #np.linspace((1e-7)**(0.5),(8e-4)**(0.5),51) #np.array([(1e-5)**(0.5)]) # stratification
delta_list = [0] #np.linspace(0, 1, 26)
f = 10**(-4) #*(1+N_list[0]**2*theta**2/10**(-8))**(0.5)
gm = 0.81
S2 = 0.25
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
# wz = dist.Field(name="wz",bases=basis)
b = dist.Field(name="b",bases=basis)
p = dist.Field(name="p",bases=basis)
# pz = dist.Field(name="pz",bases=basis)
tau_1 = dist.Field(name="tau_1")
tau_2 = dist.Field(name="tau_2")
tau_3 = dist.Field(name="tau_3")
tau_4 = dist.Field(name="tau_4")
# tau_5 = dist.Field(name="tau_5")
# tau_6 = dist.Field(name="tau_6")
# tau_p = dist.Field(name='tau_p')
omega = dist.Field(name="omega")

# Substitutions
# z = dist.local_grid(basis)
# hz = dist.Field(bases=basis)
# hz['g'] = z
gamma = dist.Field(name="gamma")
Gshear = dist.Field(name="Ghsear")
Ri = dist.Field(name="Ri")
n = dist.Field(name="n")
alpha = dist.Field(name="alpha")
k = dist.Field(name="k")
lift_basis = basis.derivative_basis(1)
lift = lambda A: d3.Lift(A,lift_basis,-1)
dx = lambda A: 1j*k*A
dt = lambda A: -1j*omega*A
dz = lambda A: d3.Differentiate(A, coord)
wz = dz(w) + lift(tau_1) + lift(tau_2)
pz = dz(p) + lift(tau_3) + lift(tau_4)

# Problem
problem = d3.EVP([u,v,w,b,p,tau_1,tau_2,tau_3,tau_4], eigenvalue=omega, namespace=locals()) #,tau_2
problem.add_equation("dt(u)-v*np.cos(theta)+Ri*dx(p)-alpha*b*np.cos(theta)= 0")
problem.add_equation("dt(v)+w+u*np.cos(theta)-np.sin(theta)*n*w= 0") 
problem.add_equation("pz-b*np.cos(theta)= 0") # n**2*dt(w)+n*np.sin(theta)*v+Ri*pz-Ri*b*np.cos(theta)
problem.add_equation("dx(u)+wz =0")
problem.add_equation("dt(b)+ Ri**(-1)*(1+alpha)*u*np.cos(theta)+(1-Ri**(-1)*n*np.tan(theta))*w*np.cos(theta)= 0") 
problem.add_equation("w(z='left')=0")
problem.add_equation("w(z='right')=0")
# problem.add_equation("integ(p)=0")
problem.add_equation("p(z='left')=0")
problem.add_equation("p(z='right')=0")
# problem.add_equation("b(z='left')=0")
# problem.add_equation("b(z='right')=0")

# Solver

evals = []
gammas = []
us = []
vs = []
ws = []
bs = []
k_list = np.arange(0,51,1)
time = [0] #np.arange(0,2*np.pi,0.1)
k_l = []
for ti in time:
    gammas5 = []
    evals5 = []
    # u5 = []
    # v5 =[]
    # w5 = []
    # b5 = []
    for Ni in N_list:
        gammas2 = []
        eval2 = []
        # u2 = []
        # v2 =[]
        # w2 = []
        # b2 = []
        for deltai in delta_list:
            gamma_list = [gm]#np.linspace(gamma_lower_limit(Ni,deltai),gamma_upper_limit(Ni,deltai),21)
            eval3 = []
            gammas3 = []
            # u3 = []
            # v3 =[]
            # w3 = []
            # b3 = []
            for gammai in gamma_list:
                eval4 = []
                gammas4 = []
                # u4 =[]
                # v4 =[]
                # w4 =[]
                # b4 =[]
                gamma['g'] = gammai
                Nii = (Ni)**2#*(1-gammai)
                Gsheari = (np.sin(theta)*Nii*gammai)/(f)
                Rii = (Nii*(1-gammai))/(Gsheari**2) #*(1-gammai*Nii*theta**2/f**2)
                Ri['g'] = Rii
                ni = f/Gsheari
                n['g'] = ni
                alphai = (Nii*(1-gammai)*np.tan(theta))/(f*Gsheari) #*(1-gammai*Nii*theta**2/f**2)**(0.5)
                alpha['g'] = alphai
                Gshear['g'] = Gsheari
                for ki in k_list:
                    ui = []
                    solver = problem.build_solver() # 

                    k['g'] = ki
                    k_l.append(ki)
                    # solver.print_subproblem_ranks(solver.subproblems)
                    # solver.tol = 1-10
                    solver.solve_dense(solver.subproblems[0], rebuild_matrices=True) #
                    # solver.solve_sparse(solver.subproblems[0],10,1j,vO=1j, rebuild_matrices=True) #

                    # solver.iteration = 1000
                    
                    # solver.solve_dense(solver.subproblems[0], rebuild_matrices=True)
                    omg = solver.eigenvalues
                    omg[np.isnan(omg)] = 0.
                    omg[np.isinf(omg)] = 0.
                    # omg[omg>0] = 0 #  TO AVOID SPURIOUS VALUES....
                    idx = np.argsort(omg.imag)
                    sorted_evals = solver.eigenvalues[idx[-1]].imag
                    eval4.append([sorted_evals])
                    solver.set_state(idx[-1],solver.subsystems[0])

                    # print(np.shape(solver.state))
                    ui = v['g'].real
                    # print(solver.eigenvectors)

                    # print(np.shape(solver.eigenvectors))
                    vi = []#solver.state['v']
                    wi = [] #solver.state['w']
                    bi = [] #solver.state['b']
                    us.append(np.array(ui))
                    # v4.append(vi)
                    # w4.append(wi)
                    # b4.append(bi)
                    gammas4.append([gammai])
                eval3.append(eval4)
                # u3.append(u4)
                # v3.append(v4)
                # w3.append(w4)
                # b3.append(b4)
                gammas3.append(gammas4)
            eval2.append(eval3)
            # u2.append(u3)
            # v2.append(v3)
            # w2.append(w3)
            # b2.append(b3)
            gammas2.append(gammas3)
        evals5.append(eval2)
        # u5.append(u2)
        # v5.append(v2)
        # w5.append(w2)
        # b5.append(b2)
        gammas5.append(gammas2)
    # us.append(u5)
    # vs.append(v5)
    # ws.append(w5)
    # bs.append(b5)
    evals.append(evals5)
    gammas.append(gammas5)
us = np.array(us)    
vs = np.array(vs) 
ws = np.array(ws) 
bs = np.array(us) 
print(np.shape(us))
evals = np.array(evals)
gammas = np.array(gammas)
print(np.shape(evals))
print(np.shape(gammas))
g_index= np.linspace(0,len(gamma_list)+1,len(gamma_list))
gr_data = xr.Dataset(data_vars={"growth_rate":(["t","N","delta","gamma_index","k","gr_i"],evals[:,:,:,:,:,:]),"gamma":(["t","N","delta","gamma_index","k"],gammas[:,:,:,:,:,0])},coords={"t":time,"N":N_list,"delta":delta_list,"gamma_index":g_index,"k":k_l,"gr_i":np.linspace(0,1,np.shape(evals)[-1])})
gr_data.to_netcdf("SI_non_dim.nc")
print(us)
# field_data = xr.Dataset(us)
# field_data.to_netcdf("SI_non_dim_field.nc")
# grid_normal = basis.global_grid(dist,scale=1).ravel()
# fig, ax = plt.subplots()
# ax.plot( us[0,:]/np.max(us[-1,:]),grid_normal,c="c",linestyle="dashed",label="k=0")
# ax.plot( us[10,:]/np.max(us[-1,:]),grid_normal,c="b",linestyle="dashed",label="k=10")
# ax.plot( us[20,:]/np.max(us[-1,:]),grid_normal,c="r",linestyle="dashed",label="k=20")
# ax.plot( us[30,:]/np.max(us[-1,:]),grid_normal,c="g",linestyle="dashed",label="k=30")
# ax.plot( us[-1,:]/np.max(us[-1,:]),grid_normal,c="k",label="k=40")
# ax.legend()
# ax.set_title("Vertical Structure of Fastest Growing Mode")
# ax.set_ylabel("z")
# ax.set_xlabel("v-velocity normalized by maximum velocity of final mode")
# plt.show()