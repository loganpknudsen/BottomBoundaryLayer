import numpy as np
import dedalus.public as d3
import xarray as xr



# Parameters
f = 10**(-4)
N_list = [81.7*f]  # np.linspace((1e-7)**(0.5),(8e-4)**(0.5),51) #np.array([(1e-5)**(0.5)]) # stratification
delta_list = [0.6] #np.linspace(0, 1, 26)
S = 7.7*f
H = 1 
nz = 256

# Basis
coord = d3.Coordinate('z')
dist = d3.Distributor(coord, dtype=np.complex128)
basis = d3.Chebyshev(coord, nz, bounds=(0, H))

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
z = dist.local_grid(basis) #dist.local_grid(basis)
one_z = dist.Field(name="one_z",bases=basis)
one_z['g'] = z
delta = dist.Field()
N = dist.Field()
Gshear = dist.Field()
Ri = dist.Field()
n = dist.Field()
t =dist.Field()
u_sz = np.cos(t)
v_sz = np.sin(t)
b_sz = 1-np.cos(t)
k = dist.Field()
dy = lambda A: 1j*k*A
dt = lambda A: -1j*omega*A
lift_basis = basis.derivative_basis(1)
lift = lambda A: d3.Lift(A,lift_basis,-1)
dz = lambda A: d3.Differentiate(A, coord)
wz = dz(w)#+lift(tau_1)+lift(tau_2)
pz = dz(p)# +lift(tau_3)+lift(tau_4)


# Problem
problem = d3.EVP([u,v,w,b,p], eigenvalue=omega, namespace=locals()) # ,tau_1,tau_2,tau_3,tau_4

problem.add_equation("dt(u)+(1+delta*u_sz)*w-delta*v_sz*one_z*dy(u)-v= 0")
problem.add_equation("dt(v)-delta*v_sz*w-delta*v_sz*one_z*dy(v)+u+Ri*dy(p)=0")
problem.add_equation("n**2*dt(w)-n**2*delta*v_sz*one_z*dy(w)+Ri*pz-Ri*b=0") 
problem.add_equation("dy(v)+wz=0")
problem.add_equation("dt(b)+(1-delta*Ri**(-1)*b_sz)*w-Ri**(-1)*v-delta*v_sz*one_z*dy(b)=0") # *gamma**(-1)
# problem.add_equation("w(z=0)=0")
# problem.add_equation("w(z="+str(H)+")=0")
# problem.add_equation("p(z=0)=0")
# problem.add_equation("p(z="+str(H)+")=0")


# Solver
solver = problem.build_solver()
evals_r = []
evals_i =[]
gammas = []
k_list = np.arange(0.01,30.1,10)
time = np.linspace(0,(2*np.pi),3) #np.arange(0,(2*np.pi+1)/(1+N_list[0]**2*theta**2*f**(-2))**(0.5),1*(1+N_list[0]**2*theta**2*f**(-2))**(-0.5)) # np.arange(0,2*np.pi,0.1)
us = []
usc = []
vs = []
vsc = []
ws = []
wsc = []
bs = []
bsc = []
xs = []
for ti in time:
    gammas5 = []
    evals5 = []
    evals_i1 = []
    ut = []
    vt = []
    wt = []
    bt = []
    xt = []
    utc = []
    vtc = []
    wtc = []
    btc = []
    t["g"]= ti
    for Ni in N_list:
        gammas2 = []
        N['g'] = Ni
        eval2 = []
        eval_i2 = []
        for deltai in delta_list:
            delta['g'] = deltai
            eval3 = []
            eval_i3 = []
            gammas3 = []
            Gsheari = (S**2)/(f)
            Rii = Ni**2/(Gsheari**2)
            Ri['g'] = Rii
            ni = f/Gsheari
            n['g'] = ni
            Gshear['g'] = Gsheari
            for ki in k_list:
                ui = []
                vi = []
                wi = []
                bi = []
                xi = []
                uic = []
                vic = []
                wic = []
                bic = []
                k['g'] = ki
                solver.solve_dense(solver.subproblems[0], rebuild_matrices=True)
                omg = solver.eigenvalues
                omg[np.isnan(omg)] = 0.
                omg[np.isinf(omg)] = 0.
                # omg[omg>0] = 0 #  TO AVOID SPURIOUS VALUES....
                idx = np.argsort(omg.imag)
                sorted_evals = solver.eigenvalues[idx[-1]].imag
                sorted_evals_i = solver.eigenvalues[idx[-1]].real
                solver.set_state(idx[-1])
                x_domain = np.linspace(0,2*np.pi/ki,nz)
                xt.append(x_domain)


                # print(np.shape(solver.state))
                ui = np.real(u['g'].reshape(nz,1)@np.exp(1j*ki*x_domain.reshape(1,nz)))
                u_max = np.max(ui)
                ui = ui/u_max
                ut.append(ui)
                vi = np.real(v['g'].reshape(nz,1)@np.exp(1j*ki*x_domain.reshape(1,nz)))
                vi = vi/u_max
                vt.append(vi)
                wi = np.real(w['g'].reshape(nz,1)@np.exp(1j*ki*x_domain.reshape(1,nz)))
                wi = wi/u_max
                wt.append(wi)
                bi = np.real(b['g'].reshape(nz,1)@np.exp(1j*ki*x_domain.reshape(1,nz)))
                b_max = np.max(bi)
                bi = bi/b_max
                bt.append(bi)
                eval3.append(sorted_evals)
                eval_i3.append(sorted_evals_i)
            eval2.append(eval3)
            eval_i2.append(eval_i3)
        evals5.append(eval2)
        evals_i1.append(eval_i2)
    us.append(ut)
    vs.append(vt)
    ws.append(wt)
    bs.append(bt)
    xs.append(xt)
    usc.append(utc)
    vsc.append(vtc)
    wsc.append(wtc)
    bsc.append(btc)
    evals_r.append(evals5)
    evals_i.append(evals_i1)
    
evals_r = np.array(evals_r)
evals_i = np.array(evals_i)
us = np.array(us)
vs = np.array(vs)
ws = np.array(ws)
bs = np.array(bs)
xs = np.array(xs)
# g_index= np.linspace(0,len(gamma_list)+1,len(gamma_list))
gr_data = xr.Dataset(data_vars={"growth_rate":(["t","N","delta","k"],evals_r[:,:,:,:]),"oscillation":(["t","N","delta","k"],evals_i[:,:,:,:])},coords={"t":time,"N":N_list,"delta":delta_list,"k":k_list})
gr_data.to_netcdf("PSI_non_dim_full_form_mid_res.nc") 
grid_normal = basis.global_grid(dist,scale=1).ravel()
field_data = xr.Dataset({"u_structure":(["t","k","z","x"],us[:,:,:,:]),"v_structure":(["t","k","z","x"],vs[:,:,:,:]),"w_structure":(["t","k","z","x"],ws[:,:,:,:]),"b_structure":(["t","k","z","x"],bs[:,:,:,:]),"x_domain":(["t","k","z"],xs)},coords={"t":time,"k":k_list,"z":grid_normal,"x":np.linspace(0,1,nz)})
field_data.to_netcdf("PSI_non_dim_field_mid_res.nc")


