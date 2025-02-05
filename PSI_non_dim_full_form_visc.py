import numpy as np
import dedalus.public as d3
import xarray as xr



# Parameters
N_list = [(1*10**(-5))**(0.5)]  # np.linspace((1e-7)**(0.5),(8e-4)**(0.5),51) #np.array([(1e-5)**(0.5)]) # stratification
theta =  5*10**(-3) #np.arctan(1*10**(-4)/N_list[0])
delta_list = [0.5] #np.linspace(0, 1, 26)
f = 10**(-4)
S2 = N_list[0]**2*theta**2/f**2
gm = (1+S2)**(-1)
H = 1 
H_1 = f*0.05/(gm*N_list[0]**2*theta)
lmbd = N_list[0]**2*theta*(1-gm)/f
visc  = 10**(-4)

# Basis
coord = d3.Coordinate('z')
dist = d3.Distributor(coord, dtype=np.complex128)
basis = d3.Chebyshev(coord, 64, bounds=(0, H))

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
omega = dist.Field(name="omega")
Ek = dist.Field(name="Ek") # Ekman Number

# Substitutions
z = dist.local_grid(basis)
# onez = dist.Field(bases=basis)
# onez['g'] = z
one_z = dist.Field(name="one_z",bases=basis)
one_z['g'] = 1-z
beta = dist.Field()
delta = dist.Field()
N = dist.Field()
gamma = dist.Field()
Gshear = dist.Field()
Ri = dist.Field()
S = dist.Field()
n = dist.Field()
t =dist.Field()
u_sz = beta**(-1)*np.sin(beta*t)
v_sz = 1 + beta**(-2)*(np.cos(beta*t)-1)
b_sz = beta**(-2)*(np.cos(beta*t)-1)
alpha = dist.Field()
k = dist.Field()
dx = lambda A: 1j*k*A
dt = lambda A: -1j*omega*A
lift_basis = basis.derivative_basis(1)
lift = lambda A: d3.Lift(A,lift_basis,-1)
dz = lambda A: d3.Differentiate(A, coord)
# wz = dz(w)+lift(tau_1)+lift(tau_2)
# pz = dz(p)+lift(tau_3)+lift(tau_4)


# Problem
problem = d3.EVP([u,uz,uzz,v,vz,vzz,w,wz,wzz,b,bz,bzz,p,tau_1,tau_2,tau_3,tau_4,tau_5,tau_6,tau_p], eigenvalue=omega, namespace=locals()) # ,tau_p

problem.add_equation("dt(u)-delta*u_sz*w+delta*u_sz*one_z*dx(u)-v*np.cos(theta)+Ri*dx(p)-alpha*b*np.cos(theta)-Ek*uzz+lift(tau_2)= 0")
problem.add_equation("dt(v)+(1-delta*v_sz)*w+delta*u_sz*one_z*dx(v)+u*np.cos(theta)-n*np.sin(theta)*w-Ek*vzz+lift(tau_4)=0")
problem.add_equation("n**2*dt(w)+n**2*delta*u_sz*one_z*dx(w)+n*np.sin(theta)*v+Ri*dz(p)+lift(tau_p)-Ri*b*np.cos(theta)-n**(2)*Ek*wzz+lift(tau_6)=0") #+lift(tau_p)
problem.add_equation("dx(u)+wz=0")
problem.add_equation("dt(b)+Ri**(-1)*(1+alpha)*u*np.cos(theta)+(1-delta*Ri**(-1)*gamma**(-1)*b_sz-Ri**(-1)*n*np.tan(theta))*w*np.cos(theta)+delta*u_sz*one_z*dx(b)-Ek*bzz=0") # +lift(tau_8)
problem.add_equation("uz- dz(u)-lift(tau_1)=0")
problem.add_equation("uzz- dz(uz)=0")
problem.add_equation("vz- dz(v)-lift(tau_3)=0")
problem.add_equation("vzz- dz(vz)=0")
problem.add_equation("wz- dz(w)-lift(tau_5)=0")
problem.add_equation("wzz- dz(wz)=0")
problem.add_equation("bz- dz(b)=0") # -lift(tau_7)
problem.add_equation("bzz- dz(bz)=0")
# Setting Boundary Values
problem.add_equation("u(z=0)=0")
problem.add_equation("u(z="+str(H)+")=0")
problem.add_equation("v(z=0)=0")
problem.add_equation("v(z="+str(H)+")=0")
problem.add_equation("w(z=0)=0")
problem.add_equation("w(z="+str(H)+")=0")
problem.add_equation("integ(p)=0")
# problem.add_equation("b(z=0)=0")
# problem.add_equation("b(z="+str(H)+")=0")


# Solver
solver = problem.build_solver()
evals_r = []
evals_i =[]
gammas = []
k_list = np.arange(0.1,32.1,4)
# phase = np.pi/2
time = np.linspace(0,(2*np.pi)*(1+N_list[0]**2*theta**2*f**(-2))**(-0.5),6) #np.arange(0,(2*np.pi+1)/(1+N_list[0]**2*theta**2*f**(-2))**(0.5),1*(1+N_list[0]**2*theta**2*f**(-2))**(-0.5)) # np.arange(0,2*np.pi,0.1)
us = []
usc = []
vs = []
vsc = []
ws = []
wsc = []
bs = []
bsc = []
for ti in time:
    gammas5 = []
    evals5 = []
    evals_i1 = []
    ut = []
    vt = []
    wt = []
    bt = []
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
            gamma_list = [gm]#np.linspace(gamma_lower_limit(Ni,deltai),gamma_upper_limit(Ni,deltai),21)
            eval3 = []
            eval_i3 = []
            gammas3 = []
            for gammai in gamma_list:
                eval4 = []
                eval_i4 = []
                gammas4 = []
                gamma['g'] = gammai
                Gsheari = (np.tan(theta)*(Ni)**2*(gammai))/(f)
                Rii = Ni**2*(1-gammai)/(Gsheari**2)
                Ek['g'] = visc/(f*H_1**2) # Ekman number calculation
                Ri['g'] = Rii
                S['g'] = (Ni*np.tan(theta))/f
                ni = f/Gsheari
                n['g'] = ni
                alphai = (Ni**2*(1-gammai)*np.tan(theta))/(f*Gsheari)
                alpha['g'] = alphai
                Gshear['g'] = Gsheari
                beta['g'] = (1+Ni**2*np.tan(theta)**2/f**2)**(0.5)
                for ki in k_list:
                    ui = []
                    vi = []
                    wi = []
                    bi = []
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

                    # print(np.shape(solver.state))
                    ui = (u['g']).real
                    ut.append(np.array(ui))
                    vi = (v['g']).real
                    vt.append(np.array(vi))
                    wi = (w['g']).real
                    wt.append(np.array(wi))
                    bi = (b['g']).real
                    bt.append(np.array(bi))
                    uic = (u['g']).imag
                    utc.append(np.array(uic))
                    vic = (v['g']).imag
                    vtc.append(np.array(vic))
                    wic = (w['g']).imag
                    wtc.append(np.array(wic))
                    bic = (b['g']).imag
                    btc.append(np.array(bic))
                    eval4.append([sorted_evals])
                    eval_i4.append([sorted_evals_i])
                    gammas4.append([gammai])
                eval3.append(eval4)
                eval_i3.append(eval_i4)
                gammas3.append(gammas4)
            eval2.append(eval3)
            eval_i2.append(eval_i3)
            gammas2.append(gammas3)
        evals5.append(eval2)
        evals_i1.append(eval_i2)
        gammas5.append(gammas2)
    us.append(ut)
    vs.append(vt)
    ws.append(wt)
    bs.append(bt)
    usc.append(utc)
    vsc.append(vtc)
    wsc.append(wtc)
    bsc.append(btc)
    evals_r.append(evals5)
    evals_i.append(evals_i1)
    gammas.append(gammas5)
    
evals_r = np.array(evals_r)
evals_i = np.array(evals_i)
gammas = np.array(gammas)
us = np.array(us)
vs = np.array(vs)
ws = np.array(ws)
bs = np.array(bs)
usc = np.array(usc)
vsc = np.array(vsc)
wsc = np.array(wsc)
bsc = np.array(bsc)
g_index= np.linspace(0,len(gamma_list)+1,len(gamma_list))
gr_data = xr.Dataset(data_vars={"growth_rate":(["t","N","delta","gamma_index","k"],evals_r[:,:,:,:,:,0]),"oscillation":(["t","N","delta","gamma_index","k"],evals_i[:,:,:,:,:,0]),"gamma":(["t","N","delta","gamma_index","k"],gammas[:,:,:,:,:,0])},coords={"t":time,"N":N_list,"delta":delta_list,"gamma_index":g_index,"k":k_list})
gr_data.to_netcdf("PSI_non_dim_full_form_visc_low_res_b.nc") 
grid_normal = basis.global_grid(dist,scale=1).ravel()
field_data = xr.Dataset({"u_structure":(["t","k","z"],us),"v_structure":(["t","k","z"],vs),"w_structure":(["t","k","z"],ws),"b_structure":(["t","k","z"],bs),"u_structure_complex":(["t","k","z"],usc),"v_structure_complex":(["t","k","z"],vsc),"w_structure_complex":(["t","k","z"],wsc),"b_structure_complex":(["t","k","z"],bsc)},coords={"t":time,"k":k_list,"z":grid_normal})
field_data.to_netcdf("PSI_non_dim_field_low_res_b.nc")


