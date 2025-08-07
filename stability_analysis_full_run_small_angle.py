import numpy as np
import scipy.integrate as sc
import xarray as xr

def m_t(t,m,fstar,delta,lmbd):
    phi = np.pi/2 
    return m+delta*lmbd/(fstar)*np.sin(t+phi)

def PSI_system(v,t,m,theta,gm,S2,delta,lmbd,fstar):
    psi, D = v
    A0 = -S2/(theta*(1+S2))*(m+(1-gm)/theta)
    A1 = m+(S2)*(theta*(1+S2))**(-1)
    dpsidt = (A0-A1*m_t(t,m,fstar,delta,lmbd))*D
    dDdt = (1+m_t(t,m,fstar,delta,lmbd)**2)**(-1)*psi
    return [dpsidt, dDdt]

max_grs = []
max_ms = []
max_frs = []
gms = []
dS = 0.5 #0.01
theta = 0.01
S_list = np.arange(dS,2+dS,dS)
f = 1e-4

tau = 2*np.pi
dt = 200
t = np.linspace(0, tau+1/dt, dt)
dm = 1
m = np.arange(-40, 40+dm, dm)
dgm = 250
ddelta = 250
delta_list = np.linspace(0,1+1/ddelta,ddelta)
for S in S_list:
    gms_2 = []
    max_grs_sub2 = []
    max_ms_sub2 = []
    max_frs_sub2 = []
    for i in delta_list:
        max_grs_sub = []
        max_ms_sub = []
        max_frs_sub = []
        S2 = S**2
        N2 = S2*f**2/theta**2
        beta = (1+S2)**(0.5)
        fstar = f*beta
        gml = (1+(1-2)*S2)/(1+S2) 
        gmu = (1+(1-4/3)*S2)/(1+S2)
        gm_list = np.linspace(gml,gmu+1/dgm,dgm)
        gms_2.append(gm_list)
        for gm in gm_list:
            lmbd = N2*theta*gm/f
            frequencies = []
            growth_rates = []
            for q in m:
                sol1 = sc.odeint(PSI_system, [1, 0], t, args=(q, theta, gm, S2, i, lmbd,fstar))
                sol2 = sc.odeint(PSI_system, [0, 1], t, args=(q, theta, gm, S2, i,lmbd,fstar)) 
                M = np.array([[sol1[-1, 0], sol1[-1, 1]], [sol2[-1, 0], sol2[-1, 1]]])
                eigs = np.log(np.linalg.eig(M)[0]+0*1j)/(tau)
                growth_rates.append(eigs.real)
                frequencies.append(eigs.imag)
            frequencies = np.array(frequencies)
            growth_rates = np.array(growth_rates)
            gr = growth_rates[np.argmax(np.abs(growth_rates[:,0])),:]
            fr = frequencies[np.argmax(np.abs(growth_rates[:,0])),:]
            m_max = m[np.argmax(np.abs(growth_rates[:,0]))]
            max_grs_sub.append(gr)
            max_ms_sub.append(m_max)
            max_frs_sub.append(fr)
        max_grs_sub2.append(max_grs_sub)
        max_ms_sub2.append(max_ms_sub)
        max_frs_sub2.append(max_frs_sub)
    max_grs.append(max_grs_sub2)
    max_ms.append(max_ms_sub2)
    max_frs.append(max_frs_sub2)
    gms.append(gms_2)
gms = np.array(gms) 
max_gr = np.array(max_grs)
max_ms = np.array(max_ms)
max_fr = np.array(max_frs)

output_file = xr.Dataset({"growth_rate":(["slope_burger_number","delta","strat_index",],np.abs(max_gr[:,:,:,0])),
            "frequency":(["slope_burger_number","delta","strat_index"],np.abs(max_fr[:,:,:,0])),
            "slope_angle":(["slope_burger_number","delta","strat_index"],max_ms),
           "strat_values":(["slope_burger_number","delta","strat_index"],gms)},
           coords = {"slope_burger_number":S_list,"delta":delta_list,"strat_index":np.linspace(0,1+1/dgm,dgm)})

output_file.to_netcdf("stability_analysis_output_small_angle_test_first_third.nc")