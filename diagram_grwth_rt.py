import numpy as np
import scipy.integrate as sc
import matplotlib.pyplot as plt
import matplotlib.colors as cl
N2 = 1e-5
f = 1e-4
delta = 0.5

def m_t(t,m,beta):
    return m-delta*beta**(-2)*(np.cos(beta*t)-1)

def PSI_system(v,t,m,alpha,Ri,gamma,beta,theta,n):
    psi, D = v
    dpsidt = (alpha*m_t(t,m,beta)-Ri*np.cos(theta))*(1-Ri**(-1)*gamma**(-1)*m)*D+(n*np.sin(theta)+np.cos(theta)*m_t(t,m,beta))*(1-m-delta-n*np.sin(theta))*D
    dDdt = (n**2+m_t(t,m,beta)**2)**(-1)*psi
    return [dpsidt, dDdt]

def max_growth_rate(gamma,Sinf):
    theta =np.arctan(Sinf*f/N2**(0.5))
    S2 = Sinf**2
    beta = (1+S2)**(0.5)
    fstar = f*np.cos(theta)*beta
    delta = 0.5
    lmbd = N2*np.tan(theta)*gamma/f
    n = f/lmbd
    Ri = N2*(1-gamma)/lmbd**2
    alpha = N2*(1-gamma)*np.tan(theta)/(f*lmbd)
    m = np.linspace(0,1.5,75)
    frequencies = []
    growth_rates = []
    for i in m:
        t = np.linspace(0,2*np.pi/beta,100)
        sol1 = sc.odeint(PSI_system,[1,0],t,args=(i,alpha,Ri,gamma,beta,theta,n,)) 
        sol2 = sc.odeint(PSI_system,[0,1],t,args=(i,alpha,Ri,gamma,beta,theta,n,)) 
        M = np.array([[sol1[-1,0],sol1[-1,1]],[sol2[-1,0],sol2[-1,1]]])
        eigs = np.log(np.linalg.eig(M)[0]+0*1j)/(2*np.pi/beta)
        frequencies.append(eigs.imag)
        growth_rates.append(eigs.real)
    frequencies = np.array(frequencies)
    growth_rates = np.array(growth_rates)
    gr = growth_rates[np.argmax(growth_rates[:,0]),:]
    fr = frequencies[np.argmax(growth_rates[:,0]),:]
    m_max = m[np.argmax(growth_rates[:,0])]
    return gr, fr, m_max

Sinf = np.linspace(0.001,2.001,20)
gamma = np.linspace(0.001,1.101,20)
gr_list = []
fr_list = []
mm_list = []
for i in Sinf:
    gr_list2 =[]
    fr_list2 =[]
    mm_list2 =[]
    print(i)
    for j in gamma:
        sigma, fr, m_max = max_growth_rate(j,i)
        gr_list2.append(sigma)
        fr_list2.append(fr)
        mm_list2.append(m_max)
    gr_list.append(gr_list2)
    fr_list.append(fr_list2)
    mm_list.append(mm_list2)

sigmas = np.array(gr_list)
freqs = np.array(fr_list)
m_nums = np.array(mm_list)

plt.contourf(Sinf,gamma,np.abs(sigmas[:,:,0]).T,norm=cl.LogNorm(vmax=1))
plt.colorbar()
plt.show()