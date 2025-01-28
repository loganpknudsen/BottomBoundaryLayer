import numpy as np
# import matplotlib.pyplot as plt
import dedalus.public as d3
import xarray as xr
# import logging
# logger = logging.getLogger(__name__)


# Parameters
N_list = [(0.4225*10**(-8)/(5e-3)**2)**(0.5)] #np.linspace((1e-7)**(0.5),(8e-4)**(0.5),51) #np.array([(1e-5)**(0.5)]) # stratification
delta_list = [0] #np.linspace(0, 1, 26)
theta =  5e-3 #(0.4225*10**(-8)/(10**(-5)))**(0.5)
f = 10**(-4)
gm = 0.9
H = f*0.01/((1-gm)*N_list[0]**2*theta)
lmbd = N_list[0]**2*theta*(1-gm)/f
print(H)
