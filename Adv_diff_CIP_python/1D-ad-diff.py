# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 16:13:50 2019

@author: Wentao
"""

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# Setting the river geometry and model parameters

chlen  = 100       # length of river reach (m)

nx = 50    # number of grid
dx = chlen/float(nx)    # size of grid
dt = 1.    # computational time step (sec)
u  = 0.5    # advection velocity
D  = 0.1   # diffusion coefficient
h0 = 1.

i_advection = 1   # 0: upwind, 1: CIP

dx2 = dx*dx
dx3 = dx2*dx

x   = np.zeros(nx+1)
h   = np.zeros(nx+1)
wh  = np.zeros(nx+1)
gh  = np.zeros(nx+1)
wgh = np.zeros(nx+1)

# initial condition for h

x[0] = 0.

for i in np.arange(1,nx+1):
    x[i] = x[i-1]+dx
    
for i in np.arange(nx+1):
    if 0.7*nx < i < 0.9*nx:
        h[i] = h0
    else:
        h[i]  = 0.5*h0
        
for i in np.arange(1,nx):
    gh[i] = (h[i+1]-h[i-1])*0.5/dx

time = 0.
optime = 0.
etime = 50.
tuk = 2.

n = 0
    
while time<etime:
    
    for i in np.arange(1,nx):
        wh[i] = h[i]+D*(h[i+1]-2.0*h[i]+h[i-1])/dx2*dt
        
    wh[0] = wh[1]
    wh[nx] = wh[nx-1]
    
    for i in np.arange(1,nx):
        wgh[i] = gh[i]+0.5*(wh[i+1]-h[i+1]-wh[i-1]+h[i-1])/dx
        
    wgh[0] = wgh[1]
    wgh[nx] = wgh[nx-1]
    
    if i_advection==0:
        for i in np.arange(1,nx):
            h[i] = wh[i]-((u+abs(u))*(wh[i]-wh[i-1])+(u-abs(u))*(wh[i+1]-wh[i]))*dt/dx*0.5
        
    else:
    
        for i in np.arange(1,nx):
            xx = -u*dt
            
            if u>0:
                isn = 1.
            else:
                isn = -1.
                
            im1 = i-int(isn)
            
            a1 = ((wgh[im1]+wgh[i])*dx*float(isn)-2.0*(wh[i]-wh[im1]))/(dx3*float(isn))
            b1 = (3.0*(wh[im1]-wh[i])+(wgh[im1]+2.0*wgh[i])*dx*float(isn))/dx2
            h[i]	= ((a1*xx+b1)*xx+wgh[i])*xx+wh[i]
            gh[i]	= (3.0*a1*xx+2.0*b1)*xx+wgh[i]
        
        
    h[0] = h[1]
    h[nx] = h[nx-1]
    gh[0] = gh[1]
    gh[nx] = gh[nx-1]
    
    if optime>tuk:
        print("Time= ", time, "sec")
        
        plt.xlim([0,chlen])
        plt.ylim([0,h0*1.2])
        plt.xlabel("Downstream distance, x (m)")
        plt.ylabel("Water depth, h (m)")
        plt.plot(x,h,color='k')
        
        nnn = str(n)
        
        plt.savefig('Figure' + nnn.zfill(4) +".jpg", dpi=300)
        
        plt.close()
        
        n += 1
        
        optime = optime-tuk
    
    optime+=dt
    time+=dt
    
