# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 17:00:06 2019

@author: Wentao
"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import csv
import os
import glob
import re


xmax  = 100       # length of river reach (m)
ymax  = 100

nx = 50    # number of grid
ny = 50
dx = xmax/float(nx)    # size of grid
dy = ymax/float(ny)    # size of grid
dt = 0.2    # computational time step (sec)
u  = 0.5    # advection velocity
v  = 0.5
D  = 0.1   # diffusion coefficient
h0 = 1.
i_cip = 1

dx2 = dx*dx
dx3 = dx2*dx
dy2 = dy*dy
dy3 = dy2*dy

x    = np.zeros((ny+1,nx+1))
y    = np.zeros((ny+1,nx+1))
h    = np.zeros((ny+1,nx+1))
wh   = np.zeros((ny+1,nx+1))
ghx  = np.zeros((ny+1,nx+1))
wghx = np.zeros((ny+1,nx+1))
ghy  = np.zeros((ny+1,nx+1))
wghy = np.zeros((ny+1,nx+1))

for j in np.arange(ny+1):
    x[j][0] = 0.
    
    for i in np.arange(1,nx+1):
        x[j][i] = x[j][i-1]+dx
        
for i in np.arange(nx+1):
    y[0][i] = 0.
    
    for j in np.arange(1,ny+1):
        y[j][i] = y[j-1][i]+dy
        
for j in np.arange(ny+1):
    for i in np.arange(nx+1):
        if 0.2*nx < i < 0.4*nx and 0.2*ny < j < 0.4*ny:
            h[j][i] = h0
        else:
            h[j][i]  = 0
            
for j in np.arange(1,ny):
    for i in np.arange(1,nx):
        ghx[j][i] = (h[j][i+1]-h[j][i-1])*0.5/dx
        ghy[j][i] = (h[j+1][i]-h[j-1][i])*0.5/dy
        
for j in np.arange(1,ny):
    ghx[j][ 0] = ghx[j][   1]
    ghx[j][ny] = ghx[j][ny-1]
    ghy[j][ 0] = ghy[j][   1]
    ghy[j][ny] = ghy[j][ny-1]
        
for i in np.arange(nx+1):
    ghx[ 0][i] = ghx[   1][i]
    ghx[ny][i] = ghx[ny-1][i]
    ghy[ 0][i] = ghy[   1][i]
    ghy[ny][i] = ghy[ny-1][i]
    
         
time = 0.
optime = 0.
etime = 150.
tuk = 2.

n = 0
    
while time<etime:
    
    for j in np.arange(1,ny):
        for i in np.arange(1,nx):
            wh[j][i] = h[j][i]+D*dt*( (h[j][i+1]-2.0*h[j][i]+h[j][i-1])/dx**2.0+(h[j+1][i]-2.0*h[j][i]+h[j-1][i])/dy**2.0 )

    for j in np.arange(1,ny):
        wh[j][ 0] = wh[j][   1]
        wh[j][ny] = wh[j][ny-1]
        
    for i in np.arange(nx+1):
        wh[ 0][i] = wh[   1][i]
        wh[ny][i] = wh[ny-1][i]
        
    for j in np.arange(1,ny):
        for i in np.arange(1,nx):
            wghx[j][i] = ghx[j][i]+(-wh[j][i-1]+wh[j][i+1]+h[j][i-1]-h[j][i+1])*0.5/dx
            wghy[j][i] = ghy[j][i]+(-wh[j-1][i]+wh[j+1][i]+h[j-1][i]-h[j+1][i])*0.5/dy
    
    for j in np.arange(1,ny):
        wghx[j][ 0] = wghx[j][   1]
        wghx[j][ny] = wghx[j][ny-1]
        wghy[j][ 0] = wghy[j][   1]
        wghy[j][ny] = wghy[j][ny-1]
        
    for i in np.arange(nx+1):
        wghx[ 0][i] = wghx[   1][i]
        wghx[ny][i] = wghx[ny-1][i]
        wghy[ 0][i] = wghy[   1][i]
        wghy[ny][i] = wghy[ny-1][i]        
    
    if i_cip==0:
    
        for j in np.arange(1,ny):
            for i in np.arange(1,nx):
                udhdx = ((u+abs(u))*(wh[j][i]-wh[j][i-1]) + (u-abs(u))*(wh[j][i+1]-wh[j][i]))*0.5/dx
                vdhdy = ((v+abs(v))*(wh[j][i]-wh[j-1][i]) + (v-abs(v))*(wh[j+1][i]-wh[j][i]))*0.5/dy
                
                h[j][i] = wh[j][i] - (udhdx+vdhdy)*dt
                
    else:
        
        for j in np.arange(1,ny):
            for i in np.arange(1,nx):
                
                xx = -u*dt
                yy = -v*dt
                isn = int(np.sign(u))
                jsn = int(np.sign(v))
                fis = float(isn)
                fjs = float(jsn)
                im1 = i-isn
                jm1 = j-jsn
                
                a1	= ((wghx[j][im1]+wghx[j][i])*dx*fis-2.0*(wh[j][i]-wh[j][im1]))/(dx3*fis)
                e1	= (3.0*(wh[j][im1]-wh[j][i])+(wghx[j][im1]+2.0*wghx[j][i])*dx*fis)/dx2
                b1	= ((wghy[jm1][i]+wghy[j][i])*dy*fjs-2.0*(wh[j][i]-wh[jm1][i]))/(dy3*fjs)
                f1	= (3.0*(wh[jm1][i]-wh[j][i])+(wghy[jm1][i]+2.0*wghy[j][i])*dy*fjs)/dy2
                tmp	= wh[j][i]-wh[jm1][i]-wh[j][im1]+wh[jm1][im1]
                tmq	= wghy[j][im1]-wghy[j][i]
                d1	= (-tmp-tmq*dy*fjs)/(dx*dy2*fis)
                c1	= (-tmp-(wghx[jm1][i]-wghx[j][i])*dx*fis)/(dx2*dy*fjs)
                g1	= (-tmq+c1*dx2)/(dx*fis)
				
                h[j][i]   = ((a1*xx+c1*yy+e1)*xx+g1*yy+wghx[j][i])*xx+((b1*yy+d1*xx+f1)*yy+wghy[j][i])*yy+wh[j][i]
                ghx[j][i] = (3.0*a1*xx+2.0*(c1*yy+e1))*xx+(d1*yy+g1)*yy+wghx[j][i]
                ghy[j][i] = (3.0*b1*yy+2.0*(d1*xx+f1))*yy+(c1*xx+g1)*xx+wghy[j][i]
           
            
    for j in np.arange(1,ny):
        h[j][ 0] = h[j][   1]
        h[j][ny] = h[j][ny-1]
        
    for i in np.arange(nx+1):
        h[ 0][i] = h[   1][i]
        h[ny][i] = h[ny-1][i]

    for j in np.arange(1,ny):
        ghx[j][ 0] = ghx[j][   1]
        ghx[j][ny] = ghx[j][ny-1]
        ghy[j][ 0] = ghy[j][   1]
        ghy[j][ny] = ghy[j][ny-1]
        
    for i in np.arange(nx+1):
        ghx[ 0][i] = ghx[   1][i]
        ghx[ny][i] = ghx[ny-1][i]
        ghy[ 0][i] = ghy[   1][i]
        ghy[ny][i] = ghy[ny-1][i]
        
    if optime>tuk:
        print("Time= ", time, "sec")
        
        plt.figure(figsize=(4,4*ymax/xmax*0.9))
        plt.xlabel('x(m)')
        plt.ylabel('y(m)')
        plt.axis([0,xmax,0.,ymax])
        
        levels = np.linspace(0,h0,21)
        labels = np.linspace(0,h0,6)
                
        plt.contourf(x, y, h, levels, cmap=cm.rainbow, extend="both")
        plt.colorbar(ticks=labels, label='h (m)')
        
        nnn = str(n)
                
        plt.savefig('Figure' + nnn.zfill(4) +".jpg", dpi=300)
        plt.close()
        
        n += 1
        
        optime = optime-tuk
    
    optime+=dt
    time+=dt
          

            
            
