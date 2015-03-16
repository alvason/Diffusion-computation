# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Diffusion computation
# https://github.com/alvason/diffusion-computation
# 
# ### Lecture002 Numerical solution for the diffusion equation

# <codecell>

'''
author: Alvason Zhenhua Li
date:   03/15/2015
'''

%matplotlib inline

import numpy as np
import matplotlib.pyplot as plt
import time
import IPython.display as idisplay
from mpl_toolkits.mplot3d.axes3d import Axes3D

AlvaFontSize = 23;
numberingFig = 0;

numberingFig = numberingFig + 1;
plt.figure(numberingFig, figsize=(12,3))
plt.axis('off')
plt.title(r'$ Diffusion \ equation \ and \ analytic \ solution $',fontsize = AlvaFontSize)
plt.text(0,1.0/2,r'$ \frac{\partial H(x,t)}{\partial t}=\xi \ \frac{\partial^2 H(x,t)}{\partial x^2} $', fontsize = 1.5*AlvaFontSize)
plt.text(0,0.0/2,r'$H(x,t) = \frac{1}{(1 + 4 \ \xi \ t)^{1/2}} e^\frac{-x^2}{1 + 4 \ \xi \ t}}$', fontsize = 1.5*AlvaFontSize)
plt.show()

# <codecell>

# define GridXX function for making 2D-grid from 1D-grid
def AlvaGridXX(gridX, totalGPoint_Y):
    gridXX = gridX;
    for n in range(totalGPoint_Y - 1):
        gridXX = np.vstack((gridXX, gridX));
    return gridXX;

# <codecell>

# Initial conditions
minX = float(0); maxX = float(3);
minT = float(0); maxT = float(1);

totalGPoint_X = int(20 + 1);
dx = (maxX - minX)/(totalGPoint_X - 1);
gridX = np.linspace(minX, maxX, totalGPoint_X); 

totalGPoint_T = int(100 + 1); 
dt = (maxT - minT)/(totalGPoint_T - 1);
gridT = np.linspace(minT, maxT, totalGPoint_T)

gridHtx = np.zeros([totalGPoint_T, totalGPoint_X])

movingRate = 1.0/10 # diffusion coefficience

tn = 0; # inital time = minT = gridT[tn = 0]
for xn in range(totalGPoint_X):
    gridHtx[tn, xn] = (1.0/(1.0+4.0*movingRate*gridT[tn]))*np.exp(-(gridX[xn]-(maxX-minX)/2.0)**2/(1.0+4.0*movingRate*gridT[tn]))

initialH = gridHtx.copy();

numberingFig = numberingFig + 1;
plt.figure(numberingFig,figsize=(10,5));     
plt.plot(gridX[:], gridHtx[:,:].T);
plt.grid(True)
plt.title(r'$ Initial \ conditions $', fontsize = AlvaFontSize);
plt.xlabel(r'$x \ (space)$', fontsize = AlvaFontSize); plt.ylabel(r'$H(x,t)$', fontsize = AlvaFontSize)
plt.text(maxX, 2.0/3, r'$H(t,x) = \frac{1}{(1 + 4 \ \xi \ t)^{1/2}} e^\frac{-x^2}{1 + 4 \ \xi \ t}}$', fontsize = AlvaFontSize);
plt.text(maxX, 1.0/3, r'$ dt = %f $'%(dt), fontsize = AlvaFontSize);
plt.text(maxX, minX, r'$ dx = %f $'%(dx), fontsize = AlvaFontSize); 
plt.show()

# <codecell>

# Analytic solution
gridHtx_A = np.zeros([totalGPoint_T, totalGPoint_X]); # Define the space for analytic values

for tn in range(totalGPoint_T):  
    for xn in range(totalGPoint_X):
        gridHtx_A[tn,xn] = (1.0/(1.0 + 4.0*movingRate*gridT[tn]))*np.exp(-(gridX[xn]-(maxX-minX)/2.0)**2/(1.0+4.0*movingRate*gridT[tn]))


numberingFig = numberingFig + 1;
plt.figure(numberingFig,figsize=(10,6));     
plt.plot(gridX[:], gridHtx_A[:,:].T);
plt.grid(True)
plt.title(r'$Analytic \ solution: (dt = %f,\ dx = %f) $'%(dt, dx), fontsize = AlvaFontSize);
plt.xlabel(r'$x \ (space)$', fontsize = AlvaFontSize); plt.ylabel(r'$H(x,t)$', fontsize = AlvaFontSize)
plt.text(maxX, 2.0/3, r'$ \frac{\partial H(x,t)}{\partial t}=\xi \ \frac{\partial^2 H(x,t)}{\partial x^2} $', fontsize = 1.5*AlvaFontSize)
plt.text(maxX, 1.0/3, r'$H(t,x) = \frac{1}{(1 + 4 \ \xi \ t)^{1/2}} e^\frac{-x^2}{1 + 4 \ \xi \ t}}$', fontsize = 1.5*AlvaFontSize);
plt.show()

# for 3D plotting
X = AlvaGridXX(gridX, totalGPoint_T); 
Y = AlvaGridXX(gridT, totalGPoint_X).T; 
Z = gridHtx_A;

numberingFig = numberingFig + 1;
plt.figure(numberingFig,figsize=(12,6)); 
plt.pcolor(X, Y, Z);
plt.title(r'$ Analytic \ solution: (dt = %f,\ dx = %f) $'%(dt, dx), fontsize = AlvaFontSize);
plt.xlabel(r'$x \ (space)$', fontsize = AlvaFontSize); plt.ylabel(r'$H(x,t)$', fontsize = AlvaFontSize);
plt.colorbar();
plt.show()

# <codecell>

numberingFig = numberingFig + 1;
plt.figure(numberingFig, figsize=(12,3))
plt.axis('off')
plt.title('Centered fomular of 2nd derivative', fontsize = AlvaFontSize)
plt.text(0, 1.0/2, r'$ \frac{\partial^2 H(t,x)}{\partial x^2} \approx \frac{H(t,x - \Delta x) - 2H(t,x) + H(t,x + \Delta x)}{(\Delta x)^2} $', fontsize = 1.5*AlvaFontSize)
plt.text(0, 0, r'$ \Longrightarrow\Delta H(t+\Delta t,x) = \Delta t (\frac{H(t,x - \Delta x) - 2H(t,x) + H(t,x + \Delta x)}{(\Delta x)^2}) $', fontsize = 1.5*AlvaFontSize)
plt.show()

# <codecell>

# Numerical solution with traditional loop
gridHtx = initialH.copy();
#print gridHtx
start_time = time.time();
for tn in range(totalGPoint_T - 1):
    for xn in range(totalGPoint_X):
        if (xn - 1) < 0: leftX = 0.0
        else: leftX = gridHtx[tn, xn - 1];
        if (xn + 1) > totalGPoint_X - 1: rightX = 0.0
        else: rightX = gridHtx[tn, xn + 1];    
        gridHtx[tn + 1, xn] = gridHtx[tn , xn] + dt*(leftX - 2.0*gridHtx[tn, xn] + rightX)/(dx)**2;
#       print 'tn = %i, xn = %i, gridHtx = %f'% (tn, xn, gridHtx[tn + 1,xn]);
            
stop_time = time.time(); 
total_time = stop_time - start_time;
print 'total computational time = %f'% (total_time);
#print gridHtx  
numberingFig = numberingFig + 1;
plt.figure(numberingFig,figsize=(10,5));     
plt.plot(gridX, gridHtx.T);
plt.grid(True)
plt.title(r'$ Numerical \ solution \ (dt = %f,\ dx = %f) $'%(dt, dx), fontsize = AlvaFontSize);
plt.xlabel(r'$x \ (space)$', fontsize = AlvaFontSize); plt.ylabel(r'$H(x,t)$', fontsize = AlvaFontSize);
plt.show()

# <codecell>

# Numerical solution with the feature of NumPy
gridHtx = initialH.copy();
Hgear = initialH.copy();
#print gridHtx

start_time = time.time();
for tn in range(totalGPoint_T - 1):
    leftX =   np.roll(Hgear[tn, :], 1); leftX[0:1] = 0.0; 
    centerX = Hgear[tn, :]; 
    rightX =  np.roll(Hgear[tn, :], -1); rightX[-1:] = 0.0;
    
    gridHtx[tn + 1, :] = Hgear[tn, :] + dt*(leftX - 2.0*centerX + rightX)/(dx)**2; 
    Hgear = gridHtx.copy()

stop_time = time.time(); 
total_time = stop_time - start_time;
print 'total computational time = %f'% (total_time);

#print gridHtx
numberingFig = numberingFig + 1;
plt.figure(numberingFig,figsize=(10,5));     
plt.plot(gridX, gridHtx.T);
plt.grid(True)
plt.title(r'$ Numerical \ solution: (dt = %f,\ dx = %f) $'%(dt, dx), fontsize = AlvaFontSize);
plt.xlabel(r'$x \ (space)$', fontsize = AlvaFontSize); plt.ylabel(r'$H(x,t)$', fontsize = AlvaFontSize);
plt.show()

numberingFig = numberingFig + 1;
plt.figure(numberingFig,figsize=(12,6)); 
plt.pcolor(X, Y, gridHtx);
plt.title(r'$ Analytic \ solution: (dt = %f,\ dx = %f) $'%(dt, dx), fontsize = AlvaFontSize);
plt.xlabel(r'$x \ (space)$', fontsize = AlvaFontSize); plt.ylabel(r'$H(x,t)$', fontsize = AlvaFontSize);
plt.colorbar();
plt.show()

# <codecell>


