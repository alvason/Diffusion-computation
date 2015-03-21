# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Diffusion computation
# https://github.com/alvason/diffusion-computation
# 
# ### Lecture002 --- Numerical solution for the diffusion equation

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
AlvaFigSize = (6, 6);
numberingFig = 0;

numberingFig = numberingFig + 1;
plt.figure(numberingFig, figsize=(12,3))
plt.axis('off')
plt.title(r'$ Diffusion \ equation \ and \ analytic \ solution $',fontsize = AlvaFontSize)
plt.text(0,1.0/2,r'$ \frac{\partial H(x,t)}{\partial t} = \
         \xi \frac{\partial^2 H(x,t)}{\partial x^2} $', fontsize = 1.5*AlvaFontSize)
plt.text(0,0.0/2,r'$H(x,t) = \frac{1}{(1 + 4 \xi \ t)^{1/2}} \
         e^\frac{-x^2}{1 + 4 \xi \ t}}$', fontsize = 1.5*AlvaFontSize)
plt.show()


# <codecell>

# define GridXX function for making 2D-grid from 1D-grid
def AlvaGridXX(gridX, totalGPoint_Y):
    gridXX = gridX;
    for n in range(totalGPoint_Y - 1):
        gridXX = np.vstack((gridXX, gridX));
    return gridXX;

# define analytic solution H(x,t)
def analyticHtx(t, x):
    analyticH = (1.0/np.sqrt(1.0 + 4.0*movingRate*t))*np.exp(-(x - (maxX - minX)/2.0 + 0
                                                               )**2/(1.0 + 4.0*movingRate*t));
    return analyticH

# <codecell>

# numerical way for 2nd order derivative
numberingFig = numberingFig + 1;
plt.figure(numberingFig, figsize=(12,3))
plt.axis('off')
plt.title('Centered fomular of 2nd derivative', fontsize = AlvaFontSize)
plt.text(0, 2.0/3, r'$ \frac{\partial^2 H(t,x)}{\partial x^2} \approx \
         \frac{H(t,x - \Delta x) - 2H(t,x) + H(t,x + \Delta x)}{(\Delta x)^2} $'
         , fontsize = AlvaFontSize)
plt.text(0, 1.0/3, r'$ \Longrightarrow\frac{\Delta H(t+\Delta t,x)}{\Delta t} \
         = \xi(\frac{H(t,x - \Delta x) - 2H(t,x) + H(t,x + \Delta x)}{(\Delta x)^2}) $'
         , fontsize = AlvaFontSize)
plt.text(0, 0, r'$ \Longrightarrow H(t+\Delta t,x) = H(t,x) + \Delta H(t+\Delta t,x) \
         = H(t,x) + \Delta t \ \xi(\frac{H(t,x - \Delta x) - 2H(t,x) \
         + H(t,x + \Delta x)}{(\Delta x)^2}) $', fontsize = AlvaFontSize)
plt.show()

# <codecell>

# Initial conditions
minX = float(0); maxX = float(40);
minT = float(0); maxT = float(40);

resolution = 40;

totalGPoint_X = int(resolution + 1);
dx = (maxX - minX)/(totalGPoint_X - 1);
gridX = np.linspace(minX, maxX, totalGPoint_X); 

totalGPoint_T = int(9*resolution + 1); 
dt = (maxT - minT)/(totalGPoint_T - 1);
gridT = np.linspace(minT, maxT, totalGPoint_T);

gridHtx = np.zeros([totalGPoint_T, totalGPoint_X]);

# for 3D plotting
X = AlvaGridXX(gridX, totalGPoint_T); 
Y = AlvaGridXX(gridT, totalGPoint_X).T; 

# diffusion coefficience (area moving rate per time)
movingRate = 1.0/10;

tn = 0; # inital time = minT = gridT[tn = 0]
for xn in range(totalGPoint_X):
    gridHtx[tn, xn] = analyticHtx(gridT[tn], gridX[xn]);
initialH = gridHtx.copy();

numberingFig = numberingFig + 1;
plt.figure(numberingFig, figsize = AlvaFigSize);     
plt.plot(gridX[:], gridHtx[:,:].T);
plt.grid(True)
plt.title(r'$ Initial \ conditions: (dt = %f,\ dx = %f) $'%(dt, dx), fontsize = AlvaFontSize);
plt.xlabel(r'$x \ (space)$', fontsize = AlvaFontSize); plt.ylabel(r'$H(x,t)$'
                                                                  , fontsize = AlvaFontSize)
plt.text(maxX, 2.0/3, r'$H(t,x) = \frac{1}{(1 + 4 \ \xi \ t)^{1/2}} \
         e^\frac{-x^2}{1 + 4 \ \xi \ t}}$', fontsize = AlvaFontSize);
plt.text(maxX, 1.0/3, r'$ dt = %f $'%(dt), fontsize = AlvaFontSize);
plt.text(maxX, minX, r'$ dx = %f $'%(dx), fontsize = AlvaFontSize); 
plt.show()


# Numerical solution with the feature of NumPy
gridHtx = initialH.copy();
Hcopy = initialH.copy();
#print gridHtx

start_time = time.time();

 

# 3 points scheme
# isolated boundary requires leftX[0] = centerX[0] ==> centered formula = (- centerX[0] + rightX[0])
# constant boundary requires leftX[0] = constant
for tn in range(totalGPoint_T - 1):
    centerX = Hcopy[tn, :]; 
    leftX = np.roll(Hcopy[tn, :], 1);
    rightX = np.roll(Hcopy[tn, :], -1);
    leftX[0] = centerX[0];
    rightX[-1] = centerX[-1]; 
    gridHtx[tn + 1, :] = Hcopy[tn, :] + dt*movingRate*(leftX - 2.0*centerX + rightX)/((dx)**2); 
    Hcopy = gridHtx.copy();
 
 
'''
# 5 points scheme    
for tn in range(totalGPoint_T - 1):
    leftX = np.roll(Hcopy[tn, :], 1); leftX[0:1] = 0.0; 
    left2X = np.roll(Hcopy[tn, :], 2); left2X[0:2] = 0.0;
    centerX = Hcopy[tn, :]; 
    rightX = np.roll(Hcopy[tn, :], -1); rightX[-1:] = 0.0;
    right2X = np.roll(Hcopy[tn, :], -2); right2X[-2:] = 0.0;
    
    gridHtx[tn + 1, 2:-2] = Hcopy[tn, 2:-2] + dt*movingRate*(-left2X[2:-2] + 16*leftX[2:-2] - 30*centerX[2:-2]
                                                       + 16*rightX[2:-2] - right2X[2:-2])/(12*(dx)**2);   
    Hcopy = gridHtx.copy();
'''

stop_time = time.time(); 
total_time = stop_time - start_time;
print 'total computational time = %f'% (total_time);

numberingFig = numberingFig + 1;
plt.figure(numberingFig, figsize = AlvaFigSize);     
plt.plot(gridX, gridHtx[0::resolution].T);
plt.grid(True)
plt.title(r'$ Numerical \ solution: (dt = %f,\ dx = %f) $'%(dt, dx), fontsize = AlvaFontSize);
plt.xlabel(r'$x \ (space)$', fontsize = AlvaFontSize);
plt.ylabel(r'$H(x,t)$', fontsize = AlvaFontSize);
plt.show()

numberingFig = numberingFig + 1;
plt.figure(numberingFig, figsize = AlvaFigSize);     
#plt.pcolormesh(X, Y, gridHtx); 
plt.contourf(X, Y, gridHtx, levels = np.arange(0,1,0.02));
plt.title(r'$ Numerical \ solution: (dt = %f,\ dx = %f) $'%(dt, dx), fontsize = AlvaFontSize);
plt.xlabel(r'$ x \ (space) $', fontsize = AlvaFontSize);
plt.ylabel(r'$ t \ (time) $', fontsize = AlvaFontSize);
plt.show()

# <codecell>

# Analytic solution
resolutionA = 500;
totalGPoint_TA = int(3*resolutionA + 1); totalGPoint_XA = int(3*resolutionA + 1);
gridX_A = np.linspace(minX, maxX, totalGPoint_XA);
gridT_A = np.linspace(minT, maxT, totalGPoint_TA);
gridHtx_A = np.zeros([totalGPoint_TA, totalGPoint_XA]); # Define the space for analytic values

for tn in range(totalGPoint_TA):  
    for xn in range(totalGPoint_XA):
        gridHtx_A[tn,xn] = analyticHtx(gridT_A[tn], gridX_A[xn])

numberingFig = numberingFig + 1;
plt.figure(numberingFig, figsize = AlvaFigSize);   
plt.plot(gridX_A[:], gridHtx_A[0::resolutionA].T);
plt.grid(True)
plt.title(r'$Analytic \ solution: (dt = %f,\ dx = %f) $'%(dt, dx), fontsize = AlvaFontSize);
plt.xlabel(r'$x \ (space)$', fontsize = AlvaFontSize);
plt.ylabel(r'$H(x,t)$', fontsize = AlvaFontSize)
plt.text(maxX, 2.0/3, r'$ \frac{\partial H(x,t)}{\partial t}=\xi \
         \frac{\partial^2 H(x,t)}{\partial x^2} $', fontsize = 1.5*AlvaFontSize)
plt.text(maxX, 1.0/3, r'$H(t,x) = \frac{1}{(1 + 4 \ \xi \ t)^{1/2}} \
         e^\frac{-x^2}{1 + 4 \xi \ t}}$', fontsize = 1.5*AlvaFontSize);
plt.show();


numberingFig = numberingFig + 1;
plt.figure(numberingFig, figsize = (20,12));   
plt.plot(gridX_A[:], gridHtx_A[0::resolutionA].T)
plt.plot(gridX[:], gridHtx[0::resolution*3].T, linestyle = 'dashed');
plt.show()

# <codecell>

dx

# <codecell>


