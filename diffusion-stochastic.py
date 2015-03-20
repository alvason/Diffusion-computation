# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Diffusion computation
# https://github.com/alvason/diffusion-computation
# 
# ### Lecture003 --- Stochastic solution for the diffusion equation

# <codecell>

'''
author: Alvason Zhenhua Li
date:   03/19/2015
'''

%matplotlib inline

import numpy as np
import matplotlib.pyplot as plt
import time
import IPython.display as idisplay
from mpl_toolkits.mplot3d.axes3d import Axes3D

AlvaFontSize = 23;
AlvaFigSize = (9, 6);
numberingFig = 0;

# <codecell>

# Avarage Many Brownian ways

minT = float(0); maxT = float(3);
totalGPoint_T = 100; 
dt = (maxT - minT)/totalGPoint_T;
gridT = np.linspace(minT, maxT, totalGPoint_T);

totalWay = 10;
GaussSeed = np.sqrt(dt)*np.random.randn(totalWay, totalGPoint_T); 

gridB = np.zeros([totalWay, totalGPoint_T]);
gridB = np.cumsum(GaussSeed, 1); 
gridB[:,0] = 0.0; # initialize B(0) = 0

Evaluate = np.exp(gridT + gridB/2)
EvaluateMean = np.mean(Evaluate, axis = 0)

numberingFig = numberingFig + 1;
plt.figure(numberingFig, figsize = AlvaFigSize); 
plt.plot(gridT, Evaluate[::1].T, linewidth = 1)
plt.plot(gridT, EvaluateMean, linewidth = 3, color = 'r')
plt.grid(True)
plt.title(r'$ Brownian \ motion: (dt = %f,\ motion = %i) $'%(dt, totalGPoint_T)
          , fontsize = AlvaFontSize);
plt.xlabel(r'$t \ (time)$', fontsize = AlvaFontSize); 
plt.ylabel(r'$ Evaluate(t) $', fontsize = AlvaFontSize);
plt.text(maxT, 1, r'$ Evaluate(B(t)) = e^{t + \frac{1}{2}B(t)} $'
         , fontsize = AlvaFontSize);
# plt.legend(('Eva ways', 'Eva Mean'), loc = (1, 1))
plt.show()

# <codecell>

# Many Brownian ways

minT = float(0); maxT = float(3);
totalGPoint_T = 100; 
dt = (maxT - minT)/totalGPoint_T;
gridT = np.linspace(minT, maxT, totalGPoint_T);

totalWay = 30;
GaussSeed = np.sqrt(dt)*np.random.randn(totalWay, totalGPoint_T); 

gridB = np.zeros([totalWay, totalGPoint_T]);
gridB = np.cumsum(GaussSeed, 1); 
gridB[:,0] = 0.0; # initialize B(0) = 0

numberingFig = numberingFig + 1;
plt.figure(numberingFig, figsize = AlvaFigSize);     
plt.plot(gridT, gridB.T);
plt.grid(True)
plt.title(r'$ Brownian \ motion: (dt = %f,\ motion = %i) $'%(dt, totalGPoint_T)
          , fontsize = AlvaFontSize);
plt.xlabel(r'$t \ (time)$', fontsize = AlvaFontSize); 
plt.ylabel(r'$ B(t) $', fontsize = AlvaFontSize);
plt.text(maxT, minT, r'$ B(t_{i+1}) = B(t_i) + (t_{i+1} - t_i)^{1/2}Gauss_{i+1} $'
         , fontsize = AlvaFontSize);
plt.show()

# <codecell>

# Brownian motion

minT = float(0); maxT = float(3);
totalGPoint_T = 100; 
dt = (maxT - minT)/totalGPoint_T;

GaussSeed = np.sqrt(dt)*np.random.randn(totalGPoint_T); 

gridT = np.linspace(minT, maxT, totalGPoint_T); 

gridB = np.zeros(totalGPoint_T);

for tn in range(totalGPoint_T - 1):
    gridB[tn + 1] = gridB[tn] + GaussSeed[tn + 1]
    
# gridB = np.cumsum(GaussSeed); 
# gridB[0] = 0.0; # initialize B(0) = 0

numberingFig = numberingFig + 1;
plt.figure(numberingFig, figsize = AlvaFigSize);     
plt.plot(gridT, gridB, marker = 'x', label = 'Brownian');
plt.plot(gridT, GaussSeed, marker = '+', label = 'Gaussian');
plt.grid(True)
plt.title(r'$ Brownian \ motion: (dt = %f,\ motion = %i) $'%(dt, totalGPoint_T)
          , fontsize = AlvaFontSize);
plt.xlabel(r'$t \ (time)$', fontsize = AlvaFontSize); 
plt.ylabel(r'$ B(t) $', fontsize = AlvaFontSize);
plt.text(maxT, minT, r'$ B(t_{i+1}) = B(t_i) + (t_{i+1} - t_i)^{1/2}Gauss_{i+1} $'
         , fontsize = AlvaFontSize);
plt.legend()
plt.show()

# <codecell>

aaa = np.zeros([2,4])
aaa.T

# <codecell>


