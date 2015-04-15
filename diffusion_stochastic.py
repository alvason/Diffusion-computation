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
AlvaFigSize = (12, 4);
numberingFig = 0;

# <codecell>

# Gaussian randomness --- Gaussian distribution
minT = float(0)
maxT = float(1000)
totalGPoint_T = int(maxT + 1)
spacingT = np.linspace(minT, maxT, num = totalGPoint_T, retstep = True)
gridT = spacingT[0]
dt = spacingT[1]
GaussSeed = np.random.randn(totalGPoint_T)

#GaussSeed = np.arange(totalGPoint_T)*10.0
#print GaussSeed
# mean = 0
sumG = 0
for i in range(totalGPoint_T):
    sumG = sumG + GaussSeed[i]
meanG = sumG/(totalGPoint_T)

# min-max sorting
def AlvaMinMax(data):
    totalData = np.size(data)
    minMaxListing = np.zeros(totalData)   
    for i in range(totalData):
        # searching the minimum in current array
        jj = 0 
        minMaxListing[i] = data[jj]
        for j in range(totalData - i):
            if data[j] < minMaxListing[i]: 
                minMaxListing[i] = data[j]
                jj = j
        # removing the minmum from current array
        data = np.delete(data, jj)
    return (minMaxListing)


# leveling by using min-max technique
def AlvaLevel(data, totalLevel, normalization = True):
    totalDataPoint = np.size(data)
    minMaxListing = AlvaMinMax(data)
    # searching minimum and maximum values
    minValue = minMaxListing[0]
    maxValue = minMaxListing[-1]
    spacingValue = np.linspace(minValue, maxValue, num = totalLevel + 1, retstep = True)        
    gridLevel = np.delete(spacingValue[0], 0)
    # catogerizing the level set
    # initialize the levelspace by a 'null' space
    levelSpace = np.zeros([2])
    numberLevel = np.zeros([totalLevel])
    jj = 0 # counting the checked number
    for i in range(totalLevel): 
        n = 0 # counting the number in each level
        for j in range(jj, totalDataPoint):
            if minMaxListing[j] <= gridLevel[i]: 
                levelSpace = np.vstack((levelSpace, [i, minMaxListing[j]]))
                n = n + 1
        numberLevel[i] = n
        jj = jj + n
    # delete the inital 'null' space
    levelSpace = np.delete(levelSpace, 0, 0) 
    if normalization == True:
        numberLevel = numberLevel/AlvaMinMax(numberLevel)[-1]
    return (gridLevel, numberLevel, levelSpace)

totalLevel = int(totalGPoint_T/10)
category = AlvaLevel(GaussSeed,totalLevel)
gridLevel = category[0]
numberLevel = category[1]
print category[2].shape
#print numberLevel

numberingFig = numberingFig + 1
plt.figure(numberingFig, figsize = AlvaFigSize)
plt.plot(gridT, GaussSeed, label = 'data')
plt.plot(gridT, AlvaMinMax(GaussSeed), color = 'red', label = 'minMaxListing')
plt.grid(True)
plt.title(r'$ Random \ motion \ (dt = %f,\ mean = %f) $'%(dt, meanG)
          , fontsize = AlvaFontSize)
plt.xlabel(r'$t \ (time)$', fontsize = AlvaFontSize)
plt.ylabel(r'$ Randomness(t) $', fontsize = AlvaFontSize)
plt.legend(loc = (1,0))
plt.show()


numberingFig = numberingFig + 1
plt.figure(numberingFig, figsize = AlvaFigSize)
plt.plot(gridLevel, numberLevel, color = 'red', marker = 'o', label = 'category')
plt.plot(gridLevel, np.exp(-gridLevel**2), label = 'Gaussian')
plt.grid(True)
plt.title(r'$ Gaussian \ distribution\ (data = %i,\ level = %i) $'%(totalGPoint_T, totalLevel)
          , fontsize = AlvaFontSize)
plt.xlabel(r'$ Output \ level$', fontsize = AlvaFontSize)
plt.ylabel(r'$ Number/level $', fontsize = AlvaFontSize)
plt.legend(loc = (1,0))
plt.show()

# <codecell>

# Avarage Many Brownian ways

minT = float(0)
maxT = float(3)
totalGPoint_T = 100
spacingT = np.linspace(minT, maxT, num = totalGPoint_T, retstep = True)
gridT = spacingT[0]
dt = spacingT[1]
totalWay = 10
GaussSeed = np.sqrt(dt)*np.random.randn(totalWay, totalGPoint_T)

gridB = np.zeros([totalWay, totalGPoint_T])
gridB = np.cumsum(GaussSeed, 1)
gridB[:,0] = 0.0 # initialize B(0) = 0

Evaluate = np.exp(gridT + gridB/2)
EvaluateMean = np.mean(Evaluate, axis = 0)

numberingFig = numberingFig + 1;
plt.figure(numberingFig, figsize = AlvaFigSize); 
plt.plot(gridT, Evaluate[::1].T, linewidth = 1)
plt.plot(gridT, EvaluateMean, linewidth = 3, color = 'red')
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
minT = float(0)
maxT = float(3)
totalGPoint_T = 100
spacingT = np.linspace(minT, maxT, num = totalGPoint_T, retstep = True)
gridT = spacingT[0]
dt = spacingT[1]
totalWay = 10
GaussSeed = np.sqrt(dt)*np.random.randn(totalGPoint_T)

# checking Gauss distribution
ddd = np.zeros(totalGPoint_T)
for tn in range(totalGPoint_T):
    ddd[tn] = GaussSeed[tn]**2;
plt.plot(gridT, ddd)

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

np.random.seed(100)

gamma=2; mu=1; Xzero=1
T=1; N=2**8; dt = float(T)/N
t=np.linspace(0,T,N+1)

dW=np.sqrt(dt)*np.random.randn(1,N)
W=np.cumsum(dW)

Xtrue=Xzero*np.exp((gamma-0.5*mu**2)*t[1:]+mu*W); Xtrue=np.insert(Xtrue,0,Xzero)
ax=plt.subplot(111)
ax.plot(t,Xtrue)

R=4; Dt=R*dt; L=float(N)/R
Xem=np.zeros(L+1); Xem[0] = Xzero

for j in xrange(1,int(L)+1):
    Winc=np.sum(dW[0][range(R*(j-1),R*j)])
    Xem[j] = Xem[j-1] + Dt*gamma*Xem[j-1] + mu*Xem[j-1]*Winc

emerr=np.abs(Xem[-1]-Xtrue[-1])
print "Error at endpoint: ", emerr

ax.plot(np.linspace(0,T,L+1),Xem,'r--*')
ax.legend(("exact","em"),loc=2)
plt.show()

# <codecell>


