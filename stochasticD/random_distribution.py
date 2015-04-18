# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Diffusion computation
# https://github.com/alvason/diffusion-computation
# 
# ### Section003 --- Stochastic solution for the diffusion equation

# <codecell>

'''
author: Alvason Zhenhua Li
date:   03/19/2015
'''

%matplotlib inline

import numpy as np
import matplotlib.pyplot as plt

import alva_machinery as alva

AlvaFontSize = 23;
AlvaFigSize = (12, 4);
numberingFig = 0;

# <codecell>

'''Gaussian randomness --- Gaussian distribution --- Standard normal distribution'''

minT = float(0)
maxT = float(1000)
totalGPoint_T = int(maxT + 1)
spacingT = np.linspace(minT, maxT, num = totalGPoint_T, retstep = True)
gridT = spacingT[0]
dt = spacingT[1]
randomSeed = np.random.standard_normal(totalGPoint_T)
#randomSeed = np.arange(totalGPoint_T)*10.0
#print randomSeed
# mean = 0
sumG = 0
for i in range(totalGPoint_T):
    sumG = sumG + randomSeed[i]
meanG = sumG/(totalGPoint_T)

# leveling by using min-max way
def AlvaLevel(data, totalLevel, normalization = True):
    totalDataPoint = np.size(data)
    minMaxListing = alva.AlvaMinMax(data)
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
        numberLevel = numberLevel/alva.AlvaMinMax(numberLevel)[-1]
    return (gridLevel, numberLevel, levelSpace)

totalLevel = int(totalGPoint_T/10)
category = AlvaLevel(randomSeed, totalLevel)
gridLevel = category[0]
numberLevel = category[1]
print category[2].shape

numberingFig = numberingFig + 1
plt.figure(numberingFig, figsize = AlvaFigSize)
plt.plot(gridT, randomSeed, label = 'data')
plt.plot(gridT, alva.AlvaMinMax(randomSeed), color = 'red', label = 'minMaxListing')
plt.grid(True)
plt.title(r'$ Random \ output \ (dt = %f,\ mean = %f) $'%(dt, meanG)
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
plt.xlabel(r'$ Output-level$', fontsize = AlvaFontSize)
plt.ylabel(r'$ Number/level $', fontsize = AlvaFontSize)
plt.legend(loc = (1,0))
plt.show()

# <codecell>

'''Poisson randomness --- Poisson distribution'''

minT = float(0)
maxT = float(1000)
totalGPoint_T = int(maxT + 1)
spacingT = np.linspace(minT, maxT, num = totalGPoint_T, retstep = True)
gridT = spacingT[0]
dt = spacingT[1]
meanP = 9
randomSeed = np.random.poisson(meanP, totalGPoint_T)

#randomSeed = np.arange(totalGPoint_T)*10.0
#print randomSeed
# calculating the mean
sumP = 0
for i in range(totalGPoint_T):
    sumP = sumP + randomSeed[i]
current_mean = sumP/(totalGPoint_T)

# leveling by using min-max way
def AlvaLevel(data, totalLevel, normalization = True):
    totalDataPoint = np.size(data)
    minMaxListing = alva.AlvaMinMax(data)
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
        numberLevel = numberLevel/alva.AlvaMinMax(numberLevel)[-1]
    return (gridLevel, numberLevel, levelSpace)

totalLevel = int(totalGPoint_T/10)
category = AlvaLevel(randomSeed,totalLevel)
gridLevel = category[0]
numberLevel = category[1]
print category[2].shape

def AlvaPoisson(m, i):
    if type(i) == np.ndarray:
        totalInput = np.size(i)
        poissonD = np.zeros(totalInput)
        for xn in range(totalInput):
            productN = 1
            for j in range(1, int(i[xn]) + 1):
                productN = productN*j
            poissonD[xn] = (m**i[xn])*np.exp(-m)/productN
    else:
        productN = 1
        for j in range(1, int(i) + 1):
            productN = productN*j
        poissonD = (m**i)*np.exp(-m)/productN
    return (poissonD)


numberingFig = numberingFig + 1
plt.figure(numberingFig, figsize = AlvaFigSize)
plt.plot(gridT, randomSeed, label = 'data')
plt.plot(gridT, alva.AlvaMinMax(randomSeed), color = 'red', label = 'minMaxListing')
plt.grid(True)
plt.title(r'$ Random \ output \ (dt = %f,\ mean = %f) $'%(dt, current_mean)
          , fontsize = AlvaFontSize)
plt.xlabel(r'$t \ (time)$', fontsize = AlvaFontSize)
plt.ylabel(r'$ Randomness(t) $', fontsize = AlvaFontSize)
plt.legend(loc = (1,0))
plt.show()

numberingFig = numberingFig + 1
plt.figure(numberingFig, figsize = AlvaFigSize)
plt.plot(gridLevel, numberLevel, color = 'red', marker = 'o', label = 'category')
plt.plot(gridLevel, AlvaPoisson(meanP, gridLevel), label = 'Poisson')
plt.grid(True)
plt.title(r'$ Poisson \ distribution\ (data = %i,\ level = %i) $'%(totalGPoint_T, totalLevel)
          , fontsize = AlvaFontSize)
plt.xlabel(r'$ Output-level$', fontsize = AlvaFontSize)
plt.ylabel(r'$ Number/level $', fontsize = AlvaFontSize)
plt.legend(loc = (1,0))
plt.show()

# <codecell>

def AlvaProduct(i):
    i = np.copy(int(i))
    productN = 1
    for j in range(1, i + 1):        
        productN = productN*j
#       print productN
    return (productN)

def AlvaPoisson(m, i):
    if type(i) == np.ndarray:
        totalInput = np.size(i)
        poissonD = np.zeros(totalInput)  
        for xn in range(totalInput):
            poissonD[xn] = (m**i[xn])*np.exp(-m)/np.prod(np.arange(1, i[xn] + 1))
    return (poissonD)

a = 20
b = 1.0/2
print a*b/2
aaa =np.arange(1, a+1)*b
print aaa
plt.figure(figsize = (12, 4))
plt.plot(aaa, AlvaPoisson(a*b/2, aaa), marker ='o', color = 'red')
plt.grid(True)
plt.show()

# <codecell>

i = 20
print ('Alva = ', AlvaProduct(i))
print ('NumP = ', np.prod(np.arange(1, i + 1)))

