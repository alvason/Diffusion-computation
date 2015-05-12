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

AlvaFontSize = 23
AlvaFigSize = (16, 7)
numberingFig = 0

# leveling by using min-max way
def AlvaLevel(data, totalLevel, normalization = True):
    totalDataPoint = np.size(data)
    minMaxListing = alva.AlvaMinMax(data)
    # searching minimum and maximum values
    minValue = minMaxListing[0]
    maxValue = minMaxListing[-1]
    spacingValue = np.linspace(minValue, maxValue, num = totalLevel + 1, retstep = True)        
    gLevel = np.delete(spacingValue[0], 0)
    # catogerizing the level set
    # initialize the levelspace by a 'null' space
    levelSpace = np.zeros([2])
    numberLevel = np.zeros([totalLevel])
    jj = 0 # counting the checked number
    for i in range(totalLevel): 
        n = 0 # counting the number in each level
        for j in range(jj, totalDataPoint):
            if minMaxListing[j] <= gLevel[i]: 
                levelSpace = np.vstack((levelSpace, [i, minMaxListing[j]]))
                n = n + 1
        numberLevel[i] = n
        jj = jj + n
    # delete the inital 'null' space
    levelSpace = np.delete(levelSpace, 0, 0) 
    if normalization == True:
        numberLevel = numberLevel/alva.AlvaMinMax(numberLevel)[-1]
    return (gLevel, numberLevel, levelSpace)

# <codecell>

'''Gaussian randomness --- Gaussian distribution --- Standard normal distribution'''

totalPoint_Input = int(1000)
gInput = np.arange(totalPoint_Input)
randomSeed = np.random.standard_normal(totalPoint_Input)

sumG = 0
for i in range(totalPoint_Input):
    sumG = sumG + randomSeed[i]
meanG = sumG/(totalPoint_Input)

totalLevel = int(totalPoint_Input/1)
category = AlvaLevel(randomSeed, totalLevel, False)
gLevel = category[0]
numberLevel = category[1]
print category[2].shape

numberingFig = numberingFig + 1
figure = plt.figure(numberingFig, figsize = AlvaFigSize)
plot1 = figure.add_subplot(1, 2, 1)
plot1.plot(gInput, randomSeed, color = 'gray', marker = 'o', label = 'data')
plot1.plot(gInput, alva.AlvaMinMax(randomSeed), color = 'red', marker = 'o', label = 'minMaxListing')
if totalPoint_Input < 100:
    plot1.set_xticks(gInput, minor = True) 
    plot1.set_yticks(randomSeed, minor = True)
    plot1.grid(True, which = 'minor')
else:
    plot1.grid(True, which = 'major')
plt.title(r'$ Random \ output \ (total-input = %i,\ mean = %f) $'%(totalPoint_Input, meanG)
          , fontsize = AlvaFontSize)
plt.xlabel(r'$ input-time $', fontsize = AlvaFontSize)
plt.ylabel(r'$ output $', fontsize = AlvaFontSize)
plt.legend(loc = (0, -0.2))

plot2 = figure.add_subplot(1, 2, 2)
plot2.plot(numberLevel, gLevel, color = 'red', marker = 'o', label = 'category')
plot2.plot(np.exp(-gLevel**2)*alva.AlvaMinMax(numberLevel)[-1], gLevel, color = 'blue', marker = 'o', label = 'Gaussian') 
if totalPoint_Input < 100:
    plot2.set_xticks(numberLevel, minor = True) 
    plot2.set_yticks(gLevel, minor = True)
    plot2.grid(True, which = 'minor')
else:
    plot2.grid(True, which = 'major')
plt.title(r'$ Gaussian \ distribution\ (data = %i,\ level = %i) $'%(totalPoint_Input, totalLevel)
          , fontsize = AlvaFontSize)
plt.xlabel(r'$ Number/level $', fontsize = AlvaFontSize)
plt.ylabel(r'$ Output-level $', fontsize = AlvaFontSize)
plt.legend(loc = (0, -0.2))

figure.tight_layout()
plt.show()

# <codecell>

'''Poisson randomness --- Poisson distribution'''
totalPoint_Input = int(100)
gInput = np.arange(totalPoint_Input)
meanP = 10
randomSeed = np.random.poisson(meanP, totalPoint_Input)

totalLevel = int(totalPoint_Input/1)
category = AlvaLevel(randomSeed, totalLevel, False)
gLevel = category[0]
numberLevel = category[1]
print category[2].shape

# calculating the mean
sumP = 0
for i in range(totalPoint_Input):
    sumP = sumP + randomSeed[i]
current_mean = sumP/(totalPoint_Input)
print ('current \ mean', current_mean)

totalLevel = int(totalPoint_Input/1)
category = AlvaLevel(randomSeed, totalLevel, False)
gLevel = category[0]
numberLevel = category[1]
print category[2].shape

numberingFig = numberingFig + 1
figure = plt.figure(numberingFig, figsize = AlvaFigSize)
plot1 = figure.add_subplot(1, 2, 1)
plot1.plot(gInput, randomSeed, color = 'gray', marker = 'o', label = 'data')
plot1.plot(gInput, alva.AlvaMinMax(randomSeed), color = 'red', marker = 'o', label = 'minMaxListing')
if totalPoint_Input < 100:
    plot1.set_xticks(gInput, minor = True) 
    plot1.set_yticks(randomSeed, minor = True)
    plot1.grid(True, which = 'minor')
else:
    plot1.grid(True, which = 'major')
plt.title(r'$ Random \ output \ (total-input = %i,\ mean = %f) $'%(totalPoint_Input, meanG)
          , fontsize = AlvaFontSize)
plt.xlabel(r'$ input-time $', fontsize = AlvaFontSize)
plt.ylabel(r'$ output $', fontsize = AlvaFontSize)
plt.legend(loc = (0, -0.2))

plot2 = figure.add_subplot(1, 2, 2)
plot2.plot(numberLevel, gLevel, color = 'red', marker = 'o', label = 'category')
plot2.plot(AlvaPoissonD(meanP, gLevel), gLevel, color = 'blue', marker = 'o', label = 'Gaussian') 
if totalPoint_Input < 100:
    plot2.set_xticks(numberLevel, minor = True) 
    plot2.set_yticks(gLevel, minor = True)
    plot2.grid(True, which = 'minor')
else:
    plot2.grid(True, which = 'major')
plt.title(r'$ Gaussian \ distribution\ (data = %i,\ level = %i) $'%(totalPoint_Input, totalLevel)
          , fontsize = AlvaFontSize)
plt.xlabel(r'$ Number/level $', fontsize = AlvaFontSize)
plt.ylabel(r'$ Output-level $', fontsize = AlvaFontSize)
plt.legend(loc = (0, -0.2))

figure.tight_layout()
plt.show()

# <codecell>

i = 30
print ('Alva = ', AlvaProduct(i))
print ('NumP = ', np.prod(np.arange(1, i + 1), dtype=np.int64))

# <codecell>

''' Gaussian Distribution '''
unitD = 1
rangeD = float(100)*unitD
meanD = rangeD/2
varianceD = rangeD
print ('mean = ', mean)
aaa =np.arange(1, rangeD)*b
print aaa
plt.figure(figsize = (12, 4))
#plt.plot(aaa, AlvaPoissonProcess(aaa, N), marker ='^', color = 'blue', label = 'Process')
#plt.plot(aaa, AlvaPoissonD(mean, aaa), marker ='o', color = 'red', label = 'Distribution')
plt.plot(aaa, np.exp(-((aaa - meanD)/varianceD)**2), marker ='+', color = 'red', label = 'Gaussian')
plt.plot(aaa, np.exp(-((aaa - meanD)**2/varianceD)), marker ='+', color = 'blue', label = 'Gaussian')
plt.xlabel(r'$ Output-level$', fontsize = AlvaFontSize)
plt.ylabel(r'$ Number/level $', fontsize = AlvaFontSize)
plt.title(r'$ Poisson \ process $', fontsize = AlvaFontSize)
plt.grid(True)
plt.legend(loc = (1, 0))
plt.show()

# <codecell>

def AlvaIntegrateArea(out_i, min_i, max_i, totalGPoint_i):
    spacing_i = np.linspace(min_i, max_i, num = totalGPoint_i, retstep = True)
    grid_i = spacing_i[0]
    dx = spacing_i[1]
    outArea = 0
    for xn in range(totalGPoint_i):
        outArea = outArea + out_i(grid_i[xn])*dx
    return (outArea)

def gaussianA(i):
    inOut = np.exp(-i**2)
    return (inOut)

ggg = AlvaIntegrateArea(gaussianA, -10, 10, 100)
print ggg
ppp = (np.pi)**(1.0/2)
print ppp

# <codecell>

gInput

# <codecell>

'''Poisson process --- Poisson distribution'''
numberingFig = numberingFig + 1
plt.figure(numberingFig, figsize=(12, 3))
plt.axis('off')
plt.title(r'$ Poisson-distribution \ equation $',fontsize = AlvaFontSize)
plt.text(0,2.0/3,r'$ P_{p}(n|N) = \frac{N!}{n!(N - n)!} p^n (1 - p)^{N - n} $', fontsize = 1.2*AlvaFontSize)
plt.text(0,1.0/3,r'$ P_{m}(n) = \frac{e^{-m} m^n}{n!}, where \ the \ mean \ (average) \ m \equiv pN  $', fontsize = 1.2*AlvaFontSize)
plt.show()

def AlvaProduct(i):
    product = 1
    for j in range(1, int(i) + 1):        
        product = product*j
#       print product
    return (product)

def AlvaPoissonProcess(i, N):
    if type(i) == np.ndarray:
        total_Input = np.size(i)
    p = np.zeros(total_Input)
    process = np.zeros(total_Input)
    for xn in range(total_Input):
        p[xn] = 0.5
        process[xn] = AlvaProduct(N)/(AlvaProduct(i[xn])*AlvaProduct(N - i[xn])
                                      ) * p[xn]**i[xn] * (1 - p[xn])**(N - i[xn])
    print ('probability = ', p)
    return (process)

def AlvaPoissonD(m, i):
    if type(i) == np.ndarray:
        total_Input = np.size(i)
        distribution = np.zeros(total_Input)  
        for xn in range(total_Input):
            distribution[xn] = (m**i[xn])*np.exp(-m)/AlvaProduct(i[xn])
    return (distribution)

N = float(100)
b = 1
mean = N*b/2
rangeN = N*b
print ('mean = ', mean)
aaa =np.arange(1, N)*b
print aaa
plt.figure(figsize = (12, 4))
#plt.plot(aaa, AlvaPoissonProcess(aaa, N), marker ='^', color = 'blue', label = 'Process')
#plt.plot(aaa, AlvaPoissonD(mean, aaa), marker ='o', color = 'red', label = 'Distribution')
plt.plot(aaa, np.exp(-((aaa - mean)/rangeN)**2), marker ='+', color = 'red', label = 'Gaussian')
plt.xlabel(r'$ Output-level$', fontsize = AlvaFontSize)
plt.ylabel(r'$ Number/level $', fontsize = AlvaFontSize)
plt.title(r'$ Poisson \ process $', fontsize = AlvaFontSize)
plt.grid(True)
plt.legend(loc = (1, 0))
plt.show()

# <codecell>


