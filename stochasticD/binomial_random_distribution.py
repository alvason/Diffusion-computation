# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Diffusion computation
# https://github.com/alvason/diffusion-computation
# 
# ### Section003 --- Stochastic solution for the diffusion equation
# ##### Random distribution --- Binomial distribution (discrete)

# <codecell>

'''
author: Alvason Zhenhua Li
date:   03/19/2015
'''

%matplotlib inline

import numpy as np
import matplotlib.pyplot as plt

import alva_machinery_diffusion as alva

AlvaFontSize = 23
AlvaFigSize = (16, 7)
numberingFig = 0

'''Binomial process --- Binomial distribution'''
numberingFig = numberingFig + 1
plt.figure(numberingFig, figsize=(12, 3))
plt.axis('off')
plt.title(r'$ Binomial-distribution \ equation $',fontsize = AlvaFontSize)
plt.text(0,2.0/3,r'$ P_{b}(n|N) = \frac{N!}{n!(N - n)!} p^n (1 - p)^{N - n} $', fontsize = 1.2*AlvaFontSize)
plt.show()

# <codecell>

def AlvaProduct(i):
    product = 1
    for j in range(1, int(i) + 1):        
        product = product*j
#       print product
    return (product)

def AlvaBinomialD(i, N):
    if type(i) == np.ndarray:
        total_Input = np.size(i)
    p = np.zeros(total_Input)
    distribution = np.zeros(total_Input)
    for xn in range(total_Input):
        p[xn] = 0.5
        distribution[xn] = AlvaProduct(N)/(AlvaProduct(i[xn])*AlvaProduct(N - i[xn])
                                      ) * p[xn]**i[xn] * (1 - p[xn])**(N - i[xn])
#    print ('probability = ', p)
    return (distribution)

N = float(100)
b = 1
mean = N*b/2
rangeN = N*b
print ('mean = ', mean)
aaa =np.arange(1, N)*b
print aaa
plt.figure(figsize = (12, 4))
plt.plot(aaa, AlvaBinomialProcess(aaa, N), marker ='^', color = 'blue', label = 'Distribution')
plt.xlabel(r'$ output-level$', fontsize = AlvaFontSize)
plt.ylabel(r'$ input/level $', fontsize = AlvaFontSize)
plt.title(r'$ Binomial \ distribution $', fontsize = AlvaFontSize)
plt.grid(True)
plt.legend(loc = (1, 0))
plt.show()

# <codecell>

'''Binomial randomness --- Binomial distribution (discrete)'''

totalPoint_Input = int(100)
gInput = np.arange(totalPoint_Input)
output_level = 10
probability_peak = 0.9
randomSeed = np.random.binomial(output_level, probability_peak, totalPoint_Input)

sumO = 0
for i in range(totalPoint_Input):
    sumO = sumO + randomSeed[i]
meanO = sumO/(totalPoint_Input)

totalLevel = int(totalPoint_Input/1)
category = alva.AlvaLevel(randomSeed, totalLevel, False)
gLevel = category[0]
numberLevel = category[1]
print category[2].shape

numberingFig = numberingFig + 1
figure = plt.figure(numberingFig, figsize = AlvaFigSize)
plot1 = figure.add_subplot(1, 2, 1)
plot1.plot(gInput, randomSeed, color = 'gray', marker = 'o', label = 'data')
plot1.plot(gInput, alva.AlvaMinMax(randomSeed), color = 'red', marker = 'o', label = 'minMaxSorting')
if totalPoint_Input < 100:
    plot1.set_xticks(gInput, minor = True) 
    plot1.set_yticks(randomSeed, minor = True)
    plot1.grid(True, which = 'minor')
else:
    plot1.grid(True, which = 'major')
plt.title(r'$ Binomial \ (total-input = %i,\ mean = %f) $'%(totalPoint_Input, meanO)
          , fontsize = AlvaFontSize)
plt.xlabel(r'$ input $', fontsize = AlvaFontSize)
plt.ylabel(r'$ output $', fontsize = AlvaFontSize)
plt.legend(loc = (0, -0.2))

plot2 = figure.add_subplot(1, 2, 2)
plot2.plot(numberLevel, gLevel, color = 'red', marker = 'o', label = 'category') 
if totalPoint_Input < 100:
    plot2.set_xticks(numberLevel, minor = True) 
    plot2.set_yticks(gLevel, minor = True)
    plot2.grid(True, which = 'minor')
else:
    plot2.grid(True, which = 'major')
plt.title(r'$ Binomial \ distribution\ (data = %i,\ level = %i) $'%(totalPoint_Input, totalLevel)
          , fontsize = AlvaFontSize)
plt.xlabel(r'$ input/level $', fontsize = AlvaFontSize)
plt.ylabel(r'$ output-level $', fontsize = AlvaFontSize)
plt.legend(loc = (0, -0.2))

figure.tight_layout()
plt.show()

# <codecell>


