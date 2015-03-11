# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

'''
author: Alvason Zhenhua Li
date:   03/11/2015
'''

%matplotlib inline

import numpy as np
import matplotlib.pyplot as plt
import time1
import IPython.display as idisplay
from mpl_toolkits.mplot3d.axes3d import Axes3D

AlvaFontSize = 23;

plt.figure(1, figsize=(12,3))
plt.axis('off')
plt.title('Diffusion Equation',fontsize = AlvaFontSize)
plt.text(0,1.0/2,r'$ \frac{\partial H(x,t)}{\partial t}=\xi \ \frac{\partial^2 H(x,t)}{\partial x^2} $', fontsize = AlvaFontSize)
plt.show()

plt.figure(2, figsize=(12,3))
plt.axis('off')
plt.title('Analytic Solution', fontsize = AlvaFontSize)
plt.text(0,1.0/2,r'$H(x,t) = \frac{1}{(1 + 4 \ \xi \ t)^{1/2}} e^\frac{-x^2}{1 + 4 \ \xi \ t}}$', fontsize = AlvaFontSize)
plt.show()

# <codecell>

def AlvaGridXX(gridX, totalGPointY):
    gridXX = gridX;
    for n in range(totalGPointY - 1):
        gridXX = np.vstack((gridXX, gridX));
    return gridXX;

minX = 0; maxX = 1;
minT = 0; maxT = 3;

totalGPointX = 30 + 1;
gridX = np.linspace(minX, maxX, totalGPointX)
totalGPointT = 6 + 1;
gridT = np.linspace(minT, maxT, totalGPointT)

X = AlvaGridXX(gridX, totalGPointT)
Y = AlvaGridXX(gridT, totalGPointX).T

plt.figure(3, figsize=(6,6));
plt.xticks(fontsize = AlvaFontSize);plt.yticks(fontsize = AlvaFontSize)
plt.plot(X, Y, linestyle='dotted',marker='o', color = 'w');
plt.title(r'$GridYX$', fontsize = AlvaFontSize); 
plt.xlabel(r'$X$', fontsize = AlvaFontSize); plt.ylabel(r'$Y$', fontsize = AlvaFontSize)
plt.show()

# <codecell>

minX = 0; maxX = 10;
minT = 0; maxT = 20;

totalGPointX = 100 + 1;
gridX = np.linspace(minX, maxX, totalGPointX)

totalGPointT = 10 + 1;
gridT = np.linspace(minT, maxT, totalGPointT)

gridHtx = np.zeros([totalGPointT, totalGPointX])

movingRate = 1.0/10.0 # diffusion coefficience

# Analytic solution
for pT in range(totalGPointT):  
    for pX in range(totalGPointX):
        gridHtx[pT,pX] = (1.0/(1.0+4.0*movingRate*gridT[pT]))*np.exp(-(gridX[pX]-(maxX-minX)/2.0)**2/(1.0+4.0*movingRate*gridT[pT]))
        
plt.figure(3,figsize=(10,5));     
plt.plot(gridX[:], gridHtx[:,:].T);
plt.grid(True)
plt.title(r'$Analytic \ solution$', fontsize = AlvaFontSize);
plt.xlabel(r'$x \ (space)$', fontsize = AlvaFontSize); plt.ylabel(r'$H(x,t)$', fontsize = AlvaFontSize)
plt.text(maxX,1.0/2,r'$H(t,x) = \frac{1}{(1 + 4 \ \xi \ t)^{1/2}} e^\frac{-x^2}{1 + 4 \ \xi \ t}}$', fontsize = AlvaFontSize)
plt.show()

# <codecell>

X = AlvaGridXX(gridX, totalGPointT); 
Y = AlvaGridXX(gridT, totalGPointX).T; 
Z = gridHtx;

view,figure = plt.subplots(1,2,figsize=(30,10));
figure[0].pcolormesh(X, Y, gridHtx); 
plt.sca(figure[0])
plt.xticks(fontsize = AlvaFontSize);plt.yticks(fontsize = AlvaFontSize)
figure[0].set_title("Analytic diffusion", fontsize = AlvaFontSize); 
figure[0].set_xlabel(r'x (space)', fontsize = AlvaFontSize);
figure[0].set_ylabel(r't (time)', fontsize = AlvaFontSize); 

figure[1].contour(X, Y, gridHtx, vmin=abs(gridHtx).min(), vmax=abs(gridHtx).max());
plt.sca(figure[1])
plt.xticks(fontsize = AlvaFontSize); plt.yticks(fontsize = AlvaFontSize)
figure[1].set_title("Analytic diffusion", fontsize = AlvaFontSize);
figure[1].set_xlabel(r'x (space)', fontsize = AlvaFontSize);
figure[1].set_ylabel(r't (time)', fontsize = AlvaFontSize);

view.tight_layout()

# <codecell>

figure = plt.figure();

figure1 = Axes3D(figure)
figure1.view_init(30, 80)

figure1.plot_surface(X, Y, Z, rstride = 1, cstride = 6, alpha = 0.6, cmap = 'jet');
plt.xlabel(r'$x \ (space)$', fontsize = AlvaFontSize); plt.ylabel(r'$t \ (time)$', fontsize = AlvaFontSize)
plt.show()

# <codecell>

minY = minT; maxY = maxT;
minZ = 0; maxZ = 1;

figure = plt.figure(figsize=(20,10));

figure1 = figure.add_subplot(1,2,1, projection='3d')
figure1.view_init(30, -60)

plt.xticks(fontsize = AlvaFontSize);plt.yticks(fontsize = AlvaFontSize)
figure1.plot_surface(X, Y, Z, rstride=10, cstride=10, alpha=0.2);
figure1.contour(X, Y, Z, zdir='x', offset=-(maxX-minX));
figure1.contour(X, Y, Z, zdir='y', offset=(maxY-minY));
figure1.contour(X, Y, Z, zdir='z', offset=-(maxZ-minZ));
figure1.set_xlabel('X');
figure1.set_xlim(minX-maxX, maxX+maxX);
figure1.set_ylabel('Y');
figure1.set_ylim(minY-maxY, maxY+maxY);
figure1.set_zlabel('Z');
figure1.set_zlim(minZ-maxZ, maxZ+maxZ);

figure1.set_xlabel(r'Space', fontsize = AlvaFontSize);
figure1.set_ylabel(r'Time', fontsize = AlvaFontSize);


figure2 = figure.add_subplot(1,2,2, projection='3d')
figure2.view_init(30, 60)

plt.xticks(fontsize = AlvaFontSize);plt.yticks(fontsize = AlvaFontSize);
figure2.plot_surface(X, Y, Z, rstride=10, cstride=10, alpha=0.2);
figure2.contour(X, Y, Z, zdir='x', offset=-(maxX-minX));
figure2.contour(X, Y, Z, zdir='y', offset=-(maxY-minY));
figure2.contour(X, Y, Z, zdir='z', offset=-(maxZ-minZ));
figure2.set_xlabel('X');
figure2.set_xlim(minX-maxX, maxX+maxX);
figure2.set_ylabel('Y');
figure2.set_ylim(minY-maxY, maxY+maxY);
figure2.set_zlabel('Z');
figure2.set_zlim(minZ-maxZ, maxZ+maxZ);

figure2.set_xlabel(r'Space', fontsize = AlvaFontSize);
figure2.set_ylabel(r'Time', fontsize = AlvaFontSize);


figure.tight_layout()

# <codecell>

gridHtx

# <codecell>

def AlvaD2nd(Htx):
    leftX = H[0]

# <codecell>

def d2nd(F, Fi):
    leftX = Fi[0,0:-2]
    centerX = -2.0*Fi[0,1:-1]
    rightX = Fi[0,2:]
    F[0,1:-1] = Fi[0,1:-1] + c0*dt*((leftX + centerX + rightX)/dx**2)      

def d2ndd(F, Fi):
    left2F = -Fi[0,0:-3]
    leftF = 16*Fi[0,0:-2]
    centerF = -30*Fi[0,1:-1]
    rightF = 16*Fi[0,2:]
    right2F = -Fi[0,3:]
    F[0,1:-1] = Fi[0,1:-1] + dt*((left2F + leftF + centerF + rightF +right2F)/(12*dx**2) + Fi[0,1:-1] - Fi[0,1:-1]**2)

Fxt = np.zeros([nt,nx])
    
tstart = time.time()
for tStep in range(nt):
    Fxt[tStep] = F
    for i in range(int(tStep/dt)):
        d2nd(F,Fi)
        Fi = F
#    print "computing F for time-step =", tStep
tfinish = time.time()

# <codecell>

X,Y = np.meshgrid(nx,nt)
Z = Fxt

plt.figure(4,figsize=(9,6))

plt.pcolor(Fxt);
plt.title("Fisher equation");
plt.xlabel(r'$Space$', fontsize=18)
plt.ylabel(r'$Time$', fontsize=18)
plt.axes().set_aspect('auto');
plt.colorbar();

plt.figure(5,figsize=(9,6))
plt.contour(Z, vmin=abs(Z).min(), vmax=abs(Z).max(), extent=[xmin,xmax,tmin,tmax])
plt.title("Fisher equation");
plt.xticks(np.arange(xmin,xmax,10)); plt.yticks(np.arange(tmin,tmax,1));
plt.xlabel(r'$Space$', fontsize=18)
plt.ylabel(r'$Time$', fontsize=18)
plt.axes().set_aspect('auto');
plt.colorbar();
plt.show()

