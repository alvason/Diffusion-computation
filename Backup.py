# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

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
print ('total computational time = %f'% total_time);
#print gridHtx  
numberingFig = numberingFig + 1;
plt.figure(numberingFig, figsize = AlvaFigSize);     
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
plt.figure(numberingFig, figsize = AlvaFigSize);     
plt.plot(gridX, gridHtx.T);
plt.grid(True)
plt.title(r'$ Numerical \ solution: (dt = %f,\ dx = %f) $'%(dt, dx), fontsize = AlvaFontSize);
plt.xlabel(r'$x \ (space)$', fontsize = AlvaFontSize); plt.ylabel(r'$H(x,t)$', fontsize = AlvaFontSize);
plt.show()

numberingFig = numberingFig + 1;
plt.figure(numberingFig, figsize = AlvaFigSize); 
plt.pcolor(X, Y, gridHtx);
plt.title(r'$ Numerical \ solution: (dt = %f,\ dx = %f) $'%(dt, dx), fontsize = AlvaFontSize);
plt.xlabel(r'$x \ (space)$', fontsize = AlvaFontSize); plt.ylabel(r'$H(x,t)$', fontsize = AlvaFontSize);
#plt.colorbar();
plt.show()

