# -*- coding: utf-8 -*-
"""
Created on Thu May 12 21:15:35 2022

@author: User
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

rc('font', family='serif')
rc('lines', linewidth=1.5)
rc('font', size=16)
plt.rc('legend',**{'fontsize':12})

M = int(input("Please enter the space steps (M): ")) #space steps
h = float(input("Please enter the time between each discrete time point (h): "))# time between each discrete time point
xmin = -1
xmax = 2
a = 1.1
lmbda = 0.05
theta = a*lmbda
k = lmbda*h #space between each discrete x point
N_ = (xmax-xmin)/(k) #timesteps
N = int(N_)
V_nm = np.zeros(shape = (N,M))  #scheme
T =  np.arange(0,1,1/N)
x = np.arange(xmin,xmax,h)


def v_NM(t,x):
    return np.e**(-10*(x-a*t- 0.1)**2)
#populate initial values at t = 0 
for i in range(0, M):
    V_nm[0][i] = v_NM(0, x[i])
  
#implement periodic boundary conditions    
for j in range(0, N):
    V_nm[j][0] = 0
    V_nm[j][M-1] = 0

#Populate the time-space grid 
for i in range(1, N):
    for j in range(1, M-1):
        V_nm[i][j] = v_NM(T[i], x[j])
        
U = V_nm  # exact solution

Um1 = np.roll(V_nm,1) #scheme at m+1
Um_1 = np.roll(V_nm,-1) #scheme at m-1

#scheme
U_n1  = V_nm - 0.5*theta*(Um1-Um_1) + 0.5*(theta**2)*(Um_1-2*V_nm+Um1) 

#Error Analysis

#Computation of error at h/2
# Here we compute the scheme again, but this time let h = h/2 then return
# norm(u(h/2)-v(h/2))

def error_of_h_2(h):
    V_nmh =  np.zeros(shape = (N,M))
    xh = np.arange(xmin,xmax,h)
    #populate initial values at t = 0 
    for i in range(0, M):
        V_nmh[0][i] = v_NM(0, xh[i])
  
    #implement periodic boundary conditions    
    for j in range(0, N):
        V_nmh[j][0] = 0
        V_nmh[j][M-1] = 0

    #Populate the time-space grid 
    for i in range(1, N):
        for j in range(1, M-1):
            V_nmh[i][j] = v_NM(T[i], x[j])
        
    Um1h = np.roll(V_nm,1) #scheme at m+1
    Um_1h = np.roll(V_nm,-1) #scheme at m-1

    #scheme
    U_n1h  = V_nmh - 0.5*theta*(Um1h-Um_1h) + 0.5*(theta**2)*(Um_1h-2*V_nmh+Um1h) 
    return np.linalg.norm(U_n1h-V_nmh)

eh = np.linalg.norm(U_n1-V_nm)
eh_2 = error_of_h_2(h/2)
errorRatio = eh/eh_2
order = 2 + errorRatio
print("")
print ('M = ', M, ', h = ',h,', Error = ',eh,', Error Ratio= ',errorRatio,', Observed Order= ', order)

#Plotting of Scheme
fig, ax = plt.subplots(figsize=(5.5,4))
for i in range(0, N):
  plt.clf()
  plt.plot(x,V_nm[i])
  plt.scatter(x,U[i], marker='o', facecolors='white', color='k')
  plt.gca().legend(('Lax-Wendroff scheme ','Exact solution'))
  plt.axis([xmin, xmax, 0, 1.4])
  plt.title('t='+str(round(T[i],3)),fontsize=16)
  plt.xlabel('x',fontsize=18)
  plt.ylabel('u',fontsize=18)
  plt.subplots_adjust(left=0.2)
  plt.subplots_adjust(bottom=0.18)
  plt.draw()
  plt.pause(0.001)

plt.show()

