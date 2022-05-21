# -*- coding: utf-8 -*-
"""
Created on Mon May 16 02:08:40 2022

@author: User
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

rc('font', family='serif')
rc('lines', linewidth=1.5)
rc('font', size=16)
plt.rc('legend',**{'fontsize':12})

M = 34 #space steps
xmin = -2
xmax = 3
a = 1
lmbda = 0.05
theta = a*lmbda
h = 0.15  # time between each discrete time point
k = lmbda*h #space between each discrete x point
N_ = (xmax-xmin)/(k) #timesteps
N = int(N_)
V_nm = np.zeros(shape = (N,M))  #scheme
T =  np.arange(0,1,0.8/N)
x = np.arange(xmin,xmax,h)

def v_NM(t,x):
    if(np.abs(x-t)>1):
        return 0
    return 1-np.abs(x-t)

#populate initial values at t = 0 
for i in range(0, M):
    V_nm[0][i] = v_NM(0, x[i])
  
#implement periodic boundary conditions    
for j in range(0, N):
    V_nm[j][0] = 0
    #V_nm[j][M-1] = 0

#Populate the time-space grid 
for i in range(1, N):
    for j in range(1, M):
        V_nm[i][j] = v_NM(T[i], x[j])
        

U = V_nm  # exact solution
Um1 = np.roll(V_nm,1) #scheme at m+1
Um_1 = np.roll(V_nm,-1) #scheme at m-1

# Compute the scheme
U_n1  = 0.5*(Um1-Um_1) + 0.5*(theta)*(-Um_1+Um1) 

#More Boundary Conditions
for j in range(0, M):
    U_n1[j][M-1] = V_nm[j][M-1]

#Plotting of Scheme
fig, ax = plt.subplots(figsize=(5.5,4))
for i in range(0, N):
  plt.clf()
  plt.plot(x,V_nm[i],color = 'red')
  plt.scatter(x,U[i], marker='o', facecolors='white', color='k')
  plt.gca().legend(('Scheme ','Exact solution'))
  plt.axis([xmin, xmax, 0, 1.4])
  plt.title('t='+str(round(T[i],3)),fontsize=16)
  plt.xlabel('x',fontsize=18)
  plt.ylabel('u',fontsize=18)
  plt.subplots_adjust(left=0.2)
  plt.subplots_adjust(bottom=0.18)
  plt.draw()
  plt.pause(0.001)