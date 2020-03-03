#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 13:58:53 2020

@author: ja.acevedo12
"""


import numpy as np
import matplotlib.pyplot as plt


a = 20
l = 900
x1 = 0
x3 = np.sqrt(0.5*a/l)
x2 = -np.sqrt(0.5*a/l)
xeq = np.array([x1,x2,x3])
def U(x, a=1, l=1):
    return -a*x**2 + l*x**4




x = np.linspace(1.5*x2,1.5*x3,1000)

xp = np.linspace(-1*x2,1.5*x2,200)

plt.plot(x,U(x,a,l))
plt.plot(xp,-a*a*0.25/l+2*a*(xp-x2)**2)
plt.scatter(xeq,U(xeq,a,l))
#%%
plt.figure()
def Pp(x,a=a,l=l):
    return 2*a*x-4*l*x**3

def Xp(p,m=1):
    return p/m

def trajectory(x0,p0,N,dt):
    
    tx1 = [x0]
    tp1 = [p0]
    
    N = int(N)
    tmax = N*dt

    for i in range(1,N):
        tp1.append(Pp(tx1[i-1], a, l)*dt + tp1[i-1])
        tx1.append(Xp(tp1[i-1]*dt + tx1[i-1]))
    
    
    tp = np.array(tp1)
    tx = np.array(tx1)
    # plt.figure()
    plt.plot(tx,tp)

    return tx,tp, tmax
tx,tp,tt = trajectory(0,0.1,1e5,1e-4)
tx,tp,tt = trajectory(np.sqrt(a*0.5000/l),0,1e5,1e-4)
tx,tp,tt = trajectory(-np.sqrt(a*0.5000/l),0,1e5,1e-4)
tx,tp,tt = trajectory(1.01*np.sqrt(a*0.5000/l),0,1e5,1e-4)

