#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2022 Jongrae.K

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import numpy as np
import matplotlib.pyplot as plt

# Set initial values & change non-SI units into the SI Units

# number of stochastic process trial
N_realize = 10

dt = 0.1 # [seconds]
time_init = 0
time_final = 120 # [seconds]
N_sample = int(time_final/dt) + 1
time = np.linspace(time_init,time_final, N_sample)

# standard deviation of the bias, sigma_beta_xyz
sigma_beta_xyz = np.array([0.01, 0.01, 0.02]) # [degrees/sqrt(s)]
sigma_beta_xyz = sigma_beta_xyz*(np.pi/180) # [rad/sqrt(s)]
sigma_eta_xyz = sigma_beta_xyz/np.sqrt(dt)

# store all beta history
beta_all = np.zeros((N_sample,3))

# multiple realization of bias noise
fig, (ax0,ax1,ax2) = plt.subplots(nrows=3,ncols=1)
for idx in range(N_realize):
    # initial beta(t)
    beta = (2*np.random.rand(3)-1)*0.03 # +/- 0.03[degrees/s]
    
    # from here all units are in SI except the plot parts
    beta = beta*(np.pi/180) # [radians/s]
    
    # main simulation loops
    for jdx in range(N_sample):
        beta_all[jdx,:] = beta
        
        eta_u = sigma_eta_xyz*np.random.randn(3)
        dbeta = eta_u*dt
        beta = beta + dbeta
        
        beta_all[jdx,:]= beta
    
    # plot all realization of beta in degrees/s
    ax0.plot(time,beta_all[:,0]*180/np.pi,'r-')
    ax0.set(xlim=(0, time_final),ylim=(-0.5,0.5))
    ax0.set_ylabel(r'$\beta_x(t) [^\circ/s]$',fontsize=18)
    ax0.set_yticklabels(np.linspace(-0.5,0.5,5),fontsize=18)
    ax0.set_xticklabels([])
    ax0.set_title('Python',fontsize=18)
    
    ax1.plot(time,beta_all[:,1]*180/np.pi,'b-')
    ax1.set_ylabel(r'$\beta_y(t) [^\circ/s]$',fontsize=18)
    ax1.set(xlim=(0, time_final),ylim=(-0.5,0.5))
    ax1.set_yticklabels(np.linspace(-0.5,0.5,5),fontsize=18)
    ax1.set_xticklabels([])
    
    ax2.plot(time,beta_all[:,2]*180/np.pi,'g-')
    ax2.set_ylabel(r'$\beta_z(t) [^\circ/s]$',fontsize=18)
    ax2.set(xlim=(0, time_final),ylim=(-0.5,0.5))
    ax2.set_xlabel('time [s]',fontsize=18);
    ax2.set_yticklabels(np.linspace(-0.5,0.5,5),fontsize=18)
    ax2.set_xticklabels(np.linspace(0,120,7),fontsize=18)

fig.set_size_inches(9,6)    
fig.savefig('bias_noise_simulation_python.pdf',dpi=250)
