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
dt = 0.05 # [seconds]
time_init = 0
time_final = 120 # [seconds]
N_sample = int(time_final/dt) + 1
time = np.linspace(time_init,time_final, N_sample)

# standard deviation of the bias, sigma_beta_xyz
sigma_beta_xyz = np.array([0.05, 0.04, 0.06]) # [degrees/sqrt(s)]
sigma_beta_xyz = sigma_beta_xyz*(np.pi/180) # [rad/sqrt(s)]
sigma_eta_xyz = sigma_beta_xyz/np.sqrt(dt)

# standard devitation of the white noise, sigma_v
sigma_v = 0.01 #[degrees/s]
sigma_v = sigma_v*(np.pi/180) #[rad/s]

# initial beta(t)
beta = (2*np.random.rand(3)-1)*0.03 # +/- 0.03[degrees/s]
beta = beta*(np.pi/180) # [radians/s]

# prepare the data store
w_all = np.zeros((N_sample,3))
w_measure_all = np.zeros((N_sample,3))

# main simulation loops
for idx in range(N_sample):
    
    time_c = time[idx]
    w_true = np.array([ 0.1*np.sin(2*np.pi*0.005*time_c), # [rad/s]
                        0.05*np.cos(2*np.pi*0.01*time_c + 0.2), #[rad/s]
                        0.02 #[rad/s]
                    ])
    
    # beta(t)
    eta_u = sigma_eta_xyz*np.random.randn(3)
    dbeta = eta_u*dt
    beta = beta + dbeta
    
    # eta_v(t)
    eta_v = sigma_v*np.random.randn(3)
    
    # w_tilde
    w_measurement = w_true + beta + eta_v
    
    # store history
    w_all[idx,:] = w_true
    w_measure_all[idx,:] = w_measurement

# plot all realization of beta in degrees/s
fig, ax = plt.subplots(nrows=1,ncols=1)
ax.plot(time,w_all*180/np.pi)
ax.plot(time,w_measure_all*180/np.pi,'--')
ax.set_ylabel(r'$[^\circ/s]$',fontsize=14);
ax.set_xlabel(r'time [s]',fontsize=14);
ax.legend((r'$\omega_x$',r'$\omega_y$',r'$\omega_z$',
    r'$\tilde{\omega}_x$',r'$\tilde{\omega}_y$',r'$\tilde{\omega}_z$'),
          fontsize=14, loc='lower left')
ax.set(xlim=(0, time_final),ylim=(-4,6))

fig.set_size_inches(9,6)    
#fig.savefig('gyro_measurement_python.pdf',dpi=250)
