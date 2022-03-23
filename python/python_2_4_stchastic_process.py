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
from numpy import linspace
import numpy.random as rp

import matplotlib.pyplot as plt

# numer of time samplng & number of stochastic process trial
N_sample = 100
N_realize = 500

# time 
dt = 0.1 # [seconds]
time_init = 0
time_final = dt*N_sample
time = linspace(time_init,time_final,N_sample)

# declare memory space for x_rand_all to include all trials
x_rand_all = np.zeros((N_realize,N_sample))

# time varying mean and sqrt(variance) at the time instance
mu_all = linspace(-2,2,N_sample)
sigma_all = linspace(0.1,1.5,N_sample)

# for a fixed time instance, generate the random numbers
# with the mean and the variance at the fixed time
for idx, (mu_t, sigma_t) in enumerate(zip(mu_all,sigma_all)):
    x_rand = mu_t+sigma_t*rp.randn(N_realize)
    x_rand_all[:,idx] = x_rand

# plot all trials with respect to the time

# Warning: this part is only executed with the small N_trial, 
# e.g., 5
# the plot takes really long and causing the computer crashed 
# with the large N_trial, e.g., 1000
if N_realize < 10:

    fig, ax = plt.subplots(nrows=1,ncols=1)
    ax.plot(time,x_rand_all.transpose(),'k-')
    ax.set_xlabel('time [s]',fontsize=14)
    ax.set_ylabel(r'$x(t)$',fontsize=14)
    ax.set(xlim=(0, time_final),ylim=(-4,6))
    
# approximate mean and variance from the realisation
# and compare with the true
mu_approx = np.mean(x_rand_all,axis=0);
sigma2_approx = np.var(x_rand_all,axis=0)
fig_ms, (ax_ms_0, ax_ms_1) = plt.subplots(nrows=2,ncols=1)
ax_ms_0.plot(time,mu_all)
ax_ms_0.plot(time,mu_approx,'r--')
ax_ms_0.set_ylabel(r'$\mu(t)$',fontsize=14)
ax_ms_0.legend(('True','Estimated'),loc='upper left', fontsize=14)

ax_ms_1.plot(time,sigma_all**2)
ax_ms_1.plot(time,sigma2_approx,'r--')
ax_ms_1.set_ylabel(r'$[\sigma(t)]^2$',fontsize=14)
ax_ms_1.set_xlabel('time [s]',fontsize=14)
ax_ms_1.legend(('True','Estimated'),loc='upper left', fontsize=14);

# esimate the pdf for each instance using N-trials at each instance
N_bin = 100
x_bin = np.linspace(-5,5,N_bin)
dx=np.mean(np.diff(x_bin))
px_all = np.zeros((N_bin-1,N_sample))
for jdx in range(N_sample):
    x_rand = x_rand_all[:,jdx]
    N_occur = np.histogram(x_rand,bins=x_bin)
    N_occur = N_occur[0]
    px_at_t = N_occur/(dx*N_realize)
    px_all[:,jdx] = px_at_t

# plot the estimated pdf
time_matrix, x_bin_matrix = np.meshgrid(time,x_bin[0:-1])

fig_3d = plt.figure()
ax_3d = plt.axes(projection='3d')
ax_3d.plot_surface(time_matrix, x_bin_matrix, px_all, rstride=1, cstride=1,
                cmap='viridis') # viridis, plasma, inferno, magma
ax_3d.set_xlabel('time [s]',fontsize=14)
ax_3d.set_ylabel(r'sampling space $x$',fontsize=14)
ax_3d.set_zlabel(r'$\hat{p}(x)$',fontsize=14)

# a=np.arange(3)
# b=np.arange(5)
# C_mat=np.reshape(np.arange(15),(5,3))
# A_mat,B_mat=np.meshgrid(a,b)
# fig=plt.figure()
# ax=plt.axes(projection='3d')
# ax.plot_surface(A_mat,B_mat,C_mat)
