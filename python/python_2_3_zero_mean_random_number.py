#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 12:23:16 2020

@author: jongrae
"""

import numpy as np
from numpy import linspace
import numpy.random as rp

import matplotlib.pyplot as plt

# true probability density function (pdf)
var_x = 1;
mean_x = 0;
Omega_x = linspace(-5,5,1000);
px = (1/(np.sqrt(2*np.pi*var_x)))*np.exp(-(Omega_x-mean_x)**2/(2*var_x));

fig, ax = plt.subplots(nrows=1,ncols=1)
ax.plot(Omega_x,px,linewidth=3)

# generate N random numbers with the mean zero and the variance 1 using
# numpy.random.randn
N_all = np.array([100,10000])
x_bin = linspace(-5,5,30)
dx=np.mean(np.diff(x_bin))
line_style = ['rs-','go-']
for N_trial, lnsty in zip(N_all,line_style):
    x_rand = rp.randn(1,N_trial)
    
    # number of occurance of x_rand in x_bin
    N_occur = np.histogram(x_rand,bins=x_bin)
    N_occur = N_occur[0]
    
    ax.plot(x_bin[0:-1]+dx/2, N_occur/(dx*N_trial),lnsty);

ax.set_xlabel(r'Random Variable x Sampling Space: $\Omega_x$',fontsize=14)
ax.set_ylabel('probability density function',fontsize=14)
ax.legend((r'True $p(x)$','N=100','N=10,000'),loc='upper right',fontsize=14)

fig.savefig('compare_mu_sgm2_true_estimated_python.pdf')