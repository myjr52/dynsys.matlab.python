#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 08:31:59 2021

@author: jongrae
"""

import numpy as np
import matplotlib.pyplot as plt

N_alpha   = 100
alpha_all = np.linspace(0.05,0.4,N_alpha)
max_mu    = np.zeros(N_alpha)

for gdx, alpha in enumerate(alpha_all):

    A = np.array([[-2]])
    B = np.array([[0.1,0.5,0.1/alpha]])
    C = np.array([[1],[1],[0]])
    D = np.array([[0,0,0],[0,0,0],[0,alpha,0]])
    
    N_omega = 300
    omega = np.logspace(-2,3,N_omega)
    mu_ub = np.zeros(N_omega)
    
    for idx in range(N_omega):
        jw = complex(0,omega[idx])
        Mjw = C@np.linalg.inv(jw-A)@B+D
        U,S,V=np.linalg.svd(Mjw)
    
        mu_ub[idx] = S.max()
    
    max_mu[gdx] = mu_ub.max()

fig1, ax = plt.subplots(nrows=1,ncols=1)
ax.plot(alpha_all,max_mu)
ax.axis([alpha_all[0],alpha_all[-1],0.4,1.6])
ax.set_ylabel(r'$\max \{\bar{\sigma}\, [M(j\omega)]\}$',fontsize=14)
ax.set_xlabel(r'$\alpha$',fontsize=14)
#fig1.savefig('max_mu_ub_alpha.pdf')
