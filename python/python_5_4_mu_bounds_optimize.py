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
