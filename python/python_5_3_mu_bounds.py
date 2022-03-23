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

alpha = 0.266

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
    

fig1, ax = plt.subplots(nrows=1,ncols=1)
ax.semilogx(omega,mu_ub)
ax.axis([1e-2,1e3,0,1.1])
ax.set_ylabel(r'$\bar{\sigma}\, [M(j\omega)]$',fontsize=14)
ax.set_xlabel(r'$\omega$ [rad/time]',fontsize=14)
