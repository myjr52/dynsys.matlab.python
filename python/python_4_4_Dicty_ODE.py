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
from scipy.integrate import solve_ivp

def Dicty_cAMP(time,state,ki_para):
    ACA, PKA, ERK2, REGA, icAMP, ecAMP, CAR1  = state

    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14 = ki_para
    
    dACA_dt   = k1*CAR1 - k2*ACA*PKA
    dPKA_dt   = k3*icAMP - k4*PKA
    dERK2_dt  = k5*CAR1 - k6*PKA*ERK2
    dREGA_dt  = k7 - k8*ERK2*REGA
    dicAMP_dt = k9*ACA - k10*REGA*icAMP
    decAMP_dt = k11*ACA - k12*ecAMP
    dCAR1_dt  = k13*ecAMP - k14*CAR1
    
    dxdt = [dACA_dt,
            dPKA_dt,
            dERK2_dt,
            dREGA_dt,
            dicAMP_dt,
            decAMP_dt,
            dCAR1_dt]
    return dxdt

init_time = 0 #[min]
final_time = 1800.0 #[min]
time_interval = [init_time, final_time]

ki_para = [2.0, 0.9, 2.5, 1.5, 0.6, 0.8, 1.0, 1.3, 0.3, 0.8, 0.7, 4.9, 23.0, 4.5]

# start with a random initial condition and simiulate long enough to
# reach the stable oscillation
init_cond = np.random.rand(7)
sol_out = solve_ivp(Dicty_cAMP, (init_time, final_time), init_cond, args=(ki_para,))
tout = sol_out.t
xout = sol_out.y

# start with the initial condtion on the oscillation trajectory
# for 30 minutes
final_time = 30.0
num_data = 1000
init_cond = xout[:,-1]
t_eval = np.linspace(init_time, final_time, num_data)
sol_out = solve_ivp(Dicty_cAMP, (init_time, final_time), init_cond, t_eval=t_eval, args=(ki_para,))
tout = sol_out.t
xout = sol_out.y

import matplotlib.pyplot as plt
figure, axis = plt.subplots(1,1)
axis.plot(tout,xout[4,:])
axis.plot(tout,xout[6,:])
axis.set_ylabel('[$\mu$M]')
axis.set_xlabel('time [min]')
axis.legend(('i-cAMP','CAR 1'),loc='upper right')
