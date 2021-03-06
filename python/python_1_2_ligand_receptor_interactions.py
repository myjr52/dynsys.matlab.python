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

from numpy import linspace
from scipy.integrate import solve_ivp


init_receptor = 0.01 #[nM]
init_ligand = 0.0415 #[nM]
init_complex = 0.0 #[kg]

init_time = 0 #[min]
final_time = 180.0 #[min]
time_interval = [init_time, final_time]

kon = 0.0972 #[1/(min nM)] 
koff = 0.24 #[1/min]
kt = 0.02 #[1/min]
ke = 0.15 #[1/min]

ft = 0.0 #[nM/min]
QR = 0.0166 #[nM/min]
R_max = 0.415 #[nM]

sim_para = [kon, koff, kt, ke, ft, QR, R_max]

init_cond = [init_receptor, init_ligand, init_complex]


num_data = int(final_time*10)
tout = linspace(init_time, final_time, num_data)


def RLC_kinetics(time,state,sim_para):
    R, L, C = state

    kon, koff, kt, ke, ft, QR, R_max = sim_para
    
    if R > R_max:
        QR = 0

    dxdt = [-kon*R*L + koff*C - kt*R + QR,
            -kon*R*L + koff*C + ft,
            kon*R*L - koff*C - ke*C]
    return dxdt

sol_out = solve_ivp(RLC_kinetics, (init_time, final_time), init_cond, args=(sim_para,))

tout = sol_out.t
xout = sol_out.y

from scipy.integrate import odeint
xout_odeint = odeint(RLC_kinetics, init_cond, linspace(init_time, final_time, num_data), args=(sim_para,),tfirst=True)

import matplotlib.pyplot as plt
plt.figure(1)
plt.plot(tout,xout[0,:])
plt.ylabel('Receptor [nM]')
plt.xlabel('time [min]')
plt.axis([0, final_time, 0, 0.5])

plt.figure(2)
plt.plot(tout,xout[1,:])
plt.ylabel('Ligand [nM]')
plt.xlabel('time [min]')
plt.axis([0, final_time, 0, 0.05])

plt.figure(3)
plt.plot(tout,xout[2,:])
plt.ylabel('Complex [nM]')
plt.xlabel('time [min]')
plt.axis([0, final_time, 0, 0.004])
