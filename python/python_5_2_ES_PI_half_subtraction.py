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

# parameters
kP = 50
kI = 5e-6
gamma_G = 8e-4
ks2 = 4e-4

#----------------- Change Here------------------------------
Kdeg = 1e-3   #1e-3    # PI-Deg gain
X3 = 1 
KF = 3  # 0.62; 3      # Pd scaling

eta = 1e-4             # integrator self-annihilation  

Pd = 1                     # Desired P
time_for_step_change = 8   # [hours]
#-----------------------------------------------------------

kon = 5e-5
kcat = 1.6*2
kdg = 8e-8
ks1 = 3

E_total = 0.2

# simulation time values
time_current = 0       # initial time 
time_final   = 3600*16 # final time [min]
tspan = [time_current, time_final]

# E-S PI Control Half Subtraction
def ES_PI_Half_Subtraction(time,state,ki_para):
    E, P, ES, X1, X2, S,DP, Xs   = state
    
    kP, kI, gamma_G, ks2, kon, kcat, kdg, ks1, Pd, KF, Kdeg, X3, eta, time_for_step_change = ki_para
    
    if time > 3600*time_for_step_change:
        Pd = Pd/2
        
    Pd = Pd*KF
    
    dE_dt = -kon*E*S + kcat*ES;
    dP_dt = kcat*ES - kdg*P - Kdeg*X3*P;
    dES_dt = kon*E*S - kcat*ES;
    dX1_dt = kI*DP - eta*X1;
    dX2_dt = -gamma_G*X2 + gamma_G*kP*DP;
    dS_dt = ks2*X1 + ks2*X2 - ks2*S;
    dDP_dt = ks1*Pd - ks1*DP*Xs - ks1*DP;
    dXs_dt = -ks1*DP*Xs + ks1*P;
    
    dxdt = [dE_dt, dP_dt, dES_dt, dX1_dt, dX2_dt, dS_dt, dDP_dt, dXs_dt]
    return dxdt


# simulation
para = [kP, kI, gamma_G, ks2, kon, kcat, kdg, ks1, Pd, KF, Kdeg, X3, eta, time_for_step_change]
state_t0    = 0.1*np.ones(8) 
state_t0[2] = E_total - state_t0[0]
state_t0[3] = 1e-3 
state_t0[4] = 1e-3

sol = solve_ivp(ES_PI_Half_Subtraction, (time_current,time_final), state_t0, args=(para,))

tout = sol.t
xout = sol.y

import matplotlib.pyplot as plt

figure, (ax0,ax1) = plt.subplots(2,1)
time_hr = tout/3600 # [hour]
P_history = xout[1,:]

ax0.plot([0, time_for_step_change,time_for_step_change, time_final/3600],
         [Pd, Pd, 0.5*Pd, 0.5*Pd], 'r--')
ax0.plot(time_hr, P_history,'b')
ax0.set_ylabel('[P(t)] [a.u.]')
ax0.legend(('desired [P]','achieved [P]'))
ax0.axis([0,time_final/3600,0,1.2])

ax1.plot(time_hr, xout[5,:])
ax1.set_xlabel('time [hour]')
ax1.set_ylabel('[S] [a.u.]')
ax1.axis([0, time_final/3600, 0, 120])