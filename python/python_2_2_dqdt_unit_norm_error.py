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
from scipy.integrate import solve_ivp

init_time = 0 # [s]
final_time = 60.0 # [s]
num_data = 1000
tout = linspace(init_time, final_time, num_data)

q0 = np.array([0,0,0,1])
RelTol_cases = np.array([1e-3, 1e-6, 1e-9])
N_sim = len(RelTol_cases)

def dqdt_attitude_kinematics(time, state):
     quat = state
     w_true = np.array([0.1*np.sin(2*np.pi*0.005*time), #[rad/s]
                        0.05*np.cos(2*np.pi*0.01*time + 0.2), #[rad/s]
                        0.02]) #[rad/s]
    
     wx=np.array([[0,           -w_true[2],     w_true[1]],
                  [w_true[2],   0,              -w_true[0]],
                  [-w_true[1],  w_true[0],      0]])
    
     Omega_13 = np.hstack((-wx,np.resize(w_true,(3,1))))
     Omega_4  = np.hstack((-w_true,0))
     Omega = np.vstack((Omega_13, Omega_4))
     
     dqdt = 0.5*(Omega@quat)
     
     return dqdt

err_idx = np.zeros((num_data,N_sim))
for idx in range(N_sim):
    sol = solve_ivp(dqdt_attitude_kinematics, (init_time, final_time), q0, t_eval=tout, rtol=RelTol_cases[idx], atol=RelTol_cases[idx]*1e-3)
    qout = sol.y
    err_idx[:,idx]=np.abs(np.sum(qout**2,0)-1)

import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot(tout,np.log(err_idx[:,0]+1e-100),'b-',
        tout,np.log(err_idx[:,1]+1e-100),'r--',
        tout,np.log(err_idx[:,2]+1e-100),'g-.')

fig.set_figheight(6) # size in inches
fig.set_figwidth(8)  # size in inches

xtick_list = np.array([0,10,20,30,40,50,60])
ax.set_xticks(xtick_list)
ax.set_xticklabels(xtick_list,fontsize=14)

ytick_list = np.array([-35,-30,-25,-20,-15,-10,-5,0])
ax.set_yticks(ytick_list)
ax.set_yticklabels(ytick_list,fontsize=14)

ax.axis((0,60,-35,0))
ax.set_xlabel('time [s]',fontsize=14)
ax.set_ylabel(r'$\log|{\bf q}^T{\bf q}-1|$',fontsize=14)

ax.legend(('RelTol=1e-3','RelTol=1e-6','RelTol=1e-9'),fontsize=14, loc='lower right')
