#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 23:11:34 2020

@author: menjkim
"""

import numpy as np
from numpy import linspace
from scipy.integrate import solve_ivp

init_time = 0 # [s]
final_time = 60.0 # [s]
num_data = 1000
tout = linspace(init_time, final_time, num_data)

q0 = np.array([0,0,0,1])

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


sol = solve_ivp(dqdt_attitude_kinematics, (init_time, final_time), q0, t_eval=tout)
qout = sol.y


import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot(tout,qout[0,:],'b-',tout,qout[1,:],'r--',tout,qout[2,:],'g-.',tout,qout[3,:],'m:')

fig.set_figheight(6) # size in inches
fig.set_figwidth(8)  # size in inches

xtick_list = np.array([0,10,20,30,40,50,60])
ax.set_xticks(xtick_list)
ax.set_xticklabels(xtick_list,fontsize=14)

ytick_list = np.array([-0.5,0.0,0.5,1.0])
ax.set_yticks(ytick_list)
ax.set_yticklabels(ytick_list,fontsize=14)

ax.legend(('q1','q2','q3','q4'),fontsize=14, loc='upper right')
ax.axis((0,60,-0.5,1.0))
ax.set_xlabel('time [s]',fontsize=14)
ax.set_ylabel('quaternion',fontsize=14)