#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 23:11:34 2020

@author: jongrae kim
"""

import numpy as np
from numpy import linspace
from scipy.integrate import solve_ivp

init_time = 0 # [s]
final_time = 10.0 # [s]
num_data = 200
tout = linspace(init_time, final_time, num_data)

J_inertia = np.array([[0.005, -0.001, 0.004],
                      [-0.001, 0.006, -0.002],
                      [0.004, -0.002, 0.004]])
J_inv = np.linalg.inv(J_inertia)
J_inv_J_inertia = np.vstack((J_inertia,J_inv))

q0 = np.array([0,0,0,1])
w0 = np.array([0,0,0])

state_0 = np.hstack((q0,w0))

def dqdt_attitude_kinematics(q_true, w_true):
    quat=q_true 

    wx=np.array([[0,           -w_true[2],     w_true[1]],
              [w_true[2],   0,              -w_true[0]],
              [-w_true[1],  w_true[0],      0]])
    
    Omega_13 = np.hstack((-wx,np.resize(w_true,(3,1))))
    Omega_4  = np.hstack((-w_true,0))
    Omega = np.vstack((Omega_13, Omega_4))
     
    dqdt = 0.5*(Omega@quat)
     
    return dqdt


def dwdt_attitude_dynamics(w_true,J_inertia,inv_J_inertia, M_torque):

    Jw = J_inertia@w_true
    Jw_dot = -np.cross(w_true,Jw) + M_torque
    
    dwdt = inv_J_inertia@Jw_dot
    
    return dwdt


def dqdt_dwdt(time,state,J_inv_J_inertia):
    
    q_current = state[0:4]
    q_current = q_current/np.linalg.norm(q_current)
    w_current = state[4::]

    J_inertia = J_inv_J_inertia[0:3,:]
    J_inv = J_inv_J_inertia[3::,:]

    M_torque = np.array([0.00001+0.0005*np.sin(2*time), 
                         -0.00002+0.0001*np.cos(0.75*time), 
                         -0.0001])
    
    dqdt = dqdt_attitude_kinematics(q_current, w_current)
    dwdt = dwdt_attitude_dynamics(w_current, J_inertia, J_inv, M_torque)
    
    dstate_dt = np.hstack((dqdt,dwdt))
    return dstate_dt

sol = solve_ivp(dqdt_dwdt, (init_time, final_time), state_0, t_eval=tout, 
                rtol=1e-6,atol=1e-9, max_step=0.01, args=(J_inv_J_inertia,))
qout = sol.y[0:4,:]
wout = sol.y[4::,:]


import matplotlib.pyplot as plt

fig, (ax,ax1) = plt.subplots(nrows=2,ncols=1)
ax.plot(tout,qout[0,:],'b-',tout,qout[1,:],'r--',tout,qout[2,:],'g-.',tout,qout[3,:],'m:')

fig.set_figheight(6) # size in inches
fig.set_figwidth(8)  # size in inches

xtick_list = np.array([0,1,2,3,4,5,6,7,8,9,10])
ax.set_xticks(xtick_list)
ax.set_xticklabels(xtick_list,fontsize=14)

ytick_list = np.array([-1.0,-0.5,0.0,0.5,1.0])
ax.set_yticks(ytick_list)
ax.set_yticklabels(ytick_list,fontsize=14)

ax.legend(('$q_1$','$q_2$','$q_3$','$q_4$'),fontsize=14, loc='upper right')
ax.axis((0,10,-1.0,1.0))
ax.set_ylabel('quaternion',fontsize=14)

ax1.plot(tout,wout[0,:],'r-',tout,wout[1,:],'b--',tout,wout[2,:],'m-.')
ax1.set_xticks(xtick_list)
ax1.set_xticklabels(xtick_list,fontsize=14)
ax1.set_xlabel('time [s]',fontsize=14)

ytick_list = np.array([-2.0,-1.0,0.0,1.0,2.0])
ax1.set_yticks(ytick_list)
ax1.set_yticklabels(ytick_list,fontsize=14)

ax1.legend(('$\omega_1$','$\omega_2$','$\omega_3$'),fontsize=14, loc='upper right')
ax1.axis((0,10,-2.0,2.0))
ax1.set_xlabel('time [s]',fontsize=14)
ax1.set_ylabel('$\omega$ [rad/s]',fontsize=14)

fig.set_size_inches(9,6)    
fig.savefig('dwdt_dqdt_solve_Mt.pdf',dpi=250)