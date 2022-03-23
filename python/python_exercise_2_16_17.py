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
final_time = 10.0 # [s]
num_data = 200
tout = linspace(init_time, final_time, num_data)

J_inertia = np.array([[0.005, -0.001, 0.004],
                      [-0.001, 0.006, -0.002],
                      [0.004, -0.002, 0.004]])
J_inv = np.linalg.inv(J_inertia)
J_inv_J_inertia = np.vstack((J_inertia,J_inv))

C_T = 8.8e-7 # motor thruster coefficient [N/(rad/s)^2]
C_D = 11.3e-8 # motor drag coefficient [Nm/(rad/s)^2]
L_arm = 0.127 # length from the centre of quadcopter to the motor [m]

Motor_para = np.array([C_T,C_D,L_arm])

q0 = np.array([0,0,0,1])
w0 = np.array([0,0,0])

#----------------------Added Part---------------------------------------
# first-order motor model
tau_motor =  0.01#[s]
w_motor = np.array([0,0,0,0]) # [initial 4-motor angular velocity]
#-----------------------------------------------------------------------

#----------------------Modified Part------------------------------------
state_0 = np.hstack((q0,w0,w_motor))
#----------------------Added Part---------------------------------------

# use global variables only for saving values
global global_motor_time_FM_all

# minimum time interval for saving values to the global
dt_save = 0.05

def propeller_motor_actuator(C_T,C_D,L_arm,w_motor):

    # this part is replaced by current motor angular velocity
    # assume perfect motor angular velocity control
    # w_motor = w_command

    F_fblr = C_T*(w_motor**2)
    tau_fblr = C_D*(w_motor**2)
    
    F_motor = -np.sum(F_fblr)
    M_motor = np.array([ L_arm*(F_fblr[2]-F_fblr[3]),
                L_arm*(F_fblr[0]-F_fblr[1]),
                np.sum(tau_fblr[0:2])-np.sum(tau_fblr[2::])])
            
    FM_Motor = np.hstack([F_motor, M_motor])
    return FM_Motor

def propeller_motor_FM2w_conversion(F_M_Desired, C_T, C_D, L_arm):

    inv_C_T = 1/C_T;
    inv_C_D = 1/C_D;
    inv_2_L_C_T = 2/(L_arm*C_T);
    
    Conv_Mat = 0.25*np.array([   
                        [-inv_C_T,     0,               inv_2_L_C_T,    inv_C_D],
                        [-inv_C_T,     0,               -inv_2_L_C_T,     inv_C_D],
                        [-inv_C_T,     inv_2_L_C_T,     0,               -inv_C_D],
                        [-inv_C_T,     -inv_2_L_C_T,    0,               -inv_C_D]])
    
    w_motor_fblr_squared_desired = Conv_Mat@F_M_Desired
    return w_motor_fblr_squared_desired

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

#------------------Updated Part: tau_motor added in the input parameters--
def dqdt_dwdt(time,state,J_inv_J_inertia,Motor_para,dt_save, tau_motor):
#-------------------------------------------------------------------------
    
    global global_motor_time_FM_all
    
    q_current = state[0:4]
    q_current = q_current/np.linalg.norm(q_current)
    w_current = state[4:7]
    
    #------------------Updated Part: motor state---
    w_motor_current = state[7::]
    #---------------------------------------------

    J_inertia = J_inv_J_inertia[0:3,:]
    J_inv = J_inv_J_inertia[3::,:]
    C_T = Motor_para[0]
    C_D = Motor_para[1]
    L_arm = Motor_para[2]

    #--------------------------------
    # Begin: this part is controller
    #--------------------------------
    M_Desired = np.array([0.00001+0.0005*np.sin(2*time), 
                         -0.00002+0.0001*np.cos(0.75*time), 
                         -0.0001])
    mg = 10.0 #[N]
    F_M_Desired = np.hstack([-mg, M_Desired])
    
    w_motor_fblr_squared_desired = propeller_motor_FM2w_conversion(F_M_Desired, 
                                                                   C_T, C_D, L_arm)
    w_motor_fblr_desired = np.sqrt(w_motor_fblr_squared_desired)
    #--------------------------------
    # End: this part is controller
    #--------------------------------

    #------------------------added part----------------------------------
    dw_motor_dt = (-w_motor_current + w_motor_fblr_desired)/tau_motor

    # Motor Force & Torque
    # desired motor angular velocity is replaced with current motor angular velocity
    FM_Motor = propeller_motor_actuator(C_T, C_D, L_arm, w_motor_current) # updated part
    M_torque = FM_Motor[1::]
    
    current_data = np.hstack((time,FM_Motor,))
    if time < 1e-200:
        global_motor_time_FM_all = current_data.reshape(1,5)
    elif time > global_motor_time_FM_all[-1,0]+dt_save:
        global_motor_time_FM_all = np.vstack((global_motor_time_FM_all,current_data,))
    
    dqdt = dqdt_attitude_kinematics(q_current, w_current)
    dwdt = dwdt_attitude_dynamics(w_current, J_inertia, J_inv, M_torque)
    
    dstate_dt = np.hstack((dqdt,dwdt,dw_motor_dt)) # updated part
    return dstate_dt

# solve ode
sol = solve_ivp(dqdt_dwdt, (init_time, final_time), state_0, t_eval=tout, 
                atol=1e-9, rtol=1e-6, max_step=0.01, 
                args=(J_inv_J_inertia, Motor_para, dt_save,tau_motor))
qout = sol.y[0:4,:]
wout = sol.y[4:7,:]
w_motor_out = sol.y[7::,:]

time_Motor = global_motor_time_FM_all[:,0]
Force_Motor = global_motor_time_FM_all[:,1]
Torque_Motor = global_motor_time_FM_all[:,2::]
del global_motor_time_FM_all

import matplotlib.pyplot as plt

# quaternion and omega plot

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

# motor force and torque plot
fig1, (ax2,ax3) = plt.subplots(nrows=2,ncols=1)
ax2.plot(time_Motor,Force_Motor)
ax2.set_ylabel('Force [N]',fontsize=14)
ax2.set_xticks(xtick_list)
ax2.set_xticklabels(xtick_list,fontsize=14)

ytick_list = np.array([-15.0,-10.0,-5.0,0.0])
ax2.set_yticks(ytick_list)
ax2.set_yticklabels(ytick_list,fontsize=14)
ax2.set_xlim([init_time, final_time])
ax2.set_ylim([-15,0])

ax3.plot(time_Motor,Torque_Motor[:,0],'r-',time_Motor,Torque_Motor[:,1],'b--',
     time_Motor,Torque_Motor[:,2],'m-.')
ax3.set_xticks(xtick_list)
ax3.set_xticklabels(xtick_list,fontsize=14)

ytick_list = np.array([-0.0006, -0.0003, 0.0, 0.0003, 0.0006])
ax3.set_yticks(ytick_list)
ax3.set_yticklabels(ytick_list,fontsize=14)
ax3.set_xlim([init_time, final_time])
ax3.set_ylim([-0.0006,0.0006])
ax3.set_ylabel('Torque [Nm]',fontsize=14)
ax3.set_xlabel('time [s]',fontsize=14)
ax3.legend(('$M_1$','$M_2$','$M_3$'),fontsize=14)

# motor angular velocity
fig2, ax1 = plt.subplots(nrows=1,ncols=1)
ax1.plot(tout, w_motor_out.T*60/(2*np.pi))
ax1.set_ylabel('Motor Angular Velocity [rpm]',fontsize=14)
ax1.set_xlabel('time [s]',fontsize=14)
ax1.legend(('$\omega_m^{forward}$','$\omega_m^{backward}$','$\omega_m^{left}$','$\omega_m^{right}$'))

#fig1.set_size_inches(9,6)    
#fig1.savefig('../figures/quad_motor_Force_Torque.pdf',dpi=250)
