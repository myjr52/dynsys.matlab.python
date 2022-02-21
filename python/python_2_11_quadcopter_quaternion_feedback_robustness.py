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
final_time = 120.0 # [s]
num_data = 1200
tout = linspace(init_time, final_time, num_data)

J_inertia = np.array([[0.005, -0.001, 0.004],
                      [-0.001, 0.006, -0.002],
                      [0.004, -0.002, 0.004]])
J_inv = np.linalg.inv(J_inertia)

C_T = 8.8e-7 # motor thruster coefficient [N/(rad/s)^2]
C_D = 11.3e-8 # motor drag coefficient [Nm/(rad/s)^2]
L_arm = 0.127 # length from the centre of quadcopter to the motor [m]

q0 = np.array([1,1,-1,1])/np.sqrt(4)
w0 = np.array([0.1,-0.2,0.1])

r0 = np.array([0,0,-50])
v0 = np.array([0,0,0])

state_0 = np.hstack((q0,w0,r0,v0))

num_MC = 300
ts_all = np.inf*np.ones(num_MC)
dJ_norm_all = np.inf*np.ones(num_MC)

#-----------------------------------------------------------------------------
# propellar & motor
#-----------------------------------------------------------------------------
def propeller_motor_actuator(C_T,C_D,L_arm,w_command):

    # assume perfect motor angular velocity control
    w_motor = w_command;

    F_fblr = C_T*(w_motor**2)
    tau_fblr = C_D*(w_motor**2)
    
    F_motor = -np.sum(F_fblr)
    M_motor = np.array([ L_arm*(F_fblr[2]-F_fblr[3]),
                L_arm*(F_fblr[0]-F_fblr[1]),
                np.sum(tau_fblr[0:2])-np.sum(tau_fblr[2::])])
            
    FM_w_Motor = np.hstack([F_motor, M_motor, w_motor])
    return FM_w_Motor

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
    w_motor_fblr_squared_desired[w_motor_fblr_squared_desired<0] = 0.0
    return w_motor_fblr_squared_desired


#-----------------------------------------------------------------------------
# quaternion feedback & altitude control
#-----------------------------------------------------------------------------
def quaternion_feedback_and_altitude_control(q_current, w_current, rv_current, 
                                             J_inertia, C_T, C_D, L_arm, C_BR, 
                                             mass_quadcopter, grv_acce):
    
    zR_desired = -30 #[m]
    zdotR_desired = 0 #[m/s]
    K_qf = 0.01*np.eye(3);
    C_qf = 0.001*np.eye(3)
    k1 = 0.1
    k2 = 0.5
    
    q_13 = q_current[0:3]
    w = w_current
    
    Fmg_R = grv_acce*mass_quadcopter #[N]
    Falt_R = k1*(zR_desired-rv_current[2])+k2*(zdotR_desired-rv_current[5])
    F_desired_R = np.array([0,0,-Fmg_R+Falt_R])
    F_desired_B = C_BR@F_desired_R;
    
    u_qf = -K_qf@q_13 - C_qf@w - np.cross(w,J_inertia@w);
    M_Desired = u_qf;
    
    F_M_desired = np.hstack((F_desired_B[2], M_Desired))
                
    w_motor_fblr_squared_desired = propeller_motor_FM2w_conversion(F_M_desired, C_T, C_D, L_arm)
   
    w_motor_fblr_desired = np.sqrt(w_motor_fblr_squared_desired)

    return w_motor_fblr_desired 


#-----------------------------------------------------------------------------
# translational & rotational kinematics & dynamics
#-----------------------------------------------------------------------------
def drdt_linear_dynamics(rv_true, mass, grv_const, motor_force_in_R):
    v = rv_true[3::]
    
    drdt = v
    dvdt = np.array([0,0,grv_const]) + motor_force_in_R/mass
    
    drv_dt = np.hstack((drdt, dvdt))
    
    return drv_dt

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

def dqdt_dwdt_drvdt(time,state,quadcopter_uav):
    
    m_quadcopter = 0.49 #[kg]
    grv_acce = 9.81 #[m/s^2]
    
    q_current = state[0:4]
    q_current = q_current/np.linalg.norm(q_current)
    
    w_current = state[4:7]
    
    rv_current = state[7::]

    J_inertia = quadcopter_uav[0]
    J_inv = quadcopter_uav[1]
    C_T = quadcopter_uav[2]
    C_D = quadcopter_uav[3]
    L_arm = quadcopter_uav[4]
    
    q_13 = q_current[0:3]
    q4 = q_current[3]
    q13x = np.array([[0, -q_13[2], q_13[1]],
                     [q_13[2], 0, -q_13[0]],
                     [-q_13[1], q_13[0], 0]])
    
    C_BR = (q4**2-q_13@q_13)*np.eye(3) + 2*(q_13.reshape((3,1))*q_13.reshape((1,3)))-2*q4*q13x
    
    #--------------------------------
    # Begin: this part is controller
    #--------------------------------
    w_motor_fblr_desired = quaternion_feedback_and_altitude_control(q_current,
        w_current, rv_current, J_inertia, C_T, C_D, L_arm, C_BR, m_quadcopter, grv_acce)
    #--------------------------------
    # End: this part is controller
    #--------------------------------

    # Motor Force & Torque
    FM_w_Motor = propeller_motor_actuator(C_T, C_D, L_arm, w_motor_fblr_desired)
    
    F_motor = np.array([0,0,FM_w_Motor[0]])
    M_torque = FM_w_Motor[1:4]
    
    motor_force_in_R = C_BR.transpose()@F_motor
    
    dqdt = dqdt_attitude_kinematics(q_current, w_current)
    dwdt = dwdt_attitude_dynamics(w_current, J_inertia, J_inv, M_torque)
    drvdt = drdt_linear_dynamics(rv_current, m_quadcopter, grv_acce, motor_force_in_R)
    
    dstate_dt = np.hstack((dqdt,dwdt,drvdt))
    return dstate_dt


#-----------------------------------------------------------------------------
# Monte-Carlo Simulation
#-----------------------------------------------------------------------------
import time
start = time.time()

for g_MC_idx in range(num_MC):
    
    not_find_dJ = True
    
    while not_find_dJ:
        
        dJ = np.diag(0.002*np.random.randn(3))
    
        J_inertia_perturbed = J_inertia + dJ
      
        pd_cond = np.min(np.linalg.eig(J_inertia_perturbed)[0])>0
        j3_cond = J_inertia_perturbed[0,0]+J_inertia_perturbed[1,1] > J_inertia_perturbed[2,2]
        j2_cond = J_inertia_perturbed[0,0]+J_inertia_perturbed[2,2] > J_inertia_perturbed[1,1]
        j1_cond = J_inertia_perturbed[1,1]+J_inertia_perturbed[2,2] > J_inertia_perturbed[0,0]

        if pd_cond and j1_cond and j2_cond and j3_cond:
            not_find_dJ = False

    dJ_norm_all[g_MC_idx] = np.linalg.norm(dJ)
    J_inv_perturbed = np.linalg.inv(J_inertia_perturbed)
    
    quadcopter_uav=(J_inertia_perturbed, J_inv_perturbed, C_T, C_D, L_arm)
    
    sol = solve_ivp(dqdt_dwdt_drvdt, (init_time, final_time), state_0, t_eval=tout, 
                    atol=1e-9, rtol=1e-6, max_step=0.01, 
                    args=(quadcopter_uav,))
    qout = sol.y[0:4,:]
    wout = sol.y[4:7,:]
    rout = sol.y[7:10,:]
    vout = sol.y[10::,:]

    q13=qout[0:3,:]
    q13_norm = np.sqrt((np.sum(q13**2,axis=0)))
    q13_ts = (q13_norm>0.01)*np.ones(num_data)
    q13_ts = np.cumsum(q13_ts)
    q13_ts = tout[q13_ts==q13_ts[-1]]
    q13_ts = q13_ts[0]
    ts_all[g_MC_idx] = q13_ts
    
    print(f'#{g_MC_idx:1d} {np.linalg.norm(dJ):6.5f}, {q13_ts:4.2f}\n')


end = time.time()
print(f'{end-start} [seconds]')
#-----------------------------------------------------------------------------
# Plot Results
#-----------------------------------------------------------------------------
import matplotlib.pyplot as plt

fig0, (ax0) = plt.subplots(nrows=1,ncols=1)
ax0.plot(dJ_norm_all*1e3, ts_all,'.')
ax0.set_ylabel(r'$t_s$ [s]',fontsize=14)
ax0.set_xlabel(r'$||\Delta J||\times 10^{-3}$ [kg m$^2$]',fontsize=14)

fig0.set_figheight(6) # size in inches
fig0.set_figwidth(8)  # size in inches

# xtick_list = np.array([0,2,4,6,8,10,12])*10
# ax.set_xticks(xtick_list)
# ax.set_xticklabels(xtick_list,fontsize=14)

# ytick_list = np.array([-1.0,0.0,1.0])
# ax.set_yticks(ytick_list)
# ax.set_yticklabels(ytick_list,fontsize=14)

# ax.legend(('$q_1$','$q_2$','$q_3$','$q_4$'),fontsize=14, loc='upper right')
# ax.axis((init_time,final_time,-1.1,1.1))
# ax.set_ylabel('quaternion',fontsize=14)

# ax1.plot(tout,wout[0,:],'r-',tout,wout[1,:],'b--',tout,wout[2,:],'m-.')
# ax1.set_xticks(xtick_list)
# ax1.set_xticklabels(xtick_list,fontsize=14)
# ax1.set_xlabel('time [s]',fontsize=14)

# ytick_list = np.array([-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0])
# ax1.set_yticks(ytick_list)
# ax1.set_yticklabels(ytick_list,fontsize=14)

# ax1.legend(('$\omega_1$','$\omega_2$','$\omega_3$'),fontsize=14, loc='upper right')
# ax1.axis((init_time,final_time,-3,3))
# ax1.set_xlabel('time [s]',fontsize=14)
# ax1.set_ylabel('$\omega$ [rad/s]',fontsize=14)











