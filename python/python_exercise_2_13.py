#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 21:57:46 2022

@author: jongrae
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def angular_velocity_true(time_c):
    w_true = np.array([ 0.01*np.sin(2*np.pi*0.005*time_c), # [rad/s]
                        0.05*np.cos(2*np.pi*0.001*time_c + 0.2), #[rad/s]
                        0.02 #[rad/s]
                       ])
    return w_true

def q2dcm(quat):

    quat = quat.squeeze()
    q13 =quat[0:3]
    q4 = quat[3]

    q13x = np.array([[ 0,          -q13[2],       q13[1]],
                     [q13[2],      0,           -q13[0]],
                     [-q13[1],      q13[0],       0]])
    
    q13 = q13.reshape(3,1)

    dcm = (q4**2-q13.T@q13)*np.eye(3) + 2*(q13@q13.T) - 2*q4*q13x
    return dcm

def dqdt_attitude_kinematics(time, state):
     quat = state
     w_true = angular_velocity_true(time) #[rad/s]
    
     wx=np.array([[0,           -w_true[2],     w_true[1]],
                  [w_true[2],   0,              -w_true[0]],
                  [-w_true[1],  w_true[0],      0]])
    
     Omega_13 = np.hstack((-wx,np.resize(w_true,(3,1))))
     Omega_4  = np.hstack((-w_true,0))
     Omega = np.vstack((Omega_13, Omega_4))
     
     dqdt = 0.5*(Omega@quat)
     
     return dqdt

def kalman_filter_attitude(x_hat_0, P0, dt_KF, 
                     rR_star_all, rB_star_measure, w_measure, 
                     sgm_v, sgm_u, sgm_star):
    
    num_KF_state = 6 # [dq1 dq2 dq3 b1 b2 b3]
    q_est = x_hat_0[0:4]; q_est = q_est.reshape(4,1)
    b_est = x_hat_0[4:];  b_est = b_est.reshape(3,1)
    
    w_hat = w_measure.reshape(3,1) - b_est
    w_hat_mag = np.sqrt(w_hat.T@w_hat)[0,0]
    w_hatx = np.array([[0,           -w_hat[2,0],         w_hat[1,0]],
                       [w_hat[2,0],     0,               -w_hat[0,0]],
                       [-w_hat[1,0],     w_hat[0,0],         0]])
    
    
    # propagate
    if w_hat_mag > 1e-12:
        dtheta_k = w_hat_mag*dt_KF
        cos_th = np.cos(dtheta_k/2)
        sin_th_over_w = np.sin(dtheta_k/2)/w_hat_mag
        
        q_Phi_1st_row = np.hstack((cos_th*np.eye(3)-sin_th_over_w*w_hatx,sin_th_over_w*w_hat))
        q_Phi_2nd_row = np.hstack((-sin_th_over_w*w_hat.T[0],cos_th))
        q_Phi = np.vstack((q_Phi_1st_row,q_Phi_2nd_row))
    
        cos_th = np.cos(w_hat_mag*dt_KF)
        sin_th = np.sin(w_hat_mag*dt_KF)
        Phi_1 = np.eye(3) - w_hatx*sin_th/w_hat_mag + (w_hatx@w_hatx)*((1-cos_th)/w_hat_mag**2)
        Phi_2 = -np.eye(3)*dt_KF + w_hatx*((1-cos_th)/w_hat_mag**2) - (w_hatx@w_hatx)*((w_hat_mag*dt_KF-sin_th)/w_hat_mag**3)
    else:
        q_Phi = np.vstack((np.eye(3)-(dt_KF/2)*w_hatx,-(dt_KF/2)*w_hat.T))
        Phi_1 = np.eye(3) - w_hatx*dt_KF
        Phi_2 = -np.eye(3)*dt_KF
        
    q_est_minus = q_Phi@q_est; q_est_minus = q_est_minus/np.linalg.norm(q_est_minus)
    b_est_minus = b_est
    dcm_BR_minus = q2dcm(q_est_minus)

    Q_1st_row = np.hstack(((sgm_v**2*dt_KF+(dt_KF**3/3)*sgm_u**2)*np.eye(3),-(dt_KF**2/2)*sgm_u**2*np.eye(3)-(dt_KF**3/6)*sgm_u**2*w_hatx))
    Q_2nd_row = np.hstack((-(dt_KF**2/2)*sgm_u**2*np.eye(3)-(dt_KF**3/6)*sgm_u**2*w_hatx,sgm_u**2*dt_KF*np.eye(3)))
    Q = np.vstack((Q_1st_row,Q_2nd_row))

    Phi_1st_row = np.hstack((Phi_1,Phi_2))
    Phi_2nd_row = np.hstack((np.zeros((3,3)),np.eye(3)))
    Phi = np.vstack((Phi_1st_row,Phi_2nd_row))
    
    P1 = Phi@P0@Phi.T + Q
    
    # update
    num_star = rB_star_measure.shape[1]
    rB_star_hat = dcm_BR_minus@rR_star_all
    H_k = np.zeros((3*num_star,6))
    R = sgm_star**2*np.eye(num_star*3)
    
    for xdx in range(num_star):
        vec = rB_star_hat[:,xdx]
        vec_x = np.array([[0, -vec[2], vec[1]],
                          [vec[2], 0, -vec[0]],
                          [-vec[1], vec[0], 0]])
        st_idx = 3*xdx
        H_k[st_idx:st_idx+3,:] = np.hstack((vec_x,np.zeros((3,3))))
    
    K_k = P1@H_k.T@np.linalg.inv(H_k@P1@H_k.T+R)
    P1 = (np.eye(num_KF_state)-K_k@H_k)@P1
    
    rB_star_mea_vec = rB_star_measure.T.reshape(3*num_star,1)
    rB_star_hat_vec = rB_star_hat.T.reshape(3*num_star,1)
    
    delta_x = K_k@(rB_star_mea_vec-rB_star_hat_vec)
    
    # quaternion & bias update
    dq_13 = 2*delta_x[0:3]; dq_13 = dq_13.reshape((3,1))
    q = q_est_minus.squeeze()
    qx = np.array([[0, -q[2], q[1]],
                   [q[2],0,-q[0]],
                   [-q[1],q[0],0]])
    
    quat_update_matrix = np.vstack((q[3]*np.eye(3)+qx,-q[0:3].T))

    q_hat_plus = q_est_minus + quat_update_matrix@dq_13; q_hat_plus = q_hat_plus/np.linalg.norm(q_hat_plus)
    b_hat_plus = b_est_minus + delta_x[3:]
    
    x_hat_1 = np.vstack((q_hat_plus,b_hat_plus))
    
    return x_hat_1, P1

#----------------------------------------------------------------------------
# Set initial values & change non-SI units into the SI Units
dt = 0.05 # [seconds]
time_init = 0
time_final = 120 # [seconds]
N_sample = int(time_final/dt) + 1
time = np.linspace(time_init,time_final, N_sample)

# standard deviation of the bias, sigma_beta_xyz
sigma_beta = 0.0005 # [degrees/sqrt(s)]
sigma_u = sigma_beta*(np.pi/180) # [rad/sqrt(s)]
sigma_eta = sigma_u/np.sqrt(dt)

# standard devitation of the white noise, sigma_v
sigma_v = 0.0001 #[degrees/s]
sigma_v = sigma_v*(np.pi/180) #[rad/s]

# initial beta(t)
beta = (2*np.random.rand(3)-1)*0.03 # +/- 0.03[degrees/s]
beta = beta*(np.pi/180) # [radians/s]

# prepare the data store
w_all = np.zeros((N_sample,3))
w_measure_all = np.zeros((N_sample,3))

# data store
# instead of calculating the exact size
# of the following matrices, use varying matrices with increasing
# time, which might not be significant but simpler to implement
w_gyr_all = []
w_hat_all = []
w_tr_all = []

q_tr_all = []
q_hat_all = []
time_all = []
pcov_all = []

q_current = np.array([0,0,0,1])

# star sensor
# star sensor reference star vectors
r1R = np.array([-0.6794, -0.3237, -0.6586]).reshape((3,1))
r2R = np.array([-0.7296,  0.5858,  0.3528]).reshape((3,1))
r3R = np.array([-0.2718,  0.6690, -0.6918]).reshape((3,1))
r4R = np.array([-0.2062, -0.3986, 0.8936]).reshape((3,1))
r5R = np.array([0.6858, -0.7274, -0.0238]).reshape((3,1))

r1R = r1R/np.sqrt(r1R.T@r1R)[0,0]
r2R = r2R/np.sqrt(r2R.T@r2R)[0,0]
r3R = r3R/np.sqrt(r3R.T@r3R)[0,0]
r4R = r4R/np.sqrt(r4R.T@r4R)[0,0]
r5R = r5R/np.sqrt(r5R.T@r5R)[0,0]

rR_star_all = np.hstack((r1R,r2R,r3R,r4R,r5R))
num_star = rR_star_all.shape[1]
sigma_star = 87.2665/3*1e-6
r_star = sigma_star**2*np.eye(num_star*3)

# Kalman filter
n_dt_KF = 2
dt_KF = n_dt_KF*dt
bias_estimate_current = np.zeros((3,1))

q_estimate_current = np.array([0, 0, 0, 1]).reshape((4,1)) + 0.0*np.random.randn(4,1)
q_estimate_current = q_estimate_current/np.linalg.norm(q_estimate_current)
x0 = np.vstack((q_estimate_current, bias_estimate_current))
p_current = 0.001*np.eye(6)
w_hat = np.array([0, 0, 0]).reshape((3,1))

# main simulation loops
for idx in range(N_sample):
    
    time_c = time[idx]
    w_true = angular_velocity_true(time_c)
    
    # beta(t)
    eta_u = sigma_eta*np.random.randn(3)
    dbeta = eta_u*dt
    beta = beta + dbeta
    
    # eta_v(t)
    eta_v = sigma_v*np.random.randn(3)
    
    # w_tilde
    w_measurement = w_true + beta + eta_v
    
    
    if np.remainder(idx,n_dt_KF)==1:
        
        # star sensor measurement
        dcm_BR = q2dcm(q_current)
        rB_star_all = dcm_BR@rR_star_all
        rB_star_measure = rB_star_all+sigma_star*np.random.randn(3,num_star)
        rB_star_measure = rB_star_measure/np.kron(np.ones((3,1)),np.sqrt(np.sum(rB_star_measure**2,0)))
        
        # kalman filter
        x_hat_1, P1 = kalman_filter_attitude(x0, p_current, dt_KF,
                            rR_star_all, rB_star_measure, w_measurement, 
                            sigma_v, sigma_u, sigma_star)
        
        x0 = x_hat_1
        p_current = P1
        
        q_estimate_current = x0[0:4]
        q_estimate_current = q_estimate_current/np.linalg.norm(q_estimate_current)  
     
        bias_estimate_current = x0[4:]
        
        w_hat = w_measurement.reshape((3,1)) - bias_estimate_current
        
        # store data to plot
        # instead of calculating the exact size
        # of the following matrices, use varying matrices with increasing
        # time, which might not be significant but simpler to implement
        time_all.append(time_c)
        w_gyr_all.append(w_measurement.squeeze())
        w_hat_all.append(w_hat.squeeze())
        w_tr_all.append(w_true.squeeze())
        
        q_tr_all.append(q_current.squeeze())
        q_hat_all.append(q_estimate_current.squeeze())
        pcov_all.append(np.diag(P1).squeeze())
        
    # integrate true dqdt to obtain true q(t): time_c -> time_c + dt
    if idx < N_sample-1:
        sol = solve_ivp(dqdt_attitude_kinematics, (time_c, time[idx+1]), q_current)
        q_current = sol.y[:,-1]


# make the lists to numpy arrays
q_tr_all=np.array(q_tr_all)
q_hat_all=np.array(q_hat_all)
dq = q_tr_all-q_hat_all

w_tr_all=np.array(w_tr_all)
w_hat_all=np.array(w_hat_all)
dw = w_tr_all - w_hat_all

pcov_all = np.array(pcov_all)

# plot results
fig0, (ax0,ax1,ax2) = plt.subplots(nrows=3,ncols=1)
ax0.plot(time_all,dq[:,0],time_all,-3*np.sqrt(pcov_all[:,0]),'r--',time_all,3*np.sqrt(pcov_all[:,0]),'r--')
ax0.axis((time_init,time_final,-4e-5,4e-5))
ax0.set_ylabel('$\delta q_1$',fontsize=14)
ax0.legend(('error','3$\sigma$ bound'),fontsize=8, loc='upper right')
ax1.plot(time_all,dq[:,1],time_all,-3*np.sqrt(pcov_all[:,1]),'r--',time_all,3*np.sqrt(pcov_all[:,1]),'r--')
ax1.axis((time_init,time_final,-4e-5,4e-5))
ax1.set_ylabel('$\delta q_2$',fontsize=14)
ax1.legend(('error','3$\sigma$ bound'),fontsize=8, loc='upper right')
ax2.plot(time_all,dq[:,2],time_all,-3*np.sqrt(pcov_all[:,2]),'r--',time_all,3*np.sqrt(pcov_all[:,2]),'r--')
ax2.axis((time_init,time_final,-4e-5,4e-5))
ax2.set_ylabel('$\delta q_3$',fontsize=14)
ax2.legend(('error','3$\sigma$ bound'),fontsize=8, loc='upper right')
ax2.set_xlabel('time [s]')


fig1, (bx0,bx1,bx2) = plt.subplots(nrows=3,ncols=1)
bx0.plot(time_all,dw[:,0],time_all,-3*np.sqrt(pcov_all[:,3]),'r--',time_all,3*np.sqrt(pcov_all[:,3]),'r--')
bx0.axis((time_init,time_final,-4e-5,4e-5))
bx0.set_ylabel('$\delta\omega_1$ [rad/s]',fontsize=14)
bx0.legend(('error','3$\sigma$ bound'),fontsize=8, loc='upper right')
bx1.plot(time_all,dw[:,1],time_all,-3*np.sqrt(pcov_all[:,4]),'r--',time_all,3*np.sqrt(pcov_all[:,4]),'r--')
bx1.axis((time_init,time_final,-4e-5,4e-5))
bx1.set_ylabel('$\delta\omega_2$ [rad/s]',fontsize=14)
bx1.legend(('error','3$\sigma$ bound'),fontsize=8, loc='upper right')
bx2.plot(time_all,dw[:,2],time_all,-3*np.sqrt(pcov_all[:,5]),'r--',time_all,3*np.sqrt(pcov_all[:,5]),'r--')
bx2.axis((time_init,time_final,-4e-5,4e-5))
bx2.set_ylabel('$\delta\omega_3$ [rad/s]',fontsize=14)
bx2.legend(('error','3$\sigma$ bound'),fontsize=8, loc='upper right')
bx2.set_xlabel('time [s]')