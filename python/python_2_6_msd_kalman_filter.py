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

# Simulation setting

m_mass = 1.0 #[kg]
k_spring = 0.5 #[N/m]
c_damper = 0.01 #[N/(m/s)]
msd_const = [m_mass, k_spring, c_damper]

init_pos = 0.0 #[m]
init_vel = 0.0 #[m/s]

init_time = 0 #[s]
final_time = 60 #[s]

Delta_t = 0.01 #[s]

time_interval = [init_time, final_time]

# process noise
num_w = int((final_time-init_time)/Delta_t)+1
sigma_beta = np.sqrt(0.5)
sigma_w =sigma_beta/np.sqrt(Delta_t)
wk_noise = sigma_w*(np.random.randn(num_w))

# sensor
sigma_v = 0.75
sample_step = 20
current_step = 0
Delta_tk_KF = Delta_t*sample_step
zk_all = 1

# kalman filter
A_KF = np.array([[1, Delta_tk_KF], 
        [-(k_spring*Delta_tk_KF)/m_mass, 1-(c_damper*Delta_tk_KF)/m_mass]])
H_KF = np.array([[1, 0]])

# kalman filter initialize
Q_KF = np.array([[0, 0], [0, (Delta_tk_KF*sigma_w)**2]])
R_KF = sigma_v**2
x_hat_plus = np.array([[0],[0]])
P_plus_KF = 0.1*np.eye(2)
x_hat_all = x_hat_plus.transpose()

# Data
x0 = np.array([init_pos, init_vel])
t0 = init_time
tf = t0 + Delta_t

tout_all = np.zeros(num_w)
xout_all = np.zeros((num_w,2))

tout_all[0] = t0
xout_all[0] = x0
xtrue_sample = x0
P_all = P_plus_KF.ravel()

# system
def msd_noisy(time,state,wk,msd_const):
    x1, x2 = state
    m, k, c = msd_const
    
    dxdt = [x2,
            -(k/m)*x1 - (c/m)*x2 + wk]
    return dxdt

# simulation
for idx in range(1,num_w):
    
    wk = wk_noise[idx]

    # RK45
    sol = solve_ivp(msd_noisy, (t0, tf), x0, args=(wk,msd_const))
    xout = sol.y.transpose()    

    tout_all[idx] = sol.t[-1]
    xout_all[idx] = xout[-1]
    
    x0 = xout[-1];
    
    #  Kalman Filter
    if current_step >= sample_step:
        
        # Prediction
        x_hat_minus = A_KF@x_hat_plus
        P_minus_KF = (A_KF@P_plus_KF)@A_KF.transpose() + Q_KF
        
        # Sensor measurement
        z_k = x0[0] + sigma_v*np.random.randn()
        current_step = 1
        zk_all = np.hstack((zk_all,z_k))
        
        # Update
        K_KF = P_minus_KF@H_KF.transpose()/(H_KF@P_plus_KF@H_KF.transpose()+R_KF);
        x_hat_plus = x_hat_minus + K_KF@(z_k - H_KF@x_hat_minus)
        P_plus_KF = (np.eye(2)-K_KF@H_KF)@P_minus_KF
        
        # Collect Kalman Filter Data
        x_hat_all = np.vstack([x_hat_all,x_hat_plus.transpose()])
        P_all = np.vstack((P_all,P_plus_KF.ravel()))
        
        # Sample the true state to compare with the estimated
        xtrue_sample = np.vstack([xtrue_sample,x0.transpose()])
        
    else:
        current_step = current_step + 1
    
    # time interval update
    t0 = tf
    tf = t0 + Delta_t


# Erro & 3-sigma bound calculation
x_error = xtrue_sample-x_hat_all
position_3sigma = 3*np.sqrt(P_all[:,0])
velocity_3sigma = 3*np.sqrt(P_all[:,-1])

# Plot
tout_measure = np.linspace(0,(zk_all.size-1)*Delta_tk_KF,zk_all.size)

import matplotlib.pyplot as plt

fig_ms, (ax_ms_0, ax_ms_1) = plt.subplots(nrows=2,ncols=1)
ax_ms_0.plot(tout_measure,zk_all,'k.')
ax_ms_0.plot(tout_all, xout_all[:,0],'b--')
ax_ms_0.plot(tout_measure,x_hat_all[:,0],'r-')
ax_ms_0.set_ylabel('position [m]',fontsize=14)
ax_ms_0.set(xlim=(0, final_time),ylim=(-10,10))
ax_ms_0.legend(('Measurement','True','Estimated'),loc='upper left')

ax_ms_1.plot(tout_all, xout_all[:,1],'b--')
ax_ms_1.plot(tout_measure,x_hat_all[:,1],'r-')
ax_ms_1.set_ylabel('velocity [m/s]',fontsize=14)
ax_ms_1.set_xlabel('time [s]',fontsize=14)
ax_ms_1.set(xlim=(0, final_time),ylim=(-10,10))
ax_ms_1.legend(('True','Estimated'),loc='upper left')

fig_ms.savefig('../figures/msd_KF_state_python.pdf',dpi=600)

fig_ns, (ax_ns_0, ax_ns_1) = plt.subplots(nrows=2,ncols=1)
ax_ns_0.plot(tout_measure,x_error[:,0])
ax_ns_0.plot(tout_measure,position_3sigma,'r--')
ax_ns_0.plot(tout_measure,-position_3sigma,'r--')
ax_ns_0.set_ylabel('[m]')
ax_ns_0.set(xlim=(0, final_time),ylim=(-2,2))
ax_ns_0.legend(('Position Estimation Error',r'3$\sigma$ bounds'),loc='upper right')

ax_ns_1.plot(tout_measure,x_error[:,1])
ax_ns_1.plot(tout_measure,velocity_3sigma,'r--')
ax_ns_1.plot(tout_measure,-velocity_3sigma,'r--')
ax_ns_1.set_ylabel('[m/s]')
ax_ns_1.set_xlabel('time [s]',fontsize=14)
ax_ns_1.set(xlim=(0, final_time),ylim=(-10,10))
ax_ns_1.legend(('Velocity Estimation Error',r'3$\sigma$ bounds'),loc='upper right')

fig_ns.savefig('../figures/msd_KF_states_error_python.pdf',dpi=600)

