#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 01:12:05 2021

@author: jongrae
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as pltwk_

m_mass = 1.0 #[kg]
k_spring = 0.5 #[N/m]
c_damper = 0.1 #[N/(m/s)]
msd_const = [m_mass, k_spring, c_damper]

init_pos = 0.0 #[m]
init_vel = 0.0 #[m/s]

init_time = 0 #[s]
final_time = 300 #[s]
Delta_t = 0.01 #[s]
time_interval = [init_time, final_time]

num_w = int((final_time-init_time)/Delta_t)+1
sigma_beta = np.sqrt(0.5)
sigma_w =sigma_beta/np.sqrt(Delta_t)
wk_noise = sigma_w*(np.random.randn(num_w))

x0 = [init_pos, init_vel]
t0 = init_time
tf = t0 + Delta_t

tout_all = np.zeros(num_w)
xout_all = np.zeros((num_w,2))

tout_all[0] = t0
xout_all[0] = x0

x0_xk_const = x0;
tout_xk_const_all = np.zeros(num_w)
xout_xk_const_all = np.zeros((num_w,2))

tout_xk_const_all[0] = t0
xout_xk_const_all[0] = x0_xk_const

def msd_noisy(time,state,wk,msd_const):
    x1, x2 = state
    m, k, c = msd_const
    
    dxdt = [x2,
            -(k/m)*x1 - (c/m)*x2 + wk]
    return dxdt

for idx in range(1,num_w):
    
    wk = wk_noise[idx]

    # RK45
    sol = solve_ivp(msd_noisy, (t0, tf), x0, args=(wk,msd_const))
    xout = sol.y.transpose()    

    tout_all[idx] = sol.t[-1]
    xout_all[idx] = xout[-1]
    
    x0 = xout[-1];
    
    # time interval update
    t0 = tf
    tf = t0 + Delta_t
