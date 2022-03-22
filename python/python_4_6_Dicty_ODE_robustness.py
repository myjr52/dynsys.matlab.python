#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 23:11:34 2020

@author: menjkim
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def Dicty_cAMP(time,state,ki_para):
    ACA, PKA, ERK2, REGA, icAMP, ecAMP, CAR1  = state

    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14 = ki_para
    
    dACA_dt   = k1*CAR1 - k2*ACA*PKA
    dPKA_dt   = k3*icAMP - k4*PKA
    dERK2_dt  = k5*CAR1 - k6*PKA*ERK2
    dREGA_dt  = k7 - k8*ERK2*REGA
    dicAMP_dt = k9*ACA - k10*REGA*icAMP
    decAMP_dt = k11*ACA - k12*ecAMP
    dCAR1_dt  = k13*ecAMP - k14*CAR1
    
    dxdt = [dACA_dt,
            dPKA_dt,
            dERK2_dt,
            dREGA_dt,
            dicAMP_dt,
            decAMP_dt,
            dCAR1_dt]
    return dxdt

# Cost function to be minimized for robustness analysis
def Dicty_x1_square_integral(delta):

    ki_para_org = np.array([2.0, 0.9, 2.5, 1.5, 0.6, 0.8, 1.0, 1.3, 0.3, 0.8, 0.7, 4.9, 23.0, 4.5])
    p_delta = 2 # [percents]
    ki_para=ki_para_org*(1+(p_delta/100)*delta)

    init_cond = np.random.rand(7)
    dt = 0.1  # [minutes]
    t0 = 600  # [minutes]
    tf = 1200 # [minutes]
    time_interval = np.linspace(0,tf,int(tf/dt)) # [min]
 
    sol_out = solve_ivp(Dicty_cAMP, (0, tf), init_cond, t_eval=time_interval, args=(ki_para,))
    xout = sol_out.y

    N_t0 = int(t0/dt) - 1

    ACA = xout[0,N_t0::]
    PKA = xout[2,N_t0::]
    CAR1 = xout[6,N_t0::]
    J_cost = np.sum((ki_para[0]*CAR1 - ki_para[1]*(ACA*PKA))**2)*dt*0.5
    
    plt.plot(ACA)
    
    return J_cost

#-------------------------------------------
# Find worst perturbation
#-------------------------------------------
delta_dim = 14;
delta = (2*np.random.rand(delta_dim)-1)

bounds = [(-1,1)]*delta_dim

from scipy import optimize
result = optimize.differential_evolution(Dicty_x1_square_integral, bounds,
                                          updating='deferred', disp=True, popsize=200, maxiter=100, workers=4)