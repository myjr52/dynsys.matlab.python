#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 15:00:20 2021

@author: jongrae
"""

import numpy as np
from scipy.integrate import solve_ivp




def robustness_analysis_MC(sim_para):

    J_inertia = 1.0
    init_time = sim_para[0]
    final_time = sim_para[1]
    state_0 = sim_para[2]
    tout = sim_para[3]
    
    for g_MC_idx in range(num_MC):
        
        dJ = np.random.randn(1)

        J_inv_perturbed = J_inertia = 
        
        quadcopter_uav=(J_inertia_perturbed,)
        
        sol = solve_ivp(dqdt_dwdt_drvdt, (init_time, final_time), state_0, t_eval=tout, 
                        atol=1e-9, rtol=1e-6, max_step=0.01, 
                        args=(quadcopter_uav,))
        qout = sol.y[0:4,:]

        
        print(f'#{g_MC_idx:1d} {np.linalg.norm(dJ):6.5f}, {q13_ts:4.2f}\n')

        return ts_all#, dJ_norm_all


#-----------------------------------------------------------------------------
# Run Parallel Monte-Carlo Simulation
#-----------------------------------------------------------------------------
from multiprocessing import Pool
with Pool(5) as p:
    print(p.map(robustness_analysis_MC, [1,2]))