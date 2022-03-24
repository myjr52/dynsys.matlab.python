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

# number of cells
num_cell = 2

# simulation time values
time_current = 0    # initial time 
time_final   = 300.0 # final time [min]
time_record  = time_current # data record time
dt_record    = 0.1  # minimum time interval for data recording
max_num_data = np.floor((time_final-time_current)/dt_record+0.5);

# kinetic parameters for the Laub-Loomis Dicty cAMP oscillation 
# network model from k1 to k14
ki_para_org = np.array([2.0, 0.9, 2.5, 1.5, 0.6, 0.8, 1.0, 1.3, 0.3, 0.8, 0.7, 4.9, 23.0, 4.5])
Cell_Vol = 3.672e-14; # [litre]
NA = 6.022e23;        # Avogadro's number
num_molecule_species = 7
num_reactions = 14

# robustness
delta_worst = np.array([-1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, -1, 1])
p_delta = 2;
ki_para=ki_para_org*(1+(p_delta/100)*delta_worst)

# nominal initial number of molecules
ACA_n   = 35403  # [# of molecules]
PKA_n   = 32888  # [# of molecules]
ERK2_n  = 11838  # [# of molecules]
REGA_n  = 27348  # [# of molecules]
icAMP_n = 15489  # [# of molecules]
ecAMP_n = 4980   # [# of molecules]
CAR1_n  = 25423  # [# of molecules]

# initial conditions for each cell
ACA     = ACA_n*np.ones(num_cell) + np.random.randint(1,ACA_n+1,num_cell)
PKA     = PKA_n*np.ones(num_cell) + np.random.randint(1,PKA_n+1,num_cell)
ERK2    = ERK2_n*np.ones(num_cell) + np.random.randint(1,ERK2_n+1,num_cell)
REGA    = REGA_n*np.ones(num_cell) + np.random.randint(1,REGA_n+1,num_cell)
icAMP   = icAMP_n*np.ones(num_cell) + np.random.randint(1,icAMP_n+1,num_cell)
CAR1    = CAR1_n*np.ones(num_cell) + np.random.randint(1,CAR1_n+1,num_cell)

# total external cAMP
ecAMP_total = num_cell*ecAMP_n + np.random.randint(1,ecAMP_n+1)
# the total external cAMP distributed equally to each cell
ecAMP   = (ecAMP_total/num_cell)*np.ones(num_cell)

# storing data
species_all = np.zeros((int(max_num_data), 2))
species_all[0,:] = np.array([time_current, ecAMP_total])
data_idx = 0

propensity_a = np.zeros((num_reactions,num_cell))

while data_idx < max_num_data-1:
   
    # each cell
    for idx_cell in range(num_cell):
        propensity_a[0,idx_cell] = ki_para[0]*CAR1[idx_cell]
        propensity_a[1,idx_cell] = ki_para[1]*ACA[idx_cell]*PKA[idx_cell]/(NA*Cell_Vol*1e-6)
        propensity_a[2,idx_cell] = ki_para[2]*icAMP[idx_cell]
        propensity_a[3,idx_cell] = ki_para[3]*PKA[idx_cell]
        propensity_a[4,idx_cell] = ki_para[4]*CAR1[idx_cell]
        propensity_a[5,idx_cell] = ki_para[5]*PKA[idx_cell]*ERK2[idx_cell]/(NA*Cell_Vol*1e-6)
        propensity_a[6,idx_cell] = ki_para[6]*(NA*Cell_Vol*1e-6)
        propensity_a[7,idx_cell] = ki_para[7]*ERK2[idx_cell]*REGA[idx_cell]/(NA*Cell_Vol*1e-6)
        propensity_a[8,idx_cell] = ki_para[8]*ACA[idx_cell]
        propensity_a[9,idx_cell] = ki_para[9]*REGA[idx_cell]*icAMP[idx_cell]/(NA*Cell_Vol*1e-6)
        propensity_a[10,idx_cell] = ki_para[10]*ACA[idx_cell]
        propensity_a[11,idx_cell] = ki_para[11]*ecAMP[idx_cell]
        propensity_a[12,idx_cell] = ki_para[12]*ecAMP[idx_cell]
        propensity_a[13,idx_cell] = ki_para[13]*CAR1[idx_cell]
    
    # determine the reaction time tau
    propensity_vec = propensity_a.T.flatten()
    sum_propensity_a = np.sum(propensity_vec)
    tau = np.random.exponential(1/sum_propensity_a)

    # determine the reaction
    normalized_propensity_a = propensity_vec/sum_propensity_a
    cumsum_propensity_a = np.cumsum(normalized_propensity_a)
    which_reaction = np.random.rand(1)
    reaction_idx = np.cumsum((cumsum_propensity_a-which_reaction)<0)
    active_reaction = reaction_idx[-1]
    reaction = np.remainder(active_reaction,num_reactions)
    reaction_cell_num = np.int(np.floor(active_reaction/num_reactions))

    # update number of molecules
    if reaction==0:
        ACA[reaction_cell_num] = ACA[reaction_cell_num] + 1
    elif reaction==1:
        ACA[reaction_cell_num] = ACA[reaction_cell_num] - 1
    elif reaction==2:
        PKA[reaction_cell_num] = PKA[reaction_cell_num] + 1
    elif reaction==3:
        PKA[reaction_cell_num] = PKA[reaction_cell_num] - 1
    elif reaction==4:
        ERK2[reaction_cell_num] = ERK2[reaction_cell_num] + 1
    elif reaction==5:
        ERK2[reaction_cell_num] = ERK2[reaction_cell_num] - 1
    elif reaction==6:
        REGA[reaction_cell_num] = REGA[reaction_cell_num] + 1
    elif reaction==7:
        REGA[reaction_cell_num] = REGA[reaction_cell_num] - 1
    elif reaction==8:
        icAMP[reaction_cell_num] = icAMP[reaction_cell_num] + 1
    elif reaction==9:
        icAMP[reaction_cell_num] = icAMP[reaction_cell_num] - 1
    elif reaction==10:
        #ecAMP[reaction_cell_num] = ecAMP[reaction_cell_num] + 1
        ecAMP_total = ecAMP_total + 1
    elif reaction==11:
        #ecAMP[reaction_cell_num] = ecAMP[reaction_cell_num] - 1
        ecAMP_total = ecAMP_total - 1
    elif reaction==12:
        CAR1[reaction_cell_num] = CAR1[reaction_cell_num] + 1
    elif reaction==13:
        CAR1[reaction_cell_num] = CAR1[reaction_cell_num] - 1
    else:
        print(reaction,'Wrong reaction number!')
    
    # distribute the total ecAMP equally to all cells,
    # where allows non-integer numbers
    ecAMP = (ecAMP_total/num_cell)*np.ones(num_cell)
    
    time_current = time_current + tau
    
    if time_record < time_current:
       data_idx = data_idx + 1
       species_all[data_idx,:] = np.array([time_current, ecAMP_total])
       time_record = time_record + dt_record
       print(time_record)
    
import matplotlib.pyplot as plt
figure, axis = plt.subplots(1,1)
tout = species_all[:,0]
ecAMP_total = species_all[:,1]
axis.plot(tout,ecAMP_total,'b-')
axis.set_ylabel('[# of molecules]', fontsize=14)
axis.set_xlabel('time [min]', fontsize=14)
axis.legend(('Total external cAMP'),loc='upper right', fontsize=14)