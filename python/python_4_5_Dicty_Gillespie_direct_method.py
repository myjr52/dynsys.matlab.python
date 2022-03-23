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

# simulation time values
time_current = 0    # initial time 
time_final   = 30.0 # final time [min]
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
p_delta = 0;
ki_para=ki_para_org*(1+(p_delta/100)*delta_worst)

# initial number of molecules
ACA   = 35403  # [# of molecules]
PKA   = 32888  # [# of molecules]
ERK2  = 11838  # [# of molecules]
REGA  = 27348  # [# of molecules]
icAMP = 15489  # [# of molecules]
ecAMP = 4980   # [# of molecules]
CAR1  = 25423  # [# of molecules]

# storing data
species_all = np.zeros((int(max_num_data), num_molecule_species+1))
species_all[0,:] = np.array([time_current, ACA, PKA, ERK2, REGA, icAMP, ecAMP, CAR1])
data_idx = 0

propensity_a = np.zeros(num_reactions)

while data_idx < max_num_data-1:
   
    propensity_a[0] = ki_para[0]*CAR1
    propensity_a[1] = ki_para[1]*ACA*PKA/(NA*Cell_Vol*1e-6)
    propensity_a[2] = ki_para[2]*icAMP
    propensity_a[3] = ki_para[3]*PKA
    propensity_a[4] = ki_para[4]*CAR1
    propensity_a[5] = ki_para[5]*PKA*ERK2/(NA*Cell_Vol*1e-6)
    propensity_a[6] = ki_para[6]*(NA*Cell_Vol*1e-6)
    propensity_a[7] = ki_para[7]*ERK2*REGA/(NA*Cell_Vol*1e-6)
    propensity_a[8] = ki_para[8]*ACA
    propensity_a[9] = ki_para[9]*REGA*icAMP/(NA*Cell_Vol*1e-6)
    propensity_a[10] = ki_para[10]*ACA
    propensity_a[11] = ki_para[11]*ecAMP
    propensity_a[12] = ki_para[12]*ecAMP
    propensity_a[13] = ki_para[13]*CAR1
    
    # determine the reaction time tau
    sum_propensity_a = np.sum(propensity_a)
    tau = np.random.exponential(1/sum_propensity_a)

    # determine the reaction
    normalized_propensity_a = propensity_a/sum_propensity_a
    cumsum_propensity_a = np.cumsum(normalized_propensity_a)
    which_reaction = np.random.rand(1)
    reaction_idx = np.cumsum((cumsum_propensity_a-which_reaction)<0)
    reaction = reaction_idx[-1]

    # update number of molecules
    if reaction==0:
        ACA = ACA + 1
    elif reaction==1:
        ACA = ACA - 1
    elif reaction==2:
        PKA = PKA + 1
    elif reaction==3:
        PKA = PKA - 1
    elif reaction==4:
        ERK2 = ERK2 + 1
    elif reaction==5:
        ERK2 = ERK2 - 1
    elif reaction==6:
        REGA = REGA + 1
    elif reaction==7:
        REGA = REGA - 1
    elif reaction==8:
        icAMP = icAMP + 1
    elif reaction==9:
        icAMP = icAMP - 1
    elif reaction==10:
        ecAMP = ecAMP + 1
    elif reaction==11:
        ecAMP = ecAMP - 1
    elif reaction==12:
        CAR1 = CAR1 + 1
    elif reaction==13:
        CAR1 = CAR1 - 1
    else:
        print(reaction,'Wrong reaction number!')
    
    time_current = time_current + tau
    
    if time_record < time_current:
       data_idx = data_idx + 1
       species_all[data_idx,:] = np.array([time_current, ACA, PKA, ERK2, REGA, icAMP, ecAMP, CAR1])
       time_record = time_record + dt_record
       print(time_record)
    
import matplotlib.pyplot as plt
figure, axis = plt.subplots(1,1)
tout = species_all[:,0]
icAMP = species_all[:,5]
CAR1 = species_all[:,7]
axis.plot(tout,icAMP,'b-')
axis.plot(tout,CAR1,'r-.')
axis.set_ylabel('[# of molecules]', fontsize=14)
axis.set_xlabel('time [min]', fontsize=14)
axis.legend(('i-cAMP','CAR 1'),loc='upper right', fontsize=14)

axis.axis([0,30,0,60000])
xtick = np.array([0,5,10,15,20,25,30])
ytick = np.array([0, 10000, 20000, 30000, 40000, 50000, 60000])
axis.set_xticks(xtick)
axis.set_yticks(ytick)
axis.set_xticklabels(xtick,fontsize=14)
axis.set_yticklabels(ytick,fontsize=14)
#figure.savefig('python_Dicty_Gillespie.pdf',bbox_inches='tight')
