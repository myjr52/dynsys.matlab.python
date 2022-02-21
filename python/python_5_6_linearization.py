#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 15:38:44 2021

@author: jongrae
"""

from sympy import symbols, Matrix

kon, kcat, kdg, Kdeg, kI, eta, gamma_G, kP, ks2, ks1, Pd, X3 = symbols('kon kcat kdg Kdeg kI eta gamma_G kP ks2 ks1 Pd X3')
E, P, ES, X1, X2, S, DP, Xs = symbols('E P ES X1 X2 S DP Xs')


dE_dt = -kon*E*S + kcat*ES;
dP_dt = kcat*ES - kdg*P - Kdeg*X3*P;
dES_dt = kon*E*S - kcat*ES;
dX1_dt = kI*DP - eta*X1;
dX2_dt = -gamma_G*X2 + gamma_G*kP*DP;
dS_dt = ks2*X1 + ks2*X2 - ks2*S;
dDP_dt = ks1*Pd - ks1*DP*Xs - ks1*DP;
dXs_dt = -ks1*DP*Xs + ks1*P;

fx = Matrix([[dE_dt], [dP_dt], [dES_dt], [dX1_dt], [dX2_dt], [dS_dt], [dDP_dt], [dXs_dt]])
state = Matrix([[E], [P], [ES], [X1], [X2], [S], [DP], [Xs]])

dfdx = fx.jacobian(state)

# Steady-state
Ess = 0.1998   
Pss = 0.4999    
ESss = 0.0002    
X1ss = 0.0558   
X2ss = 50.0415   
Sss = 50.0723   
DPss = 1.0001    
Xsss = 0.4999

dfdx_at_ss = dfdx.subs([[E,Ess],[P,Pss],[ES,ESss],[X1,X1ss],[X2,X2ss],[S,Sss],[DP,DPss],[Xs,Xsss]])

# nomial stability with the nominal parameters
kP = 50;
kI = 5e-6;
gamma_G = 8e-4;
ks2 = 4e-4;
Kdeg = 1e-3;
X3 = 1; 
KF = 3; 
eta = 1e-4;            
Pd = 1;
kon = 5e-5;
kcat = 1.6*2;
kdg = 8e-8;
ks1 = 3;


dfdx_nominal = dfdx_at_ss.subs([
     [symbols('kP'),kP], 
     [symbols('kI'),kI],
     [symbols('gamma_G'),gamma_G],
     [symbols('ks2'),ks2],
     [symbols('Kdeg'),Kdeg],
     [symbols('X3'),X3],
     [symbols('KF'),KF],
     [symbols('eta'),eta],
     [symbols('Pd'),Pd],
     [symbols('kon'),kon],
     [symbols('kcat'),kcat],
     [symbols('kdg'),kdg],
     [symbols('ks1'),ks1]
    ])

import numpy as np
dfdx_nominal_val = np.array(dfdx_nominal,dtype=np.float64)
[eig_val,eig_vec]=np.linalg.eig(dfdx_nominal_val)