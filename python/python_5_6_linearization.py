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
