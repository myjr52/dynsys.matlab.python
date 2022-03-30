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
from sympy import symbols, Matrix

kon, kcat, kdg, Kdeg, kI, eta, gamma_G, kP, ks2, ks1, Pd, X3 = symbols('kon kcat kdg Kdeg kI eta gamma_G kP ks2 ks1 Pd X3')
E, P, ES, X1, X2, S, DP, Xs = symbols('E P ES X1 X2 S DP Xs')

d1, d2, d3, d4, d5, d6, d7, d8 = symbols('d1 d2 d3 d4 d5 d6 d7 d8')
d9, d10, d11, d12, d13, d14, d15, d16 = symbols('d9 d10 d11 d12 d13 d14 d15 d16')

dE_dt = -(kon+d1)*E*S + (kcat+d2)*ES;
dP_dt = (kcat+d2)*ES - (kdg+d3)*P - (Kdeg+d4)*(X3+d5)*P;
dES_dt = (kon+d1)*E*S - (kcat+d2)*ES;
dX1_dt = (kI+d6)*DP - (eta+d7)*X1;
dX2_dt = -(gamma_G+d8)*X2 + (gamma_G+d9)*(kP+d10)*DP;
dS_dt = (ks2+d11)*X1 + (ks2+d12)*X2 - (ks2+d13)*S;
dDP_dt = (ks1+d14)*Pd - (ks1+d15)*DP*Xs - (ks1+d16)*DP;
dXs_dt = -(ks1+d15)*DP*Xs + (ks1+d16)*P;

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

# mu-analysis
Ns = 1000
eps = 1e-6

num_state = 8
num_delta = 16
A0 = dfdx_nominal_val = dfdx_nominal.subs([
    [symbols('d1'),0],
    [symbols('d2'),0],
    [symbols('d3'),0],
    [symbols('d4'),0],
    [symbols('d5'),0],
    [symbols('d6'),0],
    [symbols('d7'),0],
    [symbols('d8'),0],
    [symbols('d9'),0],
    [symbols('d10'),0],
    [symbols('d11'),0],
    [symbols('d12'),0],
    [symbols('d13'),0],
    [symbols('d14'),0],
    [symbols('d15'),0],
    [symbols('d16'),0]
    ])
A0 = np.array(A0,dtype=np.float64)

num_omega = 1
omega_all = np.hstack((0,np.logspace(-3,-1,num_omega-1)))

mu_ub = np.zeros(num_omega)

# lower bound using geometric approach
for wdx, omega in enumerate(omega_all):
    Mjw=np.linalg.inv(1j*omega*np.eye(num_state)-A0)
    
    d_lb = 1e-6
    delta_d_lb = 5e-6
    d = d_lb
    not_find_ub = True
    
    num_iteration = 1
    
    if omega==0:
        size_check = 2
    else:
        size_check = 4
    
    while not_find_ub:
        
        for idx in range(Ns):
            delta_vec = np.random.rand(num_delta)*d-d/2
            
            dfdx_nominal_val = dfdx_nominal.subs([
                [symbols('d1'),delta_vec[0]],
                [symbols('d2'),delta_vec[1]],
                [symbols('d3'),delta_vec[2]],
                [symbols('d4'),delta_vec[3]],
                [symbols('d5'),delta_vec[4]],
                [symbols('d6'),delta_vec[5]],
                [symbols('d7'),delta_vec[6]],
                [symbols('d8'),delta_vec[7]],
                [symbols('d9'),delta_vec[8]],
                [symbols('d10'),delta_vec[9]],
                [symbols('d11'),delta_vec[10]],
                [symbols('d12'),delta_vec[11]],
                [symbols('d13'),delta_vec[12]],
                [symbols('d14'),delta_vec[13]],
                [symbols('d15'),delta_vec[14]],
                [symbols('d16'),delta_vec[15]]
                ])

            Delta = np.array(dfdx_nominal_val,dtype=np.float64) - A0
            
            I_MD = np.linalg.det(np.eye(num_state)-Mjw*Delta)
            fR = np.sign(np.real(I_MD))
            fI = np.sign(np.imag(I_MD))
            
            if idx==0:
                sign_all = np.array([fR, fI])
            else:
                sign_all = np.vstack((sign_all,[fR, fI]))
                sign_all = np.unique(sign_all,axis=0)

        if sign_all.shape[0] == size_check:
            d_lb = d_lb + delta_d_lb
        else:
            not_find_ub = False
                    
        print(d_lb)
        print(sign_all)
        
        d = d_lb
        num_iteration = num_iteration*2
        
            
    mu_ub[wdx] = 2/d