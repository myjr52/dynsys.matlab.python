#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 14:50:47 2021

@author: jongrae
"""

import numpy as np
from sympy import symbols, simplify, expand

Dt, ux0, uy0, ux1, uy1, wx0, wy0, wx1, wy1 = symbols('Dt ux0 uy0 ux1 uy1 wx0 wy0 wx1 wy1 ')
xa0, ya0, vxa0, vya0, xt0, yt0, th, w_max = symbols('xa0 ya0 vxa0 vya0 xt0 yt0 th w_max')

# Dynamics
Fa = np.eye(4)+np.vstack((np.hstack((np.zeros((2,2)),Dt*np.eye(2))),np.zeros((2,4)))) 
Ga = np.vstack((np.zeros((2,2)),Dt*np.eye(2)))
Ca = np.eye(2,4)

Ft = np.eye(2)
Gt = Dt*np.eye(2)
Ct = np.eye(2)

# control inputs
u_vec_0 = np.array([[ux0], [uy0]])
w_vec_0 = np.array([[wx0], [wy0]])
u_vec_1 = np.array([[ux1], [uy1]])
w_vec_1 = np.array([[wx1], [wy1]])

# initial conditions
xa_vec_0 = np.array([[xa0], [ya0], [vxa0], [vya0]])
xt_vec_0 = np.array([[xt0], [yt0]])

# state propagation
xa_k_plus_1 = Fa@xa_vec_0    + Ga@u_vec_0
xa_k_plus_2 = Fa@xa_k_plus_1 + Ga@u_vec_1;
y_k_plus_1 = Ca@xa_k_plus_1;
y_k_plus_2 = Ca@xa_k_plus_2;

xt_k_plus_1 = Ft@xt_vec_0    + Gt@w_vec_0
xt_k_plus_2 = Ft@xt_k_plus_1 + Gt@w_vec_1
z_k_plus_1 = Ct@xt_k_plus_1
z_k_plus_2 = Ct@xt_k_plus_2

#----------------------------------------------------------------------
# calculate the cost function in the original form
#----------------------------------------------------------------------
dyz_1 = y_k_plus_1-z_k_plus_1
dyz_2 = y_k_plus_2-z_k_plus_2
J_over_dt_2 = dyz_1.T@dyz_1+dyz_2.T@dyz_2
J_over_dt_2 = simplify(expand(J_over_dt_2[0][0]))

alpha = J_over_dt_2.coeff(ux0,1)/(Dt**4)

temp = J_over_dt_2.coeff(ux0,0)
beta = temp.coeff(uy0,1)/(Dt**4)
gama = temp.coeff(uy0,0)

poly_recover = alpha*(Dt**4)*ux0 + (Dt**4)*(ux0**2) + beta*(Dt**4)*uy0 + (Dt**4)*(uy0**2) + gama

# check alpha, beta, gama are correct: the following must return zero
zero_check = float(expand(poly_recover-J_over_dt_2))
print(f"Is this zero? {zero_check:4.2f}\n")

