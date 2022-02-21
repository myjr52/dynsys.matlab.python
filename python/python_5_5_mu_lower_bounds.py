#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 12:30:13 2021

@author: jongrae
"""

import numpy as np

Ns = 1000
eps = 1e-6

d_lb = 1e-3
d_ub = 10

d_ulb = d_ub - d_lb

omega = 0
Mjw = 1/(omega*1j + 2)

num_delta = 2

def A_delta(delta_1, delta_2):
    
    return -2+delta_1+np.sin(delta_1*delta_2)

while d_ulb > eps:
    
    d = (d_lb+d_ub)/2
    
    delta_1 = np.random.rand(Ns)*d-d/2
    delta_2 = np.random.rand(Ns)*d-d/2
    
    rand_face = np.random.randint(1,num_delta+1,Ns)
    delta_1[rand_face==0]=d/2
    delta_2[rand_face==1]=d/2
    
    Delta = A_delta(delta_1,delta_2)-A_delta(0,0)
    
    if np.unique(np.sign(np.real(1-Mjw*Delta))).size == 2:
        d_ub = d
    else:
        d_lb = d
        
    d_ulb = d_ub - d_lb
    

mu_lb = 2/d_ub
