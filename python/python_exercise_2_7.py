#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 16:42:40 2022

@author: jongrae
"""
import numpy as np
from scipy.sparse import block_diag
from scipy.linalg import solve

r1R = np.array([-0.6794, -0.3237, -0.6586]).reshape((3,1))
r2R = np.array([-0.7296,  0.5858,  0.3528]).reshape((3,1))
r3R = np.array([-0.2718,  0.6690, -0.6918]).reshape((3,1))
r4R = np.array([-0.2062, -0.3986, 0.8936]).reshape((3,1))
r5R = np.array([0.6858, -0.7274, -0.0238]).reshape((3,1))


A1 = block_diag((r1R.transpose(),r1R.transpose(),r1R.transpose()))
A2 = block_diag((r2R.transpose(),r2R.transpose(),r2R.transpose()))
A3 = block_diag((r3R.transpose(),r3R.transpose(),r3R.transpose()))
A4 = block_diag((r4R.transpose(),r4R.transpose(),r4R.transpose()))
A5 = block_diag((r5R.transpose(),r5R.transpose(),r5R.transpose()))

A = np.vstack((A1.todense(),A2.todense(),A3.todense(),A4.todense(),A5.todense()))

r1B = np.array([-0.2147, -0.7985, 0.562]).reshape((3,1))
r2B = np.array([-0.7658, 0.4424, 0.4667]).reshape((3,1))
r3B = np.array([-0.8575, -0.4610, -0.228]).reshape((3,1))
r4B = np.array([0.4442, 0.6863, 0.5758]).reshape((3,1))
r5B = np.array([0.9407, -0.1845, -0.2847]).reshape((3,1))

y= np.vstack((r1B,r2B,r3B,r4B,r5B))

vec_C = solve(A.T@A,A.T@y)

C_BR = vec_C.reshape(3,3)