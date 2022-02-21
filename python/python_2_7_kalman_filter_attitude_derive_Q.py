 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 21:54:22 2021

@author: jongrae
"""

from sympy import symbols, Matrix, simplify, expand, integrate

w1, w2, w3, Dt, nv1, nv2, nv3, nu1, nu2, nu3 = symbols('w1 w2 w3 Dt nv1 nv2 nv3 nu1 nu2 nu3')
sgm2_u, sgm2_v = symbols('sgm2_u sgm2_v') # these are variance, i.e., sigma-squared

wx = Matrix([[ 0, -w3, w2], [w3, 0, -w1], [-w2, w1, 0]])

nv = Matrix([[nv1],[nv2],[nv3]])
nu = Matrix([[nu1],[nu2],[nu3]])
wc = Matrix([nv,nu])

# F & G Matrices
F = Matrix([[-wx,-Matrix.eye(3)],[Matrix.zeros(3,6)]])
G = Matrix([[-Matrix.eye(3), Matrix.zeros(3)],[Matrix.zeros(3), Matrix.eye(3)]])

# e^{Ft}
Phi = Matrix.eye(6) + F*Dt + (1/2)*(F**2)*(Dt**2) + (1/6)*(F**3)*(Dt**3) + (1/24)*(F**4)*(Dt**4)

# wd before integral
wd = Phi@wc

# E(wd wd^T)
wd_wd_T = wd@wd.transpose()
Q_cov = Matrix.zeros(6)

# Q_11: integrate from 0 to Dt
cov_wd_11 = simplify(expand(wd_wd_T[0,0]))
cov_wd_11 = cov_wd_11.subs([[nu1**2,sgm2_u],[nu2**2,sgm2_u],[nu3**2,sgm2_u],
                            [nv1**2,sgm2_v],[nv2**2,sgm2_v],[nv3**2,sgm2_v]])
cov_wd_11 = cov_wd_11.subs([[nu1,0],[nu2,0],[nu3**2,0],[nv1,0],[nv2,0],[nv3,0]])

cov_wd_11 = integrate(cov_wd_11,(Dt,0,Dt))
cov_wd_11 = simplify(expand(cov_wd_11))
cov_wd_11 = cov_wd_11.subs([[Dt**4,0],[Dt**5,0],[Dt**6,0],[Dt**7,0],[Dt**8,0],[Dt**9,0]])
cov_wd_11 = expand(cov_wd_11)
Q_cov[0,0] = cov_wd_11

# Q_12 & Q_21: integrate from 0 to Dt
cov_wd_12 = simplify(expand(wd_wd_T[0,1]))
cov_wd_12 = cov_wd_12.subs([[nu1**2,sgm2_u],[nu2**2,sgm2_u],[nu3**2,sgm2_u],
                            [nv1**2,sgm2_v],[nv2**2,sgm2_v],[nv3**2,sgm2_v]])
cov_wd_12 = cov_wd_12.subs([[nu1,0],[nu2,0],[nu3**2,0],[nv1,0],[nv2,0],[nv3,0]])

cov_wd_12 = integrate(cov_wd_12,(Dt,0,Dt))
cov_wd_12 = simplify(expand(cov_wd_12))
cov_wd_12 = cov_wd_12.subs([[Dt**4,0],[Dt**5,0],[Dt**6,0],[Dt**7,0],[Dt**8,0],[Dt**9,0]])
cov_wd_12 = expand(cov_wd_12)
Q_cov[0,1] = cov_wd_12
Q_cov[1,0] = cov_wd_12

# Q_13 & Q_31: integrate from 0 to Dt
cov_wd_13 = simplify(expand(wd_wd_T[0,2]))
cov_wd_13 = cov_wd_13.subs([[nu1**2,sgm2_u],[nu2**2,sgm2_u],[nu3**2,sgm2_u],
                            [nv1**2,sgm2_v],[nv2**2,sgm2_v],[nv3**2,sgm2_v]])
cov_wd_13 = cov_wd_13.subs([[nu1,0],[nu2,0],[nu3**2,0],[nv1,0],[nv2,0],[nv3,0]])

cov_wd_13 = integrate(cov_wd_13,(Dt,0,Dt))
cov_wd_13 = simplify(expand(cov_wd_13))
cov_wd_13 = cov_wd_13.subs([[Dt**4,0],[Dt**5,0],[Dt**6,0],[Dt**7,0],[Dt**8,0],[Dt**9,0]])
cov_wd_13 = expand(cov_wd_13)
Q_cov[0,2] = cov_wd_13
Q_cov[2,0] = cov_wd_13

# Q_22: integrate from 0 to Dt
cov_wd_22 = simplify(expand(wd_wd_T[1,1]))
cov_wd_22 = cov_wd_22.subs([[nu1**2,sgm2_u],[nu2**2,sgm2_u],[nu3**2,sgm2_u],
                            [nv1**2,sgm2_v],[nv2**2,sgm2_v],[nv3**2,sgm2_v]])
cov_wd_22 = cov_wd_22.subs([[nu1,0],[nu2,0],[nu3**2,0],[nv1,0],[nv2,0],[nv3,0]])

cov_wd_22 = integrate(cov_wd_22,(Dt,0,Dt))
cov_wd_22 = simplify(expand(cov_wd_22))
cov_wd_22 = cov_wd_22.subs([[Dt**4,0],[Dt**5,0],[Dt**6,0],[Dt**7,0],[Dt**8,0],[Dt**9,0]])
cov_wd_22 = expand(cov_wd_22)
Q_cov[1,1] = cov_wd_22
Q_cov[1,1] = cov_wd_22

# Q_23: integrate from 0 to Dt
cov_wd_23 = simplify(expand(wd_wd_T[1,2]))
cov_wd_23 = cov_wd_22.subs([[nu1**2,sgm2_u],[nu2**2,sgm2_u],[nu3**2,sgm2_u],
                            [nv1**2,sgm2_v],[nv2**2,sgm2_v],[nv3**2,sgm2_v]])
cov_wd_23 = cov_wd_22.subs([[nu1,0],[nu2,0],[nu3**2,0],[nv1,0],[nv2,0],[nv3,0]])

cov_wd_23 = integrate(cov_wd_23,(Dt,0,Dt))
cov_wd_23 = simplify(expand(cov_wd_23))
cov_wd_23 = cov_wd_23.subs([[Dt**4,0],[Dt**5,0],[Dt**6,0],[Dt**7,0],[Dt**8,0],[Dt**9,0]])
cov_wd_23 = expand(cov_wd_23)
Q_cov[1,2] = cov_wd_23
Q_cov[1,2] = cov_wd_23

# Q_33: integrate from 0 to Dt
cov_wd_33 = simplify(expand(wd_wd_T[2,2]))
cov_wd_33 = cov_wd_33.subs([[nu1**2,sgm2_u],[nu2**2,sgm2_u],[nu3**2,sgm2_u],
                            [nv1**2,sgm2_v],[nv2**2,sgm2_v],[nv3**2,sgm2_v]])
cov_wd_33 = cov_wd_33.subs([[nu1,0],[nu2,0],[nu3**2,0],[nv1,0],[nv2,0],[nv3,0]])

cov_wd_33 = integrate(cov_wd_33,(Dt,0,Dt))
cov_wd_33 = simplify(expand(cov_wd_33))
cov_wd_33 = cov_wd_33.subs([[Dt**4,0],[Dt**5,0],[Dt**6,0],[Dt**7,0],[Dt**8,0],[Dt**9,0]])
cov_wd_33 = expand(cov_wd_33)
Q_cov[2,2] = cov_wd_33
Q_cov[2,2] = cov_wd_33

