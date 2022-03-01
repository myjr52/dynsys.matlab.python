#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: jongrae
"""

import numpy as np

#-------------------------------------------------------------------------
# BEGIN of QUEST
#-------------------------------------------------------------------------
def QUEST_Core(rB,rR,a_w):
    
    w = rB
    v = rR
    
    max_iter_sw = False
    
    num_obs = v.shape[1]
    
    B = np.zeros((3,3))
    z = np.zeros((3,1))
    
    for idx in range(num_obs):
        B = B + a_w[idx]*np.reshape(w[::,idx],(3,1))@np.reshape(v[::,idx],(1,3))
        z = z + a_w[idx]*np.cross(w[::,idx],v[::,idx]).reshape((3,1))
    
    S = B + B.transpose()
    sgm = B[0,0]+B[1,1]+B[2,2]
    
    kappa = (S[1,1]*S[2,2]-S[1,2]**2) + (S[0,0]*S[2,2]-S[0,2]**2) + (S[0,0]*S[1,1]-S[0,1]**2)
    
    delta = np.linalg.det(S)
    
    a = sgm**2 - kappa
    b = (sgm**2 + z.transpose()@z)[0,0]
    c = (delta  +z.transpose()@S@z)[0,0]
    d = (z.transpose()@S@S@z)[0,0]
    
    #--------------------------------------
    # Newton-Raphson Method to find a root
    #--------------------------------------
    lamda_c = 10 # initial guess greater than 1
    tol_NR = True
    max_num_itr = 50
    cur_itr = 0
    
    while tol_NR:
        lamda_p =  lamda_c
        
        f_c = lamda_c**4 - (a+b)*lamda_c**2 - c*lamda_c + (a*b+c*sgm-d)
        
        dfdx = 4*lamda_c**3 - 2*(a+b)*lamda_c - c
        lamda_c = lamda_c - f_c/dfdx
        if abs(lamda_c - lamda_p) < 1e-6:
            tol_NR = False
    
        cur_itr = cur_itr + 1   
        
        if cur_itr > max_num_itr:
            max_iter_sw = True
            break
    
    lamda = lamda_c
    gama = np.abs(np.linalg.det((sgm+lamda)*np.eye(3)-S))

    Y = np.linalg.inv((sgm+lamda)*np.eye(3)-S)@z

    q_est = 1/np.sqrt(1+np.sum(Y*Y))*np.append(Y,1)
    q_est = q_est/np.linalg.norm(q_est)
    
    return q_est, gama, max_iter_sw
#-------------------------------------------------------------------------
# END of QUEST
#-------------------------------------------------------------------------

def q2dcm(quat): 

    q13 =quat[0:3]
    q4 = quat[3];
    
    q13x = np.array([   [ 0,          -q13[2],       q13[1]],
                        [ q13[2],      0,           -q13[0]],
                        [-q13[1],      q13[0],        0]])
    
    
    q13 = np.reshape(q13,(3,1))
    dcm = (q4**2-np.squeeze(q13.T@q13))*np.eye(3) + 2*q13@q13.T - 2*q4*q13x
    
    return dcm


#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
# generate random quaternion
eig_ax = 2*(np.random.rand(3)-0.5)
eig_ax = eig_ax/np.linalg.norm(eig_ax) #random axis normalized

eig_ang = np.pi*2*(np.random.rand(1)-0.5) # angle between +/- pi [rad]

# random quaternion
q13 = eig_ax*np.sin(eig_ang/2)
q4 = np.cos(eig_ang/2)
quat = np.hstack((q13,q4))

# random quaternion to dcm and the dcm to quaternion
dcm = q2dcm(quat)

# 5 stars in the reference frame
r1R = np.array([-0.6794, -0.3237, -0.6586]).reshape((3,1)); r1R = r1R/np.linalg.norm(r1R)
r2R = np.array([-0.7296,  0.5858,  0.3528]).reshape((3,1)); r2R = r2R/np.linalg.norm(r2R)
r3R = np.array([-0.2718,  0.6690, -0.6918]).reshape((3,1)); r3R = r3R/np.linalg.norm(r3R)
r4R = np.array([-0.2062, -0.3986,  0.8936]).reshape((3,1)); r4R = r4R/np.linalg.norm(r4R)
r5R = np.array([0.6858,  -0.7274, -0.0238]).reshape((3,1)); r5R = r5R/np.linalg.norm(r5R)

# 5 starts in the body frame
r1B = dcm@r1R
r2B = dcm@r2R
r3B = dcm@r3R
r4B = dcm@r4R
r5B = dcm@r5R
        
rB_all = np.hstack((r1B,r2B,r3B,r4B,r5B))
rR_all = np.hstack((r1R,r2R,r3R,r4R,r5R))
a_i = np.ones(rB_all.shape[1])

q_est, gama, sw = QUEST_Core(rB_all,rR_all,a_i)
dcm_est = q2dcm(q_est)

print('quaternion error =', np.linalg.norm(quat-q_est))
print('dcm error = ', np.linalg.norm(dcm-dcm_est))
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------