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
    q4 = quat[3]
    
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

# assume random quaternion and the dcm are the satellite attitude
# with respect to the reference
q13 = eig_ax*np.sin(eig_ang/2)
q4 = np.cos(eig_ang/2)
quat_BR = np.hstack((q13,q4))
dcm_BR = q2dcm(quat_BR)

# dcm from body frame to sensor frame
dcm_SB = np.array([[1, 0, 0], [0, -1, 0],[0, 0, -1]])

# generate random stars in reference frame
total_num_star = 200       
rR_all = 2*np.random.rand(3,total_num_star)-1
rR_all = rR_all/np.sqrt(sum(rR_all**2)) # normalize it

rB_all = dcm_BR@rR_all
rS_all = dcm_SB@rB_all

# star sensor direction in the sensor frame
star_sensor_in_S = np.array([[0], [0], [1]])
star_sensor_fov = 12*np.pi/180 # 12 degrees in [radian]

# angles between the starts and the star sensor
th_all = np.arccos(rS_all.T@star_sensor_in_S)
star_seen = th_all < star_sensor_fov
rB_seen = rB_all[:,star_seen.squeeze()]
        

if rB_seen.shape[1] > 1:
    rR_seen = rR_all[:,star_seen.squeeze()]
    
    if rB_seen.shape[1] == 2:
        rB_3 = np.cross(rB_seen[:,0],rB_seen[:,1])
        rB_seen = np.hstack((rB_seen,rB_3.reshape((3,1))))
        
        
        rR_3 = np.cross(rR_seen[:,0],rR_seen[:,1])
        rR_seen = np.hstack((rR_seen,rR_3.reshape((3,1))))
    
    a_i = np.ones(rB_seen.shape[1])
    q_est, gama, sw = QUEST_Core(rB_seen,rR_seen,a_i);
    dcm_est = q2dcm(q_est);
    
    print('quaternion error = ',np.linalg.norm(quat_BR-q_est))
    print('dcm error = ',np.linalg.norm(dcm_BR-dcm_est))

else:
    print('No enough number of stars are seen')


#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
