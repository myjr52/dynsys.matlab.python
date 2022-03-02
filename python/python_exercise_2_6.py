#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 13:16:12 2022

@author: jongrae
"""

import numpy as np


def q2dcm(quat): 

    q13 =quat[0:3]
    q4 = quat[3];
    
    q13x = np.array([   [ 0,          -q13[2],       q13[1]],
                        [ q13[2],      0,           -q13[0]],
                        [-q13[1],      q13[0],        0]])
    
    
    q13 = np.reshape(q13,(3,1))
    dcm = (q4**2-np.squeeze(q13.T@q13))*np.eye(3) + 2*q13@q13.T - 2*q4*q13x
    
    return dcm


def dcm2q(dcm): 
    
    quat = np.zeros(4)

    a1 = (1 + dcm[0,0] - dcm[1,1] - dcm[2,2])/4
    a2 = (1 + dcm[1,1] - dcm[0,0] - dcm[2,2])/4
    a3 = (1 + dcm[2,2] - dcm[0,0] - dcm[1,1])/4
    a4 = (1 + dcm[0,0] + dcm[1,1] + dcm[2,2])/4

    a = np.array([a1, a2, a3, a4])
    a_idx = a.argmax()
    quat[a_idx] = np.sqrt(a.max())

    if a_idx==0:
        quat[1] = (dcm[0,1]+dcm[1,0])/(4*quat[0]);
        quat[2] = (dcm[0,2]+dcm[2,0])/(4*quat[0]);
        quat[3] = (dcm[1,2]-dcm[2,1])/(4*quat[0]);
    elif a_idx==1:
        quat[0] = (dcm[0,1]+dcm[1,0])/(4*quat[1]);
        quat[2] = (dcm[1,2]+dcm[2,1])/(4*quat[1]);
        quat[3] = (dcm[2,0]-dcm[0,2])/(4*quat[1]);
    elif a_idx==2:
        quat[0] = (dcm[0,2]+dcm[2,0])/(4*quat[2]);
        quat[1] = (dcm[1,2]+dcm[2,1])/(4*quat[2]);
        quat[3] = (dcm[0,1]-dcm[1,0])/(4*quat[2]);
    elif a_idx==3:
        quat[0] = (dcm[1,2]-dcm[2,1])/(4*quat[3]);
        quat[1] = (dcm[2,0]-dcm[0,2])/(4*quat[3]);
        quat[2] = (dcm[0,1]-dcm[1,0])/(4*quat[3]);

    if quat[3] < 0:
        quat = - quat

    return quat


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
quat_from_dcm = dcm2q(dcm)

# compare quaternion
print(np.linalg.norm(quat - quat_from_dcm))

# dcm from the converted quaternion
dcm2 = q2dcm(quat_from_dcm);

# compare dcm
print(np.linalg.norm(dcm-dcm2))