#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 08:35:38 2021

@author: jongrae
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as spsg

kP = 20
kI = 2.5e-4
gamma_G = 8e-4
ks2 = 4e-4

A_PI = np.array([[0, 0, 0], [0, -gamma_G, 0], [ks2, ks2, -ks2]])
B_PI = np.array([[kI], [gamma_G*kP], [0]])
C_PI = np.array([[0, 0, 1]])
D_PI = np.array([[0]])

sys_PI = spsg.lti(A_PI,B_PI,C_PI,D_PI)

A_true = np.array([[0]])
B_true = np.array([[1]])
C_true = np.array([[kI]])
D_true = np.array([[kP]])

sys_true_PI = spsg.lti(A_true,B_true,C_true,D_true)


# bode plots
freq = np.logspace(-7,-3,1000) # [rad/time]
ww1, mm1, pp1 = spsg.bode(sys_PI,w=freq)
ww2, mm2, pp2 = spsg.bode(sys_true_PI,w=freq)

fig1, (ax1,ax2) = plt.subplots(nrows=2,ncols=1)
ax1.semilogx(ww1,mm1,'b-',ww2,mm2,'r--')
ax2.semilogx(ww1,pp1,'b-',ww2,pp2,'r--')

ax1.legend(('approximated PI','true PI'),fontsize=14)
ax2.legend(('approximated PI','true PI'),fontsize=14)

ax1.axis([1e-7,1e-3,0,65])
ax2.axis([1e-7,1e-3,-150,0.0])

ax1.set_ylabel('Magnitude [dB]',fontsize=14)
ax2.set_ylabel('Phase [$\circ$]',fontsize=14)

ax2.set_xlabel('Frequency [rad/time]',fontsize=14)

# step response and impulse response
time_sim = np.linspace(0,30000,30000)
ts1, ys1 = spsg.step(sys_PI,T=time_sim)
ts2, ys2 = spsg.step(sys_true_PI,T=time_sim)

tp1, yp1 = spsg.impulse(sys_PI,T=time_sim)
tp2, yp2 = spsg.impulse(sys_true_PI,T=time_sim)

fig2, (ax3,ax4) = plt.subplots(nrows=2,ncols=1)
ax3.plot(ts1/60,ys1,'b-',ts2/60,ys2,'r--')
ax4.plot(tp1/60,yp1,'b-',tp2/60,yp2,'r--')

ax3.legend(('approximated PI','true PI'),fontsize=14)
ax4.legend(('approximated PI','true PI'),fontsize=14)

ax3.axis([0,500,0,30])
ax4.axis([0,500,0,0.005])

ax3.set_ylabel('[a.u.]',fontsize=14)
ax4.set_ylabel('[a.u.]',fontsize=14)

ax4.set_xlabel('time [minutes]',fontsize=14)

xtick_list = np.array([0,100,200,300,400,500])
ax3.set_xticks(xtick_list)
ax3.set_xticklabels(xtick_list,fontsize=14)
ax4.set_xticks(xtick_list)
ax4.set_xticklabels(xtick_list,fontsize=14)

ytick_list1 = np.array([0,10,20,30])
ax3.set_yticks(ytick_list1)
ax3.set_yticklabels(ytick_list1,fontsize=14)

ytick_list2 = np.array([0,0.0025,0.005])
ax4.set_yticks(ytick_list2)
ax4.set_yticklabels(ytick_list2,fontsize=14)

ax3.set_title('Step Response',fontsize=14)
ax4.set_title('Impulse Response',fontsize=14)