#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 23:11:34 2020

@author: menjkim
"""

from numpy import linspace
from scipy.integrate import solve_ivp

grv_const = 9.81 # [m/s^2]
init_pos = 0.0 # [m]
init_vel = 0.5 # [m/s]
init_mass = 5.0 #[kg]

init_cond = [init_pos, init_vel, init_mass]

init_time = 0 # [s]
final_time = 5.0 # [s]
num_data = 100
tout = linspace(init_time, final_time, num_data)


def free_falling_obj(time, state, grv_const):
     x1, x2, x3 = state
     dxdt = [x2,
          grv_const + (x3-2)*(x2/x3),
          -x3 + 2]
     return dxdt


sol = solve_ivp(free_falling_obj, (init_time, final_time), init_cond, t_eval=tout, args=(grv_const,))
xout = sol.y

import matplotlib.pyplot as plt
plt.figure(1)
plt.plot(tout,xout[0,:])
plt.ylabel('position [m]');
plt.xlabel('time [s]');

plt.figure(2);
plt.plot(tout,xout[1,:])
plt.ylabel('velocity [m/s]');
plt.xlabel('time [s]');

plt.figure(3);
plt.plot(tout,xout[2,:])
plt.ylabel('m(t) [kg]');
plt.xlabel('time [s]');