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
plt.ylabel('position [m]')
plt.xlabel('time [s]')

plt.figure(2);
plt.plot(tout,xout[1,:])
plt.ylabel('velocity [m/s]')
plt.xlabel('time [s]')

plt.figure(3)
plt.plot(tout,xout[2,:])
plt.ylabel('m(t) [kg]')
plt.xlabel('time [s]')

from matplotlib.gridspec import GridSpec
fig = plt.figure(2)
gs = GridSpec(2,2,figure=fig)

ax1 = fig.add_subplot(gs[0,0])
ax1.plot(tout,xout[0,:])
ax1.set_ylabel('position [m]')
ax1.set_xlabel('time [s]')

ax2 = fig.add_subplot(gs[0,1])
ax2.plot(tout,xout[1,:])
ax2.set_ylabel('velocity [m/s]')
ax2.set_xlabel('time [s]')

ax3 = fig.add_subplot(gs[1,:])
ax3.plot(tout,xout[2,:])
ax3.set_ylabel('m(t) [kg]')
ax3.set_xlabel('time [s]')
