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
import matplotlib.pyplot as plt
from matplotlib import path

from numpy import sqrt
import time

###############################################
## input variables
###############################################
# initial uav position
xa0 = (2*np.random.rand(1)-1)*100 #[m] 
ya0 = (2*np.random.rand(1)-1)*100 #[m]
xa0 = xa0[0]
ya0 = ya0[0]

# initial uav velocity
tha0 = np.random.rand(1)*2*np.pi #[radian]
tha0 = tha0[0]
current_speed = 25 #[m/s]
vxa0 = current_speed*np.cos(tha0)
vya0 = current_speed*np.sin(tha0)

# control acceleration input magnitude constraints
ux_max = 10 # [m/s^2]
ux_min = -1 # [m/s^2]
uy_max = 2  # [m/s^2]
uy_min = -2 # [m/s^2]

# uav minimum & maximum speed
v_min = 20  #[m/s]
v_max = 40  #[m/s]

# initial target position
xt0 = (2*np.random.rand(1)-1)*200  #[m] 
yt0 = (2*np.random.rand(1)-1)*200  #[m]

xt0 = xt0[0]
yt0 = yt0[0]

# target maximum speed
w_max = 60*1e3/3600 #[m/s]

# uav flying path curvature constraint
r_min = 400 #[m]

# number of samples for the control search on the boundary
n_sample = 50

# time interval for the cost approximation
Dt = 2 # [seconds]

###############################################
## funciton: optimal target tracking control
############################################### 
#-------------------
def  uav_optimal_tracking_control(state_aircraft,state_tracking):

    xa0 = state_aircraft[0] 
    ya0 = state_aircraft[1]
    vxa0 = state_aircraft[2]
    vya0 = state_aircraft[3]
    ux_min = state_aircraft[4]
    ux_max = state_aircraft[5]
    uy_min = state_aircraft[6]
    uy_max = state_aircraft[7]
    v_min = state_aircraft[8] 
    v_max = state_aircraft[9]
    
    current_speed = np.sqrt(vxa0**2+vya0**2)
    
    xt0 = state_tracking[0]
    yt0 = state_tracking[1]
    w_max = state_tracking[2]
    r_min = state_tracking[3] 
    n_sample = state_tracking[4]
    Dt = state_tracking[5]

    J_cost_uxuy0_function = lambda Dt,ux0,uy0,vxa0,vya0,w_max,xa0,xt0,ya0,yt0: 1.0*Dt**4*ux0**2 + 1.0*Dt**4*uy0**2 + 4.0*Dt**3*ux0*vxa0 + 4.0*Dt**3*uy0*vya0 + 2.0*Dt**2*ux0*xa0 - 2.0*Dt**2*ux0*xt0 - 0.333333333333333*Dt**2*ux0*(2.0*Dt**3*ux0*w_max + 6.0*Dt**2*vxa0*w_max + 4.0*Dt*w_max*xa0 - 4*Dt*w_max*xt0)/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 2.0*Dt**2*uy0*ya0 - 2.0*Dt**2*uy0*yt0 - 0.333333333333333*Dt**2*uy0*(2.0*Dt**3*uy0*w_max + 6.0*Dt**2*vya0*w_max + 4.0*Dt*w_max*ya0 - 4*Dt*w_max*yt0)/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 5.0*Dt**2*vxa0**2 + 5.0*Dt**2*vya0**2 + 6.0*Dt*vxa0*xa0 - 6.0*Dt*vxa0*xt0 - 1.0*Dt*vxa0*(2.0*Dt**3*ux0*w_max + 6.0*Dt**2*vxa0*w_max + 4.0*Dt*w_max*xa0 - 4*Dt*w_max*xt0)/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 6.0*Dt*vya0*ya0 - 6.0*Dt*vya0*yt0 - 1.0*Dt*vya0*(2.0*Dt**3*uy0*w_max + 6.0*Dt**2*vya0*w_max + 4.0*Dt*w_max*ya0 - 4*Dt*w_max*yt0)/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 2.0*xa0**2 - 4.0*xa0*xt0 - 0.666666666666667*xa0*(2.0*Dt**3*ux0*w_max + 6.0*Dt**2*vxa0*w_max + 4.0*Dt*w_max*xa0 - 4*Dt*w_max*xt0)/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 2*xt0**2 + 0.666666666666667*xt0*(2.0*Dt**3*ux0*w_max + 6.0*Dt**2*vxa0*w_max + 4.0*Dt*w_max*xa0 - 4*Dt*w_max*xt0)/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 2.0*ya0**2 - 4.0*ya0*yt0 - 0.666666666666667*ya0*(2.0*Dt**3*uy0*w_max + 6.0*Dt**2*vya0*w_max + 4.0*Dt*w_max*ya0 - 4*Dt*w_max*yt0)/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 2*yt0**2 + 0.666666666666667*yt0*(2.0*Dt**3*uy0*w_max + 6.0*Dt**2*vya0*w_max + 4.0*Dt*w_max*ya0 - 4*Dt*w_max*yt0)/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 2.0*(0.333333333333333*Dt**3*ux0*w_max + Dt**2*vxa0*w_max + 0.666666666666667*Dt*w_max*xa0 - 0.666666666666667*Dt*w_max*xt0)**2/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 2.0*(0.333333333333333*Dt**3*uy0*w_max + Dt**2*vya0*w_max + 0.666666666666667*Dt*w_max*ya0 - 0.666666666666667*Dt*w_max*yt0)**2/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)
        
    # body and global coordinates transfrom dcm
    th_flight = np.arctan2(vya0,vxa0);
    dcm_from_body_to_global = np.array([[np.cos(th_flight), -np.sin(th_flight)], 
                                        [np.sin(th_flight), np.cos(th_flight)]])

    # check the curvature constraint in the body frame
    u_curvature = current_speed**2/r_min
    if u_curvature < uy_max:
        # active constraint & replace the uy bound
        uy_max = u_curvature
        uy_min = -u_curvature
    
    # active the maximum velocity constraint
    vmax_active = False
    if current_speed/Dt+ux_max > v_max/Dt:
        ux_max = -current_speed/Dt+np.sqrt((v_max/Dt)**2-uy_max**2)
        vmax_active = True
    
    # active the minimum velocity constraint
    vmin_active = False
    if current_speed/Dt+ux_min < v_min/Dt:
        ux_min = -current_speed/Dt+np.sqrt((v_min/Dt)**2-uy_max**2)
        vmin_active = True

    # find the optimal solution along the boundary
    n_sample = 50
    ux_sample = np.linspace(ux_min, ux_max, n_sample)
    upper_line = np.vstack((ux_sample,np.ones(n_sample)*uy_max))
    lower_line = np.vstack((ux_sample,np.ones(n_sample)*uy_min))
    
    if vmax_active:
        ux_sample = np.linspace(ux_max,(v_max-current_speed)/Dt,n_sample)
        uy_sample = np.sqrt((v_max/Dt)**2-(ux_sample+current_speed/Dt)**2)
        right_line = np.vstack((np.hstack((ux_sample,np.flip(ux_sample))), np.hstack((uy_sample,-np.flip(uy_sample)))))
    else:
        uy_sample = np.linspace(uy_min,uy_max,n_sample)
        right_line = np.vstack((np.ones(n_sample)*ux_max, uy_sample))
        
    if vmin_active:
        ux_sample = np.linspace(ux_min,(v_min-current_speed)/Dt,n_sample)
        uy_sample = np.sqrt((v_min/Dt)**2-(ux_sample+current_speed/Dt)**2)
        left_line = np.vstack((np.hstack((ux_sample,np.flip(ux_sample))), np.hstack((uy_sample,-1*np.flip(uy_sample)))))
    else:
        uy_sample = np.linspace(uy_min,uy_max,n_sample)
        left_line = np.vstack((np.ones(n_sample)*ux_min, uy_sample))
    
    
    all_samples_in_body_frame = np.hstack((upper_line,lower_line,right_line,left_line))
    all_samples_in_global_frame = dcm_from_body_to_global@all_samples_in_body_frame
    ux0_sample = all_samples_in_global_frame[0,:]
    uy0_sample = all_samples_in_global_frame[1,:]
    J_val = J_cost_uxuy0_function(Dt,ux0_sample,uy0_sample,vxa0,vya0,w_max,xa0,xt0,ya0,yt0)
    
    J_val_opt = J_val.min()
    opt_idx = J_val.argmin()
    uxy_opt_body = all_samples_in_body_frame[:,opt_idx]
    uxy_opt_global = all_samples_in_global_frame[:,opt_idx]
    
    # check the cost function inside the constraint
    polygon_points = np.vstack((ux0_sample,uy0_sample))
    polygon_center = polygon_points.mean(axis=1)
    pc_vector = polygon_points-polygon_center[:,np.newaxis]
    th_pc = np.arctan2(pc_vector[1,:],pc_vector[0,:])
    idx_pc = th_pc.argsort()
    polygon_points = polygon_points[:,idx_pc]
    polygon = path.Path(polygon_points.transpose())
    
    n_inside_sample = 1000;
    x_sample = polygon_points[0,:].min() \
        + (polygon_points[0,:].max()-polygon_points[0,:].min())*np.random.rand(n_inside_sample)
    y_sample = polygon_points[1,:].min() \
        + (polygon_points[1,:].max()-polygon_points[1,:].min())*np.random.rand(n_inside_sample)
    xy_sample=np.vstack((x_sample,y_sample)).transpose()
    
    in_out = polygon.contains_points(xy_sample)
    x_sample = x_sample[in_out]
    y_sample = y_sample[in_out]
    J_val_inside = J_cost_uxuy0_function(Dt,x_sample,y_sample,vxa0,vya0,w_max,xa0,xt0,ya0,yt0)
    J_val_inside = J_val_inside[J_val_inside<J_val_opt]

    if J_val_inside.shape[0]!=0:
        J_val_opt = J_val_inside.min()
        min_idx = J_val_inside.argmin()
          
        dJdux0 = lambda Dt,ux0,uy0,vxa0,vya0,w_max,xa0,xt0,ya0,yt0: -0.666666666666667*Dt**5*ux0*w_max/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 2.0*Dt**4*ux0 + 0.111111111111111*Dt**4*ux0*(0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)*(2.0*Dt**3*ux0*w_max + 6.0*Dt**2*vxa0*w_max + 4.0*Dt*w_max*xa0 - 4*Dt*w_max*xt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**(3/2) + 0.111111111111111*Dt**4*uy0*(0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)*(2.0*Dt**3*uy0*w_max + 6.0*Dt**2*vya0*w_max + 4.0*Dt*w_max*ya0 - 4*Dt*w_max*yt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**(3/2) - 2.0*Dt**4*vxa0*w_max/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 4.0*Dt**3*vxa0 + 0.333333333333333*Dt**3*vxa0*(0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)*(2.0*Dt**3*ux0*w_max + 6.0*Dt**2*vxa0*w_max + 4.0*Dt*w_max*xa0 - 4*Dt*w_max*xt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**(3/2) + 0.333333333333333*Dt**3*vya0*(0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)*(2.0*Dt**3*uy0*w_max + 6.0*Dt**2*vya0*w_max + 4.0*Dt*w_max*ya0 - 4*Dt*w_max*yt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**(3/2) - 1.33333333333333*Dt**3*w_max*xa0/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 1.33333333333333*Dt**3*w_max*xt0/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 1.33333333333333*Dt**3*w_max*(0.333333333333333*Dt**3*ux0*w_max + Dt**2*vxa0*w_max + 0.666666666666667*Dt*w_max*xa0 - 0.666666666666667*Dt*w_max*xt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 2.0*Dt**2*xa0 + 0.222222222222222*Dt**2*xa0*(0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)*(2.0*Dt**3*ux0*w_max + 6.0*Dt**2*vxa0*w_max + 4.0*Dt*w_max*xa0 - 4*Dt*w_max*xt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**(3/2) - 2.0*Dt**2*xt0 - 0.222222222222222*Dt**2*xt0*(0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)*(2.0*Dt**3*ux0*w_max + 6.0*Dt**2*vxa0*w_max + 4.0*Dt*w_max*xa0 - 4*Dt*w_max*xt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**(3/2) + 0.222222222222222*Dt**2*ya0*(0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)*(2.0*Dt**3*uy0*w_max + 6.0*Dt**2*vya0*w_max + 4.0*Dt*w_max*ya0 - 4*Dt*w_max*yt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**(3/2) - 0.222222222222222*Dt**2*yt0*(0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)*(2.0*Dt**3*uy0*w_max + 6.0*Dt**2*vya0*w_max + 4.0*Dt*w_max*ya0 - 4*Dt*w_max*yt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**(3/2) - 1.33333333333333*Dt**2*(0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)*(0.333333333333333*Dt**3*ux0*w_max + Dt**2*vxa0*w_max + 0.666666666666667*Dt*w_max*xa0 - 0.666666666666667*Dt*w_max*xt0)**2/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**2 - 1.33333333333333*Dt**2*(0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)*(0.333333333333333*Dt**3*uy0*w_max + Dt**2*vya0*w_max + 0.666666666666667*Dt*w_max*ya0 - 0.666666666666667*Dt*w_max*yt0)**2/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**2 - 0.333333333333333*Dt**2*(2.0*Dt**3*ux0*w_max + 6.0*Dt**2*vxa0*w_max + 4.0*Dt*w_max*xa0 - 4*Dt*w_max*xt0)/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)
        dJduy0 = lambda Dt,ux0,uy0,vxa0,vya0,w_max,xa0,xt0,ya0,yt0: -0.666666666666667*Dt**5*uy0*w_max/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 0.111111111111111*Dt**4*ux0*(0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)*(2.0*Dt**3*ux0*w_max + 6.0*Dt**2*vxa0*w_max + 4.0*Dt*w_max*xa0 - 4*Dt*w_max*xt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**(3/2) + 2.0*Dt**4*uy0 + 0.111111111111111*Dt**4*uy0*(0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)*(2.0*Dt**3*uy0*w_max + 6.0*Dt**2*vya0*w_max + 4.0*Dt*w_max*ya0 - 4*Dt*w_max*yt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**(3/2) - 2.0*Dt**4*vya0*w_max/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 0.333333333333333*Dt**3*vxa0*(0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)*(2.0*Dt**3*ux0*w_max + 6.0*Dt**2*vxa0*w_max + 4.0*Dt*w_max*xa0 - 4*Dt*w_max*xt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**(3/2) + 4.0*Dt**3*vya0 + 0.333333333333333*Dt**3*vya0*(0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)*(2.0*Dt**3*uy0*w_max + 6.0*Dt**2*vya0*w_max + 4.0*Dt*w_max*ya0 - 4*Dt*w_max*yt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**(3/2) - 1.33333333333333*Dt**3*w_max*ya0/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 1.33333333333333*Dt**3*w_max*yt0/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 1.33333333333333*Dt**3*w_max*(0.333333333333333*Dt**3*uy0*w_max + Dt**2*vya0*w_max + 0.666666666666667*Dt*w_max*ya0 - 0.666666666666667*Dt*w_max*yt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2) + 0.222222222222222*Dt**2*xa0*(0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)*(2.0*Dt**3*ux0*w_max + 6.0*Dt**2*vxa0*w_max + 4.0*Dt*w_max*xa0 - 4*Dt*w_max*xt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**(3/2) - 0.222222222222222*Dt**2*xt0*(0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)*(2.0*Dt**3*ux0*w_max + 6.0*Dt**2*vxa0*w_max + 4.0*Dt*w_max*xa0 - 4*Dt*w_max*xt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**(3/2) + 2.0*Dt**2*ya0 + 0.222222222222222*Dt**2*ya0*(0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)*(2.0*Dt**3*uy0*w_max + 6.0*Dt**2*vya0*w_max + 4.0*Dt*w_max*ya0 - 4*Dt*w_max*yt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**(3/2) - 2.0*Dt**2*yt0 - 0.222222222222222*Dt**2*yt0*(0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)*(2.0*Dt**3*uy0*w_max + 6.0*Dt**2*vya0*w_max + 4.0*Dt*w_max*ya0 - 4*Dt*w_max*yt0)/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**(3/2) - 1.33333333333333*Dt**2*(0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)*(0.333333333333333*Dt**3*ux0*w_max + Dt**2*vxa0*w_max + 0.666666666666667*Dt*w_max*xa0 - 0.666666666666667*Dt*w_max*xt0)**2/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**2 - 1.33333333333333*Dt**2*(0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)*(0.333333333333333*Dt**3*uy0*w_max + Dt**2*vya0*w_max + 0.666666666666667*Dt*w_max*ya0 - 0.666666666666667*Dt*w_max*yt0)**2/((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)**2 - 0.333333333333333*Dt**2*(2.0*Dt**3*uy0*w_max + 6.0*Dt**2*vya0*w_max + 4.0*Dt*w_max*ya0 - 4*Dt*w_max*yt0)/sqrt((0.333333333333333*Dt**2*ux0 + Dt*vxa0 + 0.666666666666667*xa0 - 0.666666666666667*xt0)**2 + (0.333333333333333*Dt**2*uy0 + Dt*vya0 + 0.666666666666667*ya0 - 0.666666666666667*yt0)**2)
            
        s_amj = 0.01
        alpha_amj = s_amj; beta_amj = 0.5; sigma_amj = 1e-5
        u_xy_current = np.array([ux0_sample[min_idx], uy0_sample[min_idx]])
        J_current = J_cost_uxuy0_function(Dt,u_xy_current[0],u_xy_current[1],vxa0,vya0,w_max,xa0,xt0,ya0,yt0)
        dJdu = np.array([dJdux0(Dt,u_xy_current[0],u_xy_current[1],vxa0,vya0,w_max,xa0,xt0,ya0,yt0),
                         dJduy0(Dt,u_xy_current[0],u_xy_current[1],vxa0,vya0,w_max,xa0,xt0,ya0,yt0)])
        while True:
            u_xy_update = u_xy_current - alpha_amj*dJdu
            J_update = J_cost_uxuy0_function(Dt,u_xy_update[0],u_xy_update[1],vxa0,vya0,w_max,xa0,xt0,ya0,yt0)
            if J_update < (J_current + sigma_amj*alpha_amj*np.sum(dJdu**2)):
                if np.linalg.norm(u_xy_current-u_xy_update)<1e-6:
                    break
                alpha_amj = s_amj
                J_current = J_cost_uxuy0_function(Dt,u_xy_update[0],u_xy_update[1],vxa0,vya0,w_max,xa0,xt0,ya0,yt0)
                dJdu = np.array([dJdux0(Dt,u_xy_update[0],u_xy_update[1],vxa0,vya0,w_max,xa0,xt0,ya0,yt0),
                                 dJduy0(Dt,u_xy_update[0],u_xy_update[1],vxa0,vya0,w_max,xa0,xt0,ya0,yt0)])
                u_xy_current = u_xy_update
            else:
                alpha_amj = beta_amj*alpha_amj
      
        uxy_opt_global = u_xy_current
        uxy_opt_body = dcm_from_body_to_global.T@uxy_opt_global
        
    return uxy_opt_global, uxy_opt_body
        
###############################################
## main part of monte-carlo simulations
###############################################
N_sim = 1800

# optimal control input error
opt_costh_err = np.zeros(N_sim-1)
opt_mag_err = np.zeros(N_sim-1)


# aircraft & target dynamics
Fa = np.eye(4)+np.vstack((np.hstack((np.zeros((2,2)),Dt*np.eye(2))),np.zeros((2,4)))) 
Ga = np.vstack((np.zeros((2,2)),Dt*np.eye(2)))

Ft = np.eye(2)
Gt = Dt*np.eye(2)

state_uav = np.array([xa0, ya0, vxa0, vya0]).squeeze()
state_target = np.array([xt0, yt0]).squeeze()

uav_pos = np.array([xa0, ya0])
target_pos = np.array([xt0, yt0])

uxy_opt_global = np.array([0, 0])
w_zero_one = np.zeros((2,2))


fig, ax = plt.subplots(nrows=1,ncols=1)
ax.set(xlim=(-2000, 2000), ylim=(-2000, 2000))
line1, = ax.plot(xa0, ya0,'bo')
line2, = ax.plot(xt0, yt0,'ro')

for idx_sim in np.arange(N_sim):
    
    print(f'{idx_sim:3d}/{N_sim:3d}\n')

    # uav optimal input
    state_aircraft = np.array([xa0,ya0,vxa0,vya0,ux_min,ux_max,uy_min,uy_max,v_min,v_max])
    state_tracking = np.array([xt0,yt0,w_max,r_min,n_sample,Dt])
    
    uxy_opt_global, uxy_opt_body = uav_optimal_tracking_control(state_aircraft,state_tracking)

    # target maximum random direction input
    rand_th = 2*np.pi*np.random.rand(1)
    wx0 = w_max*np.cos(rand_th)  
    wy0 = w_max*np.sin(rand_th)
    
    # state propagation
    state_uav = Fa@state_uav + Ga@uxy_opt_global
    state_target = Ft@state_target + Gt@(np.array([wx0, wy0]).squeeze())
    
    # reset state variables
    xa0 = state_uav[0]
    ya0 = state_uav[1]
    vxa0 = state_uav[2]
    vya0 = state_uav[3]
    xt0 = state_target[0]
    yt0 = state_target[1]
    
    line1.set_xdata(xa0)
    line1.set_ydata(ya0)
    line2.set_xdata(xt0)
    line2.set_ydata(yt0)
    
    fig.canvas.draw()
    fig.canvas.flush_events()
  
    time.sleep(0.01)

    
#     % update the plot
#     figure(1); 
#     addpoints(uav_line,xa0,ya0);
#     hold on;
#     plot(xt0,yt0,'r.');
#     hold off;
#     %F(idx_sim) = getframe(gcf);
#     pause(0.01);
