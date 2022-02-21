#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 14:47:53 2021

@author: jongrae
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import path

from sympy import symbols, simplify, expand
from sympy import cos, sin, sqrt, diff
from sympy.utilities.lambdify import lambdify

import time
from scipy.optimize import minimize, fsolve

ux0, uy0, ux1, uy1, wx0, wy0, wx1, wy1 = symbols('ux0 uy0 ux1 uy1 wx0 wy0 wx1 wy1', real=True)
xa0, ya0, vxa0, vya0, xt0, yt0, th = symbols('xa0 ya0 vxa0 vya0 xt0 yt0 th', real=True)

Dt, w_max = symbols('Dt w_max', real=True, positive=True)

#----------------------------------------------
# Dynamics
#----------------------------------------------
Fa = np.eye(4)+np.vstack((np.hstack((np.zeros((2,2)),Dt*np.eye(2))),np.zeros((2,4)))) 
Ga = np.vstack((np.zeros((2,2)),Dt*np.eye(2)))
Ca = np.eye(2,4)

Ft = np.eye(2)
Gt = Dt*np.eye(2)
Ct = np.eye(2)

# control inputs
u_vec_0 = np.array([[ux0], [uy0]])
w_vec_0 = np.array([[wx0], [wy0]])
u_vec_1 = np.array([[ux1], [uy1]])
w_vec_1 = np.array([[wx1], [wy1]])

# initial conditions
xa_vec_0 = np.array([[xa0], [ya0], [vxa0], [vya0]])
xt_vec_0 = np.array([[xt0], [yt0]])

# state propagation
xa_k_plus_1 = Fa@xa_vec_0    + Ga@u_vec_0
xa_k_plus_2 = Fa@xa_k_plus_1 + Ga@u_vec_1;
y_k_plus_1 = Ca@xa_k_plus_1;
y_k_plus_2 = Ca@xa_k_plus_2;

xt_k_plus_1 = Ft@xt_vec_0    + Gt@w_vec_0
xt_k_plus_2 = Ft@xt_k_plus_1 + Gt@w_vec_1
z_k_plus_1 = Ct@xt_k_plus_1
z_k_plus_2 = Ct@xt_k_plus_2

#----------------------------------------------
# calculate the cost function with the worst target manoeuvre
#----------------------------------------------
xa1 = y_k_plus_1[0][0]
ya1 = y_k_plus_1[1][0]
xa2 = y_k_plus_2[0][0]
ya2 = y_k_plus_2[1][0]

r_T0A1 = np.array([[xt0 - xa1], [yt0 - ya1]])
Delta_rt_0 = np.array([[Dt*w_max*cos(th)], [Dt*w_max*sin(th)]])
r_A2A1 = np.array([[xa2-xa1], [ya2-ya1]])

ell_1 = r_T0A1 + Delta_rt_0
ell_1_squared = (ell_1.T@ell_1)[0][0]

r_T1A2 = ell_1 - r_A2A1
ell_2_squared = (r_T1A2.T@r_T1A2)[0][0]

J_cost_worst = expand(ell_1_squared + ell_2_squared)
dJdth_worst = expand(simplify(diff(J_cost_worst,th)))
coeff_cos = dJdth_worst.coeff(cos(th))
coeff_sin = dJdth_worst.coeff(sin(th))

# calculate the worst coast function
a_triangle = -coeff_cos # do not forget the minus sign
b_triangle = coeff_sin
c_triangle = simplify(sqrt(a_triangle**2 + b_triangle**2))
check_a_b = float(expand(-a_triangle*cos(th)+b_triangle*sin(th)-dJdth_worst))
print(f'Check [-a*cos(th)+b*sin(th)]-dJdth_worst equal to zero? {check_a_b:4.2f}')

J_cost_worst = J_cost_worst.subs(sin(th),a_triangle/c_triangle)
J_cost_worst = J_cost_worst.subs(cos(th),b_triangle/c_triangle)

dJdux0 = diff(J_cost_worst,ux0)
dJduy0 = diff(J_cost_worst,uy0)

#---------------------------------------------------
# evaluate the cost function for test scenario values
#---------------------------------------------------

# initial target position
xt0_v = (2*np.random.rand(1)-1)*200*0+150  #[m] 
yt0_v = (2*np.random.rand(1)-1)*200*0  #[m]

# initial uav position
xa0_v = (2*np.random.rand(1)-1)*100*0 #[m] 
ya0_v = (2*np.random.rand(1)-1)*100*0 #[m]

# initial uav velocity
tha0 = np.random.rand(1)*2*np.pi*0 #[radian]
current_speed = 25 #[m/s]
vxa0_v = current_speed*np.cos(tha0)
vya0_v = current_speed*np.sin(tha0)

# uav minimum & maximum speed
v_min = 20  #[m/s]
v_max = 40  #[m/s]

# time interval for the cost approximation
Dt_v = 2 # [seconds]

# target maximum speed
w_max_v = 60*1e3/3600 #[m/s]

# uav flying path curvature constraint
r_min = 400 #[m]

# control acceleration input magnitude constraints
ux_max = 10 # [m/s^2]
ux_min = -1 # [m/s^2]
uy_max = 2  # [m/s^2]
uy_min = -2 # [m/s^2]

ux_max_org = ux_max
ux_min_org = ux_min
uy_max_org = uy_max
uy_min_org = uy_min

# evaluate the cost function over the ux0-uy0 control input
num_idx = 20
num_jdx = 19
min_max_u_plot = 20
ux_all = np.linspace(-min_max_u_plot,min_max_u_plot,num_idx)
uy_all = np.linspace(-min_max_u_plot,min_max_u_plot,num_jdx)

values = [(Dt,Dt_v), (xa0,xa0_v[0]), (ya0,ya0_v[0]), (vxa0,vxa0_v[0]), (vya0,vya0_v[0]), 
          (xt0,xt0_v[0]), (yt0,yt0_v[0]), (w_max,w_max_v)]

J_cost_uxuy0 = J_cost_worst.subs(values)
J_cost_uxuy0_function = lambdify([ux0,uy0],J_cost_uxuy0)

UX0,UY0=np.meshgrid(ux_all,uy_all)
J_cost_worst_val=J_cost_uxuy0_function(UX0,UY0)

# replace symbols by values
Dt = Dt_v
xa0 = xa0_v[0]
ya0 = ya0_v[0]
vxa0 = vxa0_v[0]
vya0 = vya0_v[0]
xt0 = xt0_v[0]
yt0 = yt0_v[0]
w_max = w_max_v

#########################################################
## Optimal control input
#---------------------------
# find optimal control
#---------------------------

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

th_flight = np.arctan2(vya0,vxa0)
dcm_from_body_to_global = np.array([
    [np.cos(th_flight), -np.sin(th_flight)], 
    [np.sin(th_flight), np.cos(th_flight)]])

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
J_val = J_cost_uxuy0_function(ux0_sample,uy0_sample)

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
J_val_inside = J_cost_uxuy0_function(x_sample,y_sample)
J_val_inside = J_val_inside[J_val_inside<J_val_opt]


if J_val_inside.shape[0]!=0:
    J_val_opt = J_val_inside.min()
    min_idx = J_val_inside.argmin()
    
    t0 = time.time()
    J_cost_minimize=lambda x: J_cost_uxuy0_function(x[0],x[1])
    sol_opt=minimize(J_cost_minimize,[ux0_sample[min_idx],uy0_sample[min_idx]],method='BFGS')
    uxy_opt_global_1 = sol_opt.x
    tf = time.time() - t0
    print(f'minimization: {tf:10.8f} [s]\n')
    
    t0 = time.time()
    dJdux0_fun = lambdify([ux0,uy0],dJdux0.subs(values))
    dJduy0_fun = lambdify([ux0,uy0],dJduy0.subs(values))
    dJduxy=lambda x: np.array([dJdux0_fun(x[0],x[1]), dJduy0_fun(x[0],x[1])])
    uxy_opt_global_2 = fsolve(dJduxy,[ux0_sample[min_idx],uy0_sample[min_idx]])
    tf = time.time() - t0
    print(f'fsolve: {tf:10.8f} [s]\n')
    
    t0 = time.time()
    s_amj = 0.01
    alpha_amj = s_amj; beta_amj = 0.5; sigma_amj = 1e-5
    u_xy_current = np.array([ux0_sample[min_idx], uy0_sample[min_idx]])
    J_current = J_cost_minimize(u_xy_current)
    dJdu = dJduxy(u_xy_current)
    while True:
        u_xy_update = u_xy_current - alpha_amj*dJdu
        J_update = J_cost_minimize(u_xy_update)
        if J_update < (J_current + sigma_amj*alpha_amj*np.sum(dJdu**2)):
            if np.linalg.norm(u_xy_current-u_xy_update)<1e-6:
                break
            alpha_amj = s_amj
            J_current = J_cost_minimize(u_xy_update)
            dJdu = dJduxy(u_xy_update)
            u_xy_current = u_xy_update
        else:
            alpha_amj = beta_amj*alpha_amj
        
    tf = time.time() - t0
    print(f'Gradient Descent with Armijo\'s Rule: {tf:10.8f} [s]\n')
   
    uxy_opt_global = u_xy_current
    uxy_opt_body = dcm_from_body_to_global.T@uxy_opt_global
    
    print(uxy_opt_global_1)
    print(uxy_opt_global_2)
    print(uxy_opt_global)
    

##############################################################
## plot
#-------------------------------------
# plot the cost and the constraints
#-------------------------------------
fig, ax = plt.subplots(nrows=1,ncols=1)
ax.contourf(UX0,UY0,J_cost_worst_val)
ax.set_ylabel(r'$u_y(0)$',fontsize=14);
ax.set_xlabel(r'$u_x(0)$',fontsize=14);
ax.set(xlim=(-min_max_u_plot, min_max_u_plot), ylim=(-min_max_u_plot, min_max_u_plot))
ax.set_aspect('equal','box')

# velocity constraints
xc = -vxa0/Dt
yc = -vya0/Dt
th_plot = np.linspace(0,2*np.pi,100)
ax.plot(xc+v_min/Dt*np.cos(th_plot),yc+v_min/Dt*np.sin(th_plot),'r-.')
ax.plot(xc+v_max/Dt*np.cos(th_plot),yc+v_max/Dt*np.sin(th_plot),'r');

# uav flying direction arrow
v_a0_normalise = np.array([vxa0, vya0])
v_a0_normalise = 10*v_a0_normalise/np.linalg.norm(v_a0_normalise)
v_a_arrow=ax.quiver(0,0,v_a0_normalise[0],v_a0_normalise[1],scale=50,color='magenta')
ax.text(v_a0_normalise[0],v_a0_normalise[1],'aircarft velocity (magenta)')

# target location direction arrow
dr_t0a0 = 0.3*(np.array([xt0, yt0]) - np.array([xa0, ya0]))
dr_t0a0 = 15*dr_t0a0/np.linalg.norm(dr_t0a0)
r_target_arrow=ax.quiver(0,0,dr_t0a0[0],dr_t0a0[1],scale=50,color='cyan',)
ax.text(dr_t0a0[0],dr_t0a0[1],'target direction (cyan)')

# curvature limit
cd_cvt = ((vxa0**2+vya0**2)**(1.5))/r_min
if np.abs(vxa0) > 0.5*v_max:
    m_cvt = vya0/vxa0
    cd_cvt = cd_cvt/vxa0
    ux_cvt = np.linspace(-30,30,100)
    uy_cvt = m_cvt*ux_cvt
    
    ux_cvt_line_1 = ux_cvt
    uy_cvt_line_1 = uy_cvt + cd_cvt
    ux_cvt_line_2 = ux_cvt
    uy_cvt_line_2 = uy_cvt - cd_cvt
else:
    n_cvt = vxa0/vya0
    cd_cvt = cd_cvt/vya0
    uy_cvt = np.linspace(-30,30,100)
    ux_cvt = n_cvt*uy_cvt
    
    ux_cvt_line_1 = ux_cvt + cd_cvt
    uy_cvt_line_1 = uy_cvt
    
    ux_cvt_line_2 = ux_cvt - cd_cvt
    uy_cvt_line_2 = uy_cvt

ax.plot(ux_cvt_line_1,uy_cvt_line_1,'b',linewidth=2)
ax.plot(ux_cvt_line_2,uy_cvt_line_2,'b',linewidth=2)

ax.plot(all_samples_in_global_frame[0,:],all_samples_in_global_frame[1,:],'k.');
ax.quiver(0,0,uxy_opt_global[0],uxy_opt_global[1],color='red')
ax.text(uxy_opt_global[0],uxy_opt_global[1],'optimal control input (red)')

# draw control input acceleration magnitude contraint box
ru = np.array([ux_max_org, uy_max_org])
rl = np.array([ux_max_org, uy_min_org])
lu = np.array([ux_min_org, uy_max_org])
ll = np.array([ux_min_org, uy_min_org])

ru = dcm_from_body_to_global@ru
rl = dcm_from_body_to_global@rl
lu = dcm_from_body_to_global@lu
ll = dcm_from_body_to_global@ll

ax.plot([ru[0], rl[0]], [ru[1], rl[1]],'g-.',linewidth=3)
ax.plot([ru[0], lu[0]], [ru[1], lu[1]],'g-.',linewidth=3)
ax.plot([lu[0], ll[0]], [lu[1], ll[1]],'g-.',linewidth=3)
ax.plot([rl[0], ll[0]], [rl[1], ll[1]],'g-.',linewidth=3)


ax.legend([r'$v_{\rm min}$',r'$v_{\rm max}$', 'curvature constraints', 'curvature constraints',
          'control input samples','control input acceleration bound'], 
          title='contours: cost function',bbox_to_anchor=(1.05, 1), loc='upper left',)
          