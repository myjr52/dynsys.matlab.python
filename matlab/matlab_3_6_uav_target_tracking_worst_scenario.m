% MIT License
% 
% Copyright (c) 2022 Jongrae.K
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

clear;

%----------------------------------------------------------------------
% define symbols
%----------------------------------------------------------------------
% define time interval
syms Dt real;

% define symbols for aircraft's and target's control inputs
syms ux0 uy0 ux1 uy1 real;
syms wx0 wy0 wx1 wy1 real;

% define symbols for the initial conditions
syms xa0 ya0 vxa0 vya0 real;
syms xt0 yt0 real; 

% target characteristics 
syms th w_max real;

%----------------------------------------------------------------------
% Dynamics
%----------------------------------------------------------------------
% aircraft & target dynamics
Fa = eye(4) + [zeros(2) Dt*eye(2); zeros(2,4)];
Ga = [zeros(2); Dt*eye(2)];
Ca = eye(2,4);

Ft = eye(2);
Gt = Dt*eye(2);
Ct = eye(2);

% control inputs
u_vec_0 = [ux0 uy0]';
w_vec_0 = [wx0 wy0]';
u_vec_1 = [ux1 uy1]';
w_vec_1 = [wx1 wy1]';

% initial conditions
xa_vec_0 = [xa0 ya0 vxa0 vya0]';
xt_vec_0 = [xt0 yt0]';

% state propagation
xa_k_plus_1 = Fa*xa_vec_0    + Ga*u_vec_0;
xa_k_plus_2 = Fa*xa_k_plus_1 + Ga*u_vec_1;
y_k_plus_1 = Ca*xa_k_plus_1;
y_k_plus_2 = Ca*xa_k_plus_2;

xt_k_plus_1 = Ft*xt_vec_0    + Gt*w_vec_0;
xt_k_plus_2 = Ft*xt_k_plus_1 + Gt*w_vec_1;
z_k_plus_1 = Ct*xt_k_plus_1;
z_k_plus_2 = Ct*xt_k_plus_2;

%----------------------------------------------------------------------
% calculate the cost function with the worst target manoeuvre
%----------------------------------------------------------------------
xa1 = y_k_plus_1(1);
ya1 = y_k_plus_1(2);
xa2 = y_k_plus_2(1);
ya2 = y_k_plus_2(2);

r_T0A1 = [xt0 - xa1; yt0 - ya1];
Delta_rt_0 = [Dt*w_max*cos(th); Dt*w_max*sin(th)];
r_A2A1 = [xa2-xa1; ya2-ya1];

ell_1 = r_T0A1 + Delta_rt_0;
ell_1_squared = ell_1(:)'*ell_1(:);

r_T1A2 = ell_1 - r_A2A1;
ell_2_squared = r_T1A2(:)'*r_T1A2(:);% + (Dt*w_max)^2;

J_cost_worst = (ell_1_squared + ell_2_squared);
dJdth_worst = simplify(diff(J_cost_worst,th));

% calculate the worst cost function
a_triangle = coeffs(dJdth_worst,cos(th));
a_triangle = -a_triangle(2); % do not forget the minus
b_triangle = coeffs(dJdth_worst,sin(th));
b_triangle = b_triangle(2);
check_a_b = expand(-a_triangle*cos(th)+b_triangle*sin(th)-dJdth_worst);
fprintf('Check [-a*cos(th)+b*sin(th)]-dJdth_worst equal to zero? %4.2f\n', check_a_b);
 
c_triangle = sqrt(a_triangle^2 + b_triangle^2);
J_cost_worst = eval(J_cost_worst);
J_cost_worst = subs(J_cost_worst,sin(th),a_triangle/c_triangle);
J_cost_worst = subs(J_cost_worst,cos(th),b_triangle/c_triangle); 
J_cost_worst = simplify(expand(J_cost_worst));

dJdux0 = simplify(diff(J_cost_worst,ux0));
dJduy0 = simplify(diff(J_cost_worst,uy0));
%----------------------------------------------------------------------
% evaluate the cost function for test scenario values
%----------------------------------------------------------------------

% initial target position
xt0 = (2*rand(1)-1)*200;  %[m] 
yt0 = (2*rand(1)-1)*200;  %[m]

% initial uav position
xa0 = (2*rand(1)-1)*100; %[m] 
ya0 = (2*rand(1)-1)*100; %[m]

% initial uav velocity
tha0 = rand(1)*2*pi; %[radian]
current_speed = 25; %[m/s]
vxa0 = current_speed*cos(tha0); 
vya0 = current_speed*sin(tha0);

% uav minimum & maximum speed
v_min = 20; v_max = 40; %[m/s]

% time interval for the cost approximation
Dt = 2; % [seconds]

% target maximum speed
w_max = 60*1e3/3600; %[m/s]

% uav flying path curvature constraint
r_min = 400; %[m]

% control acceleration input magnitude constraints
ux_max = 10; % [m/s^2]
ux_min = -1; % [m/s^2]
uy_max = 2;  % [m/s^2]
uy_min = -2; % [m/s^2]

ux_max_org = ux_max;
ux_min_org = ux_min;
uy_max_org = uy_max;
uy_min_org = uy_min;

% evaluate the cost function over the ux0-uy0 control input
num_idx = 20;
num_jdx = 19;
min_max_u_plot = 20;
ux_all = linspace(-min_max_u_plot,min_max_u_plot,num_idx);
uy_all = linspace(-min_max_u_plot,min_max_u_plot,num_jdx);

% ux_all = linspace(26.8,28,num_idx*10);
% uy_all = linspace(17,19,num_jdx*10);


J_cost_uxuy0 = eval(J_cost_worst);
J_cost_uxuy0_function = matlabFunction(J_cost_uxuy0);

[UX0,UY0]=meshgrid(ux_all,uy_all);
J_cost_worst_val=J_cost_uxuy0_function(UX0,UY0);


% dJdux0=eval(dJdux0);
% dJduy0=eval(dJduy0);
% uxy_opt_no_constraint = solve([dJdux0==0, dJduy0==0],[ux0 uy0]);


%% Optimal control input
%-------------------------------------------------------------------------
% find optimal control
%-------------------------------------------------------------------------

% check the curvature constraint in the body frame
u_curvature = current_speed^2/r_min;
if u_curvature < uy_max
    % active constraint & replace the uy bound
    uy_max = u_curvature;
    uy_min = -u_curvature;
end

% active the maximum velocity constraint
vmax_active = false;
if current_speed/Dt+ux_max > v_max/Dt
    ux_max = -current_speed/Dt+sqrt((v_max/Dt)^2-uy_max^2);
    vmax_active = true;
end

% active the minimum velocity constraint
vmin_active = false;
if current_speed/Dt+ux_min < v_min/Dt
    ux_min = -current_speed/Dt+sqrt((v_min/Dt)^2-uy_max^2);
    vmin_active = true;
end

th_flight = atan2(vya0,vxa0);
dcm_from_body_to_global = [cos(th_flight) -sin(th_flight); sin(th_flight) cos(th_flight)];

% find the optimal solution along the boundary
n_sample = 50;
ux_sample = linspace(ux_min, ux_max, n_sample);
upper_line = [ux_sample; ones(1,n_sample)*uy_max];
lower_line = [ux_sample; ones(1,n_sample)*uy_min];

if vmax_active
    ux_sample = linspace(ux_max,(v_max-current_speed)/Dt,n_sample);
    uy_sample = sqrt((v_max/Dt)^2-(ux_sample+current_speed/Dt).^2);
    right_line = [ux_sample ux_sample(end-1:-1:1); uy_sample -uy_sample(end-1:-1:1)];
else
    uy_sample = linspace(uy_min,uy_max,n_sample);
    right_line = [ones(1,n_sample)*ux_max; uy_sample];
end

if vmin_active
    ux_sample = linspace(ux_min,(v_min-current_speed)/Dt,n_sample);
    uy_sample = sqrt((v_min/Dt)^2-(ux_sample+current_speed/Dt).^2);
    left_line = [ux_sample ux_sample(end-1:-1:1); uy_sample -uy_sample(end-1:-1:1)];
else
    uy_sample = linspace(uy_min,uy_max,n_sample);
    left_line = [ones(1,n_sample)*ux_min; uy_sample];
end

all_samples_in_body_frame = [upper_line lower_line right_line left_line];
all_samples_in_global_frame = dcm_from_body_to_global*all_samples_in_body_frame;
ux0_sample = all_samples_in_global_frame(1,:);
uy0_sample = all_samples_in_global_frame(2,:);

J_val = J_cost_uxuy0_function(ux0_sample,uy0_sample);

[J_val_opt,opt_idx]=min(J_val);
uxy_opt_body = all_samples_in_body_frame(:,opt_idx);
uxy_opt_global = all_samples_in_global_frame(:,opt_idx);

% check the cost function inside the constraint
polygon_points = [ux0_sample(:)'; uy0_sample(:)'];
polygon_center = mean(polygon_points,2);
pc_vector = polygon_points-polygon_center ;
th_pc = atan2(pc_vector(2,:),pc_vector(1,:));
[~, idx_pc] = sort(th_pc);
polygon_points = polygon_points(:,idx_pc);
polygon_points = [polygon_points polygon_points(:,1)];

n_inside_sample = 1000;
x_sample = min(polygon_points(1,:)) + ...
    (max(polygon_points(1,:))-min(polygon_points(1,:)))*rand(1,n_inside_sample);
y_sample = min(polygon_points(2,:)) + ...
    (max(polygon_points(2,:))-min(polygon_points(2,:)))*rand(1,n_inside_sample);

[in,on] = inpolygon(x_sample,y_sample,polygon_points(1,:),polygon_points(2,:));
x_sample = x_sample(in);
y_sample = y_sample(in);
J_val_inside = J_cost_uxuy0_function(x_sample,y_sample);
J_val_inside = J_val_inside(J_val_inside<J_val_opt);

if ~isempty(J_val_inside)
    [J_val_opt,min_idx] = min(J_val_inside);
    
    tic
    J_cost_minimize=@(x)J_cost_uxuy0_function(x(1),x(2));
    uxy_opt_global_1 = fminunc(J_cost_minimize,[ux0_sample(min_idx) uy0_sample(min_idx)]);
    toc
    
    tic
    dJdux0_fun=matlabFunction(eval(dJdux0));
    dJduy0_fun=matlabFunction(eval(dJduy0));
    dJduxy=@(x)[dJdux0_fun(x(1),x(2)); dJduy0_fun(x(1),x(2))];
    uxy_opt_global_2 = fsolve(dJduxy,[ux0_sample(min_idx) uy0_sample(min_idx)]);
    toc
    
    tic
    s_amj = 0.01;
    alpha_amj = s_amj; beta_amj = 0.5; sigma_amj = 1e-5;
    u_xy_current = [ux0_sample(min_idx) uy0_sample(min_idx)];
    J_current = J_cost_minimize(u_xy_current);
    dJdu = dJduxy(u_xy_current);
    while true
        u_xy_update = u_xy_current - alpha_amj*dJdu(:)';
        J_update = J_cost_minimize(u_xy_update);
        if J_update < (J_current + sigma_amj*alpha_amj*sum(dJdu.^2))
            if norm(u_xy_current-u_xy_update)<1e-6
                break
            end
            alpha_amj = s_amj;
            J_current = J_cost_minimize(u_xy_update);
            dJdu = dJduxy(u_xy_update);
            u_xy_current = u_xy_update;
        else
            alpha_amj = beta_amj*alpha_amj;
        end
    end
    toc
   
    uxy_opt_global = u_xy_current(:);
    uxy_opt_body = dcm_from_body_to_global'*uxy_opt_global;
    
    [   uxy_opt_global_1(:)';
        uxy_opt_global_2(:)';
        uxy_opt_global(:)']
    
end

%% plot
%----------------------------------------------------------------------
% plot the cost and the constraints
%----------------------------------------------------------------------
figure(1); clf;

% the worst cost function values on the ux-uy space 
contourf(ux_all,uy_all,J_cost_worst_val); 
hold on;
axis equal;

% velocity constraints
xc = -vxa0/Dt; yc = -vya0/Dt;
th_plot = linspace(0,2*pi,100);
plot(xc+v_min/Dt*cos(th_plot),yc+v_min/Dt*sin(th_plot),'r-.','LineWidth',1);
plot(xc+v_max/Dt*cos(th_plot),yc+v_max/Dt*sin(th_plot),'r','LineWidth',1);

% uav flying direction arrow
v_a0_normalise = [vxa0 vya0];
v_a0_normalise = 10*v_a0_normalise/norm(v_a0_normalise);
quiver(0,0,v_a0_normalise(1),v_a0_normalise(2),'m','LineWidth',1,'MaxHeadSize',1.5);

% target location direction arrow
dr_t0a0 = 0.3*([xt0 yt0] - [xa0 ya0]);
dr_t0a0 = 15*dr_t0a0/norm(dr_t0a0);
quiver(0,0,dr_t0a0(1),dr_t0a0(2),'c--','LineWidth',1,'MaxHeadSize',1.5);

% curvature limit
cd_cvt = ((vxa0^2+vya0^2)^(1.5))/r_min;
if abs(vxa0) > 0.5*v_max
    m_cvt = vya0/vxa0;
    cd_cvt = cd_cvt/vxa0;
    ux_cvt = linspace(-30,30,100);
    uy_cvt = m_cvt*ux_cvt;
    
    ux_cvt_line_1 = ux_cvt;
    uy_cvt_line_1 = uy_cvt + cd_cvt;
    ux_cvt_line_2 = ux_cvt;
    uy_cvt_line_2 = uy_cvt - cd_cvt;
else
    n_cvt = vxa0/vya0;
    cd_cvt = cd_cvt/vya0;
    uy_cvt = linspace(-30,30,100);
    ux_cvt = n_cvt*uy_cvt;
    
    ux_cvt_line_1 = ux_cvt + cd_cvt;
    uy_cvt_line_1 = uy_cvt;
    
    ux_cvt_line_2 = ux_cvt - cd_cvt;
    uy_cvt_line_2 = uy_cvt;
end

plot(ux_cvt_line_1,uy_cvt_line_1,'b','LineWidth',2);
plot(ux_cvt_line_2,uy_cvt_line_2,'b','LineWidth',2);

plot(all_samples_in_global_frame(1,:),all_samples_in_global_frame(2,:),'k.');
quiver(0,0,uxy_opt_global(1),uxy_opt_global(2),'r','LineWidth',1,'MaxHeadSize',1);

% draw control input acceleration magnitude contraint box
ru = [ux_max_org; uy_max_org];
rl = [ux_max_org; uy_min_org];
lu = [ux_min_org; uy_max_org];
ll = [ux_min_org; uy_min_org];

ru = dcm_from_body_to_global*ru;
rl = dcm_from_body_to_global*rl;
lu = dcm_from_body_to_global*lu;
ll = dcm_from_body_to_global*ll;

rec = [ru(:)' rl(:)';
       ru(:)' lu(:)';
       ll(:)' rl(:)';
       ll(:)' lu(:)'];
   
x_rec=[rec(:,1) rec(:,3)];
y_rec=[rec(:,2) rec(:,4)];

% draw later for legend
plot(x_rec',y_rec','g-.','LineWidth',2);

legend('cost function','$v_{\rm min}$','$v_{\rm max}$','uav velocity direction', ...
    'target location direction','curvature constraints','curvature constraints', ...
    'control input samples', 'optimal control input', ...
    'control input acceleration bound', ...
    'Location','northeastoutside','Interpreter','latex');

axis([ux_all(1) ux_all(end) uy_all(1) uy_all(end)]);

% % check the second derivative if it is convex
% dJ2=[eval(diff(dJdux0,ux0)) eval(diff(dJdux0,uy0)); 
%      eval(diff(dJduy0,ux0)) eval(diff(dJduy0,uy0))];
% 
% ux_all = linspace(-20,20,20);
% uy_all = linspace(-20,20,19); 
% eig_dJ2 = zeros(length(ux_all),length(uy_all));
% for idx=1:length(ux_all)
%     idx
%     for jdx=1:length(uy_all)
%         eig_dJ2(idx,jdx) = min(eig(eval(subs(dJ2,[ux0 uy0],[ux_all(idx) uy_all(jdx)]))));
%     end
% end
% figure; clf; contourf(ux_all,uy_all, eig_dJ2'); colorbar;


