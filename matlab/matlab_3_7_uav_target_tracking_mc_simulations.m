clear;

%% input variables

% initial uav position
xa0 = (2*rand(1)-1)*100; %[m]
ya0 = (2*rand(1)-1)*100; %[m]

% initial uav velocity
tha0 = rand(1)*2*pi; %[radian]
current_speed = 25; %[m/s]
vxa0 = current_speed*cos(tha0);
vya0 = current_speed*sin(tha0);

% control acceleration input magnitude constraints
ux_max = 10; % [m/s^2]
ux_min = -1; % [m/s^2]
uy_max = 2;  % [m/s^2]
uy_min = -2; % [m/s^2]

% uav minimum & maximum speed
v_max = 40; % [m/s]
v_min = 20; % [m/s]

% initial target position
xt0 = (2*rand(1)-1)*200;  %[m]
yt0 = (2*rand(1)-1)*200;  %[m]

% target maximum speed
w_max = 60*1e3/3600; %[m/s]

% uav flying path curvature constraint
r_min = 400; %[m]

% number of samples for the control search on the boundary
n_sample = 100;

% time interval for the cost approximation
Dt = 2; % [seconds]
%--------------------------------------------------------------


%% main part of monte-carlo simulations
%--------------------------------------------------------------
N_sim = 180;

% optimal control input error
opt_costh_err = zeros(1,N_sim-1);
opt_mag_err = zeros(1,N_sim-1);

% aircraft & target dynamics
Fa = eye(4) + [zeros(2) Dt*eye(2); zeros(2,4)];
Ga = [zeros(2); Dt*eye(2)];

Ft = eye(2);
Gt = Dt*eye(2);

state_uav = [xa0 ya0 vxa0 vya0]';
state_target = [xt0 yt0]';

uav_pos = [xa0 ya0];
target_pos = [xt0 yt0];

uxy_opt_global = [0 0];
w_zero_one = zeros(2,2);

figure(1); clf;
uav_line = animatedline('Color','b');
axis([-1 1 -1 1]*2000);
axis equal

%F(N_sim) = struct('cdata',[],'colormap',[]);

for idx_sim = 1:N_sim
    
    fprintf('%d/%d\n',idx_sim,N_sim);
    % save previous state
    state_minus = [xa0 ya0 vxa0 vya0 xt0 yt0 uxy_opt_global(:)'];
    
    % uav optimal input
    state_aircraft_tracking.aircraft = [xa0 ya0 vxa0 vya0 ux_min ux_max uy_min uy_max v_min v_max];
    state_aircraft_tracking.tracking = [xt0 yt0 w_max r_min n_sample Dt];
    [uxy_opt_global, uxy_opt_body, ux0_sample, uy0_sample, dcm_from_body_to_global] = ...
        uav_optimal_tracking_control(state_aircraft_tracking);
    
%     % use the trivial immediate target direction acceleration
%     % make the control input constraint polygon
%     polygon_points = [ux0_sample(:)'; uy0_sample(:)'];
%     polygon_center = mean(polygon_points,2);
%     pc_vector = polygon_points-polygon_center ;
%     th_pc = atan2(pc_vector(2,:),pc_vector(1,:));
%     [th_pc, idx_pc] = sort(th_pc);
%     polygon_points = polygon_points(:,idx_pc);
%     polygon_points = [polygon_points polygon_points(:,1)];
%     %plot(polygon_points(1,:), polygon_points(2,:), '.-r');
%     
%     target_direction = [xt0 yt0] - [xa0 ya0];
%     % make the length long enough to go outside the polygon
%     target_direction = (2*ux_max)*target_direction/norm(target_direction);
%     % fine the intersection
%     [uxi,uyj] = polyxpoly([0 target_direction(1)],[0 target_direction(2)], ...
%         polygon_points(1,:),polygon_points(2,:));
%     if isempty(uxi)
%         uxi=0; uyj=0;
%     end
%     uxy_opt_global = [uxi uyj]';
%     uxy_opt_body = dcm_from_body_to_global'*uxy_opt_global;
    
    % target maximum random direction input
    rand_th = 2*pi*rand(1);
    wx0 = w_max*cos(rand_th); 
    wy0 = w_max*sin(rand_th);
    
    % state propagation
    state_uav = Fa*state_uav + Ga*uxy_opt_global(:);
    state_target = Ft*state_target + Gt*[wx0; wy0];
    
    % reset state variables
    xa0 = state_uav(1);
    ya0 = state_uav(2);
    vxa0 = state_uav(3);
    vya0 = state_uav(4);    rand_th = 2*pi*rand(1);
    wx0 = w_max*cos(rand_th); 
    wy0 = w_max*sin(rand_th);
    xt0 = state_target(1);
    yt0 = state_target(2);
    
    % update the plot
    figure(1); 
    addpoints(uav_line,xa0,ya0);
    hold on;
    plot(xt0,yt0,'r.');
    hold off;
    %F(idx_sim) = getframe(gcf);
    pause(0.01);
    
    if mod(idx_sim,600)==0
        keyboard;
        figure(1); clf;
        uav_line = animatedline('Color','b');
        axis([-1 1 -1 1]*2000);
        axis equal
    end
    
    % save previous and current input
    w_zero_one = [w_zero_one(2,:); wx0 wy0];
    
    if idx_sim > 2
        xa0_minus = state_minus(1);
        ya0_minus = state_minus(2);
        vxa0_minus = state_minus(3);
        vya0_minus = state_minus(4);
        xt0_minus = state_minus(5);
        yt0_minus = state_minus(6);
        uxy0_minus = state_minus(7:8);
        
        state_aircraft_target.aircraft = [xa0_minus ya0_minus vxa0_minus vya0_minus Dt];
        state_aircraft_target.target = [xt0_minus yt0_minus w_zero_one(1,:) w_zero_one(2,:)];
        uxy_opt_true = calculate_original_cost(state_aircraft_target, ux0_sample, uy0_sample);
        
        uu_dot = uxy0_minus(:)'*uxy_opt_true(:);
        uxy0_minus_norm = norm(uxy0_minus);
        uxy_opt_true_norm = norm(uxy_opt_true);
        opt_costh_err(idx_sim-1) = uu_dot/(uxy0_minus_norm*uxy_opt_true_norm);
        opt_mag_err(idx_sim-1) = uxy0_minus_norm/uxy_opt_true_norm;
    end
    
end
    rand_th = 2*pi*rand(1);
    wx0 = w_max*cos(rand_th); 
    wy0 = w_max*sin(rand_th);
figure(2); clf;
subplot(211);
h1=histogram(opt_costh_err);
h1.Normalization='probability';
subplot(212);
h2=histogram(opt_mag_err);
h2.Normalization='probability';
%--------------------------------------------------------------


%% funciton: optimal target tracking control
%--------------------------------------------------------------
function [uxy_opt_global, uxy_opt_body, ux0_sample, uy0_sample, dcm_from_body_to_global] = ...
    uav_optimal_tracking_control(state_aircraft_tracking)

    aircraft = state_aircraft_tracking.aircraft;
    tracking = state_aircraft_tracking.tracking;

    xa0 = aircraft(1); 
    ya0 = aircraft(2);
    vxa0 = aircraft(3);
    vya0 = aircraft(4);
    ux_min = aircraft(5);
    ux_max = aircraft(6);
    uy_min = aircraft(7);
    uy_max = aircraft(8);
    v_min = aircraft(9); 
    v_max = aircraft(10);
    
    current_speed = sqrt(vxa0^2+vya0^2);
    
    xt0 = tracking(1);
    yt0 = tracking(2);
    w_max = tracking(3);
    r_min = tracking(4); 
    n_sample = tracking(5);
    Dt = tracking(6);
    
    J_cost_uxuy0_function = @(Dt,ux0,uy0,vxa0,vya0,w_max,xa0,xt0,ya0,yt0)xa0.*xt0.*-4.0-ya0.*yt0.*4.0+xa0.^2.*2.0+xt0.^2.*2.0+Dt.^4.*ux0.^2+Dt.^4.*uy0.^2+ya0.^2.*2.0+yt0.^2.*2.0+Dt.^2.*vxa0.^2.*5.0+Dt.^2.*vya0.^2.*5.0+Dt.^2.*w_max.^2.*2.0-abs(Dt).*abs(w_max).*sqrt(xa0.*xt0.*-8.0-ya0.*yt0.*8.0+xa0.^2.*4.0+xt0.^2.*4.0+Dt.^4.*ux0.^2+Dt.^4.*uy0.^2+ya0.^2.*4.0+yt0.^2.*4.0+Dt.^2.*vxa0.^2.*9.0+Dt.^2.*vya0.^2.*9.0+Dt.*vxa0.*xa0.*1.2e+1-Dt.*vxa0.*xt0.*1.2e+1+Dt.*vya0.*ya0.*1.2e+1-Dt.*vya0.*yt0.*1.2e+1+Dt.^3.*ux0.*vxa0.*6.0+Dt.^3.*uy0.*vya0.*6.0+Dt.^2.*ux0.*xa0.*4.0-Dt.^2.*ux0.*xt0.*4.0+Dt.^2.*uy0.*ya0.*4.0-Dt.^2.*uy0.*yt0.*4.0).*2.0+Dt.*vxa0.*xa0.*6.0-Dt.*vxa0.*xt0.*6.0+Dt.*vya0.*ya0.*6.0-Dt.*vya0.*yt0.*6.0+Dt.^3.*ux0.*vxa0.*4.0+Dt.^3.*uy0.*vya0.*4.0+Dt.^2.*ux0.*xa0.*2.0-Dt.^2.*ux0.*xt0.*2.0+Dt.^2.*uy0.*ya0.*2.0-Dt.^2.*uy0.*yt0.*2.0;
    
    % body and global coordinates transfrom dcm
    th_flight = atan2(vya0,vxa0);
    dcm_from_body_to_global = [cos(th_flight) -sin(th_flight); sin(th_flight) cos(th_flight)];

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

    % find the optimal solution along the boundary
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

    J_val = J_cost_uxuy0_function(Dt,ux0_sample,uy0_sample,vxa0,vya0,w_max,xa0,xt0,ya0,yt0);

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

    [in,~] = inpolygon(x_sample,y_sample,polygon_points(1,:),polygon_points(2,:));
    x_sample = x_sample(in);
    y_sample = y_sample(in);
    J_val_inside = J_cost_uxuy0_function(Dt,x_sample,y_sample,vxa0,vya0,w_max,xa0,xt0,ya0,yt0);
    J_val_inside = J_val_inside(J_val_inside<J_val_opt);

    if ~isempty(J_val_inside)
        [~,min_idx] = min(J_val_inside);
        J_cost_minimize=@(x)J_cost_uxuy0_function(Dt,x(1),x(2),vxa0,vya0,w_max,xa0,xt0,ya0,yt0);
        
        dJdux0_fun=@(Dt,ux0,uy0,vxa0,vya0,w_max,xa0,xt0,ya0,yt0)Dt.^4.*ux0.*2.0+Dt.^3.*vxa0.*4.0+Dt.^2.*xa0.*2.0-Dt.^2.*xt0.*2.0-Dt.^2.*abs(Dt).*abs(w_max).*(xa0.*2.0-xt0.*2.0+Dt.*vxa0.*3.0+Dt.^2.*ux0).*1.0./sqrt(xa0.*xt0.*-8.0-ya0.*yt0.*8.0+xa0.^2.*4.0+xt0.^2.*4.0+Dt.^4.*ux0.^2+Dt.^4.*uy0.^2+ya0.^2.*4.0+yt0.^2.*4.0+Dt.^2.*vxa0.^2.*9.0+Dt.^2.*vya0.^2.*9.0+Dt.*vxa0.*xa0.*1.2e+1-Dt.*vxa0.*xt0.*1.2e+1+Dt.*vya0.*ya0.*1.2e+1-Dt.*vya0.*yt0.*1.2e+1+Dt.^3.*ux0.*vxa0.*6.0+Dt.^3.*uy0.*vya0.*6.0+Dt.^2.*ux0.*xa0.*4.0-Dt.^2.*ux0.*xt0.*4.0+Dt.^2.*uy0.*ya0.*4.0-Dt.^2.*uy0.*yt0.*4.0).*2.0;
        dJduy0_fun=@(Dt,ux0,uy0,vxa0,vya0,w_max,xa0,xt0,ya0,yt0)Dt.^4.*uy0.*2.0+Dt.^3.*vya0.*4.0+Dt.^2.*ya0.*2.0-Dt.^2.*yt0.*2.0-Dt.^2.*abs(Dt).*abs(w_max).*(ya0.*2.0-yt0.*2.0+Dt.*vya0.*3.0+Dt.^2.*uy0).*1.0./sqrt(xa0.*xt0.*-8.0-ya0.*yt0.*8.0+xa0.^2.*4.0+xt0.^2.*4.0+Dt.^4.*ux0.^2+Dt.^4.*uy0.^2+ya0.^2.*4.0+yt0.^2.*4.0+Dt.^2.*vxa0.^2.*9.0+Dt.^2.*vya0.^2.*9.0+Dt.*vxa0.*xa0.*1.2e+1-Dt.*vxa0.*xt0.*1.2e+1+Dt.*vya0.*ya0.*1.2e+1-Dt.*vya0.*yt0.*1.2e+1+Dt.^3.*ux0.*vxa0.*6.0+Dt.^3.*uy0.*vya0.*6.0+Dt.^2.*ux0.*xa0.*4.0-Dt.^2.*ux0.*xt0.*4.0+Dt.^2.*uy0.*ya0.*4.0-Dt.^2.*uy0.*yt0.*4.0).*2.0;
        
        dJduxy=@(x)[dJdux0_fun(Dt,x(1),x(2),vxa0,vya0,w_max,xa0,xt0,ya0,yt0); 
                    dJduy0_fun(Dt,x(1),x(2),vxa0,vya0,w_max,xa0,xt0,ya0,yt0)];
        
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
   
        uxy_opt_global = u_xy_current(:);
        uxy_opt_body = dcm_from_body_to_global'*uxy_opt_global;

        fprintf('******optimal inside the constraints************\n');
        fprintf('(ux,uy)*  = (%5.4f, %5.4f)\n',uxy_opt_body(1),uxy_opt_body(2));
        fprintf('************************************************\n');
    end
    
end % of the function


%% function: original cost function
function u_opt_true = calculate_original_cost(state_aircraft_target, ux0_sample, uy0_sample)
    
    aircraft = state_aircraft_target.aircraft;
    target = state_aircraft_target.target;

    xa0 = aircraft(1); 
    ya0 = aircraft(2);
    vxa0 = aircraft(3);
    vya0 = aircraft(4);
    Dt = aircraft(5);
    
    xt0 = target(1);
    yt0 = target(2);
    wx0 = target(3);
    wy0 = target(4);
    wx1 = target(5);
    wy1 = target(6);
    
    % global solution wihtout the constraint
    alpha = -(2*Dt^3*wx0 - 4*Dt^3*vxa0 + 2*Dt^3*wx1 - 2*Dt^2*xa0 + 2*Dt^2*xt0)/Dt^4;
    beta = -(2*Dt^3*wy0 - 4*Dt^3*vya0 + 2*Dt^3*wy1 - 2*Dt^2*ya0 + 2*Dt^2*yt0)/Dt^4;
    ux_opt = -alpha/2;
    uy_opt = -beta/2;

    % check the global solution if it is within the constraint
    polygon_points = [ux0_sample(:)'; uy0_sample(:)']; 
    polygon_center = mean(polygon_points,2); 
    pc_vector = polygon_points-polygon_center ; 
    th_pc = atan2(pc_vector(2,:),pc_vector(1,:));
    [~, idx_pc] = sort(th_pc);
    polygon_points = polygon_points(:,idx_pc); 
    polygon_points = [polygon_points polygon_points(:,1)]; 
  
    [in,on] = inpolygon(ux_opt,uy_opt,polygon_points(1,:),polygon_points(2,:));
    if ~in && ~on
        J_val_func = @(Dt,ux0,uy0,vxa0,vya0,wx0,wx1,wy0,wy1,xa0,xt0,ya0,yt0)xa0.*xt0.*-4.0-ya0.*yt0.*4.0+xa0.^2.*2.0+xt0.^2.*2.0+Dt.^4.*ux0.^2+Dt.^4.*uy0.^2+ya0.^2.*2.0+yt0.^2.*2.0+Dt.^2.*vxa0.^2.*5.0+Dt.^2.*vya0.^2.*5.0+Dt.^2.*wx0.^2.*2.0+Dt.^2.*wx1.^2+Dt.^2.*wy0.^2.*2.0+Dt.^2.*wy1.^2+Dt.*vxa0.*xa0.*6.0-Dt.*vxa0.*xt0.*6.0-Dt.*wx0.*xa0.*4.0-Dt.*wx1.*xa0.*2.0+Dt.*vya0.*ya0.*6.0+Dt.*wx0.*xt0.*4.0+Dt.*wx1.*xt0.*2.0-Dt.*vya0.*yt0.*6.0-Dt.*wy0.*ya0.*4.0-Dt.*wy1.*ya0.*2.0+Dt.*wy0.*yt0.*4.0+Dt.*wy1.*yt0.*2.0+Dt.^3.*ux0.*vxa0.*4.0+Dt.^3.*uy0.*vya0.*4.0-Dt.^3.*ux0.*wx0.*2.0-Dt.^3.*ux0.*wx1.*2.0-Dt.^3.*uy0.*wy0.*2.0-Dt.^3.*uy0.*wy1.*2.0+Dt.^2.*ux0.*xa0.*2.0-Dt.^2.*ux0.*xt0.*2.0-Dt.^2.*vxa0.*wx0.*6.0-Dt.^2.*vxa0.*wx1.*4.0-Dt.^2.*vya0.*wy0.*6.0-Dt.^2.*vya0.*wy1.*4.0+Dt.^2.*uy0.*ya0.*2.0-Dt.^2.*uy0.*yt0.*2.0+Dt.^2.*wx0.*wx1.*2.0+Dt.^2.*wy0.*wy1.*2.0;
        J_cost_val=J_val_func(Dt,ux0_sample,uy0_sample,vxa0,vya0,wx0,wx1,wy0,wy1,xa0,xt0,ya0,yt0);
        [~,opt_idx]=min(J_cost_val);
        u_opt_true = [ux0_sample(opt_idx) uy0_sample(opt_idx)];
    else
        u_opt_true = [ux_opt uy_opt]';
    end
    
        
end % of the function
