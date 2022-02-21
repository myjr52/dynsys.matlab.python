clear;

%% simulation parameters
r_vehicle = [0 5]; % initial vehicle position (x,y) [m]
v_vehicle = [0 0]; % initial vehicle velocity (vx,vy) [m/s]

r_desired = [20 5]; % desired vehicle position (xdst, ydst) [m]
r_obs = [5 8; 10 5; 15 8; 15 2]; % obstacles (xobs,yobs) [m]
rho_o_i = 2.4; % obstacle radius [m]

ka = 0.5;
kr = 100;
c_damping = 5;
F_attractive_max = 10;

time_interval = [0 50]; % [s]
state_0 = [r_vehicle v_vehicle];

% Perturbation Force
dt = 0.1;
wk_mag = 0.01;
N_noise = floor(diff(time_interval)/dt);
wk_time = linspace(time_interval(1),time_interval(2),N_noise);
wk_noise = wk_mag*(2*rand(N_noise,2)-1);

sim_para = {r_desired, r_obs, rho_o_i, ka, kr, c_damping, F_attractive_max, wk_time, wk_noise};


%% simulation
ode_options = odeset('RelTol',1e-2,'AbsTol',1e-3, 'MaxStep', 0.1);
[tout,state_out] = ode45(@(time,state) drvdt_potential_field(time,state,sim_para), ...
    time_interval, state_0, ode_options);

rout = state_out(:,1:2);
vout = state_out(:,3:4);

%% plot
th = 0:0.01:2*pi;
xcircle = rho_o_i*cos(th);
ycircle = rho_o_i*sin(th);

figure; clf;
hold on;
for idx = 1:size(r_obs,1)
   plot(r_obs(idx,1)+xcircle, r_obs(idx,2)+ycircle,'r');
end
plot(rout(:,1),rout(:,2),'b','LineWidth',2);
set(gca,'FontSize',14);
ylabel('y [m]');
xlabel('x [m]');
axis equal;

%% functions
function dstate_dt = drvdt_potential_field(time,state,sim_para)

    % states
    x_vehicle = state(1);
    y_vehicle = state(2);
    v_current = state(3:4);
    
    % simulation setting
    xy_dst = sim_para{1};
    xy_obs = sim_para{2};
    rho_o_i = sim_para{3};
    ka = sim_para{4};
    kr = sim_para{5};
    c_damping = sim_para{6};
    Famax = sim_para{7}*ones(2,1);
    wk_time = sim_para{8};
    wk_noise = sim_para{9};
    
    num_obs = size(xy_obs,1);
    
    % desired position
    x_dst = xy_dst(1);
    y_dst = xy_dst(2);
    
    % attaractive & damping force
    Fa = -ka*[(x_vehicle-x_dst); (y_vehicle-y_dst)];
    Fa = sign(Fa).*min([abs(Fa(:)) Famax(:)],[],2);
    Fd = -c_damping*v_current(:);
    
    % repulsive force
    Fr = [0; 0]; 
    for idx=1:num_obs
        
        x_ost = xy_obs(idx,1);
        y_ost = xy_obs(idx,2);
        rho_r_i = sqrt((x_vehicle-x_ost)^2+(y_vehicle-y_ost)^2);
        if rho_r_i > rho_o_i
            Frx_idx = 0;
            Fry_idx = 0;
        else
            Frx_idx = kr*(x_vehicle-x_ost)/(rho_r_i^3);
            Fry_idx = kr*(y_vehicle-y_ost)/(rho_r_i^3);
        end
        
        Fr(1) = Fr(1) + Frx_idx;
        Fr(2) = Fr(2) + Fry_idx;
    end
    
    wk = interp1(wk_time,wk_noise,time);
    
    F_sum = Fa(:) + Fr(:) + Fd(:) + wk(:);
    
    drdt = v_current(:);
    dvdt = F_sum;
    dstate_dt = [drdt; dvdt];

end