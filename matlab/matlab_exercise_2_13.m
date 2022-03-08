clear;

%% Set initial values & change non-SI units into the SI Units
dt = 0.05; % [seconds]
time_init = 0;
time_final = 120;
time = time_init:dt:time_final;
N_sample = length(time);

%% gyroscope
% standard deviation of the bias, sigma_beta_xyz
sigma_beta = 0.0005; % [(degrees/s)/sqrt(s)]
sigma_u = sigma_beta*(pi/180); % [(rad/s)/sqrt(s)]
sigma_eta = sigma_u/sqrt(dt);

% standard devitation of the white noise, sigma_v
sigma_v = 0.0001; %[degrees/s]
sigma_v = sigma_v*(pi/180); %[rad/s]

% initial beta(t)
beta = (2*rand(3,1)-1)*0.03; % +/- 0.03[degrees/s]
beta = beta*(pi/180); % [radians/s]

%% data store
% prepare the data store
w_all = zeros(N_sample,3);
w_measure_all = zeros(N_sample,3);

% instead of calculating the exact size
% of the following matrices, use varying matrices with increasing
% time, which might not be significant but simpler to implement
w_gyr_all = [];
w_hat_all = [];
w_tr_all = [];

q_tr_all = [];
q_hat_all = [];
time_all = [];
pcov_all = [];

q_current = [0 0 0 1]'; % current quaternion

%% star sensor
% star sensor reference star vectors
r1R = [-0.6794 -0.3237 -0.6586]'; r1R = r1R/norm(r1R);
r2R = [-0.7296  0.5858  0.3528]'; r2R = r2R/norm(r2R);
r3R = [-0.2718  0.6690 -0.6918]'; r3R = r3R/norm(r3R);
r4R = [-0.2062 -0.3986  0.8936]'; r4R = r4R/norm(r4R); 
r5R = [ 0.6858 -0.7274 -0.0238]'; r5R = r5R/norm(r5R);

rR_star_all = [r1R r2R r3R r4R r5R];
num_star = size(rR_star_all,2);
sigma_star = 87.2665/3*1e-6; 
r_star = sigma_star^2*eye(num_star*3);

%% Kalman filter
n_dt_KF =2;   % n_dt_KF * (gyro sensor simulation time interval) 
dt_KF = n_dt_KF*dt; % Kalman filter sampling time interval
bias_estimate_current = [0 0 0];
q_estimate_current = [0 0 0 1] + 0.01*randn(1,4); 
q_estimate_current = q_estimate_current/norm(q_estimate_current);
x0 = [q_estimate_current bias_estimate_current];
p_current = 0.001*eye(6);
w_hat = [0 0 0];

%% main simulation loop
for idx=1:N_sample
   
    time_c = time(idx);
    w_true = angular_velocity_true(time_c);
    w_true = w_true(:);
    
    % beta(t)
    eta_u = sigma_eta*randn(3,1);
    dbeta = eta_u*dt;
    beta = beta + dbeta;
    
    % eta_v(t)
    eta_v = sigma_v*randn(3,1);
    
    % w_tilde
    w_measurement = w_true + beta + eta_v;
    
    % run Kalman filter
    if rem(idx,n_dt_KF)==1
        
        % star_sensor measurement
        dcm_BR = q2dcm(q_current);
        rB_star_all = dcm_BR*rR_star_all;
        rB_star_measure = rB_star_all+sigma_star*randn(3,num_star);
        rB_star_measure = rB_star_measure./(ones(3,1)*sqrt(sum(rB_star_measure.^2)));
        
        % kalman filter
        [x_hat_1, P1] = kalman_filter_attitude(x0, p_current, dt_KF, ...
            rR_star_all, rB_star_measure, w_measurement, sigma_v, sigma_u, sigma_star);
        
        x0 = x_hat_1;
        p_current = P1;
        
        q_estimate_current    = x0(1:4); q_estimate_current = q_estimate_current(:);
        bias_estimate_current = x0(5:7); bias_estimate_current = bias_estimate_current(:); 
        
        w_hat = w_measurement(:) - bias_estimate_current(:);
        
        % store data to plot
        % instead of calculating the exact size
        % of the following matrices, use varying matrices with increasing
        % time, which might not be significant but simpler to implement
        time_all = [time_all; time_c];
        w_gyr_all = [w_gyr_all; w_measurement(:)'];
        w_hat_all = [w_hat_all; w_hat(:)'];
        w_tr_all = [w_tr_all; w_true(:)'];
        
        q_tr_all = [q_tr_all; q_current(:)'];
        q_hat_all = [q_hat_all; q_estimate_current(:)'];
        pcov_all = [pcov_all; diag(P1)'];
    end
    
    % integrate true dqdt to obtain true q(t): time_c -> time_c + dt
    if idx < N_sample
       [tout,qout] = ode45(@(time,state) dqdt_attitude_kinematics(time,state), ...
                           [time_c time(idx+1)], q_current);
       q_current = qout(end,:);
    end
    
   
    
end

% plot results

figure;
subplot(311);
plot(time_all,q_tr_all(:,1)-q_hat_all(:,1));
axis([0 time_final -4e-5 4e-5]);
hold on;
plot(time_all,3*sqrt(pcov_all(:,1)),'r--');
plot(time_all,-3*sqrt(pcov_all(:,1)),'r--');
ylabel('\delta q_1');
legend('error','3\sigma bound');
subplot(312);
plot(time_all,q_tr_all(:,2)-q_hat_all(:,2));
axis([0 time_final -4e-5 4e-5]);
hold on;
plot(time_all,3*sqrt(pcov_all(:,2)),'r--');
plot(time_all,-3*sqrt(pcov_all(:,2)),'r--');
ylabel('\delta q_2');
legend('error','3\sigma bound');
subplot(313);
plot(time_all,q_tr_all(:,3)-q_hat_all(:,3));
axis([0 time_final -4e-5 4e-5]);
hold on;
plot(time_all,3*sqrt(pcov_all(:,3)),'r--');
plot(time_all,-3*sqrt(pcov_all(:,3)),'r--');
ylabel('\delta q_3');
legend('error','3\sigma bound');
xlabel('time [s]');

figure;
subplot(311);
plot(time_all,w_tr_all(:,1)-w_hat_all(:,1));
axis([0 time_final -4e-5 4e-5]);
hold on;
plot(time_all,3*sqrt(pcov_all(:,4)),'r--');
plot(time_all,-3*sqrt(pcov_all(:,4)),'r--');
ylabel('\delta \omega_1 [rad/s]');
legend('error','3\sigma bound');
subplot(312);
plot(time_all,w_tr_all(:,2)-w_hat_all(:,2));
axis([0 time_final -4e-5 4e-5]);
hold on;
plot(time_all,3*sqrt(pcov_all(:,5)),'r--');
plot(time_all,-3*sqrt(pcov_all(:,5)),'r--');
ylabel('\delta \omega_2 [rad/s]');
legend('error','3\sigma bound');
subplot(313);
plot(time_all,w_tr_all(:,3)-w_hat_all(:,3));
axis([0 time_final -4e-5 4e-5]);
hold on;
plot(time_all,3*sqrt(pcov_all(:,6)),'r--');
plot(time_all,-3*sqrt(pcov_all(:,6)),'r--');
ylabel('\delta \omega_3 [rad/s]');
legend('error','3\sigma bound');
xlabel('time [s]');


% figure;
% subplot(411);
% plot(time_all,q_tr_all(:,1),'k');
% hold on;
% plot(time_all,q_hat_all(:,1),'g-.');
% ylabel('q_1');
% legend('true','estimate');
% subplot(412);
% plot(time_all,q_tr_all(:,2),'k');
% hold on;
% plot(time_all,q_hat_all(:,2),'g-.');
% ylabel('q_2');
% legend('true','estimate');
% subplot(413);
% plot(time_all,q_tr_all(:,3),'k');
% hold on;
% plot(time_all,q_hat_all(:,3),'g-.');
% ylabel('q_3');
% legend('true','estimate');
% subplot(414);
% plot(time_all,q_tr_all(:,4),'k');
% hold on;
% plot(time_all,q_hat_all(:,4),'g-.');
% ylabel('q_4');
% xlabel('time [s]');
% legend('true','estimate');
% 
% figure;
% subplot(311);
% plot(time_all,w_tr_all(:,1),'k');
% hold on;
% plot(time_all,w_gyr_all(:,1),'r--');
% plot(time_all,w_hat_all(:,1),'g-.');
% ylabel('\omega_x [rad/s]');
% legend('true','measurement','estimate');
% subplot(312);
% plot(time_all,w_tr_all(:,2),'k');
% hold on;
% plot(time_all,w_gyr_all(:,2),'r--');
% plot(time_all,w_hat_all(:,2),'g-.');
% ylabel('\omega_y [rad/s]');
% legend('true','measurement','estimate');
% subplot(313);
% plot(time_all,w_tr_all(:,3),'k');
% hold on;
% plot(time_all,w_gyr_all(:,3),'r--');
% plot(time_all,w_hat_all(:,3),'g-.');
% xlabel('time [s]');
% ylabel('\omega_z [rad/s]');
% legend('true','measurement','estimate');

%--------------------------------------------------------------------------
function w_true = angular_velocity_true(time)
    w_true(1) = 0.01*sin(2*pi*0.005*time); % [rad/s]
    w_true(2) = 0.05*cos(2*pi*0.001*time + 0.2); %[rad/s]
    w_true(3) = 0.02; %[rad/s]
end

function dcm = q2dcm(quat) 
% be careful matlab aerospace toolbox has quat2dcm function

q13 =quat(1:3); q13 = q13(:);
q4 = quat(4);

q13x = [ 0          -q13(3)       q13(2);
         q13(3)      0           -q13(1);
        -q13(2)      q13(1)       0];

dcm = (q4^2-q13'*q13)*eye(3) + 2*(q13*q13') - 2*q4*q13x;
end

function dqdt = dqdt_attitude_kinematics(time,state)
    q_true = state(:);
    w_true = angular_velocity_true(time);
    w_true = w_true(:);
    
    wx = [  0           -w_true(3)  w_true(2);
            w_true(3)   0           -w_true(1);
            -w_true(2)  w_true(1)   0];
        
    Omega = [   -wx         w_true;
                -w_true'    0];
            
    dqdt = 0.5*Omega*q_true;
end


function [x_hat_1, P1] = kalman_filter_attitude(x_hat_0, P0, dt_KF, rR_star_all, rB_star_measure, w_measure, sgm_v, sgm_u, sgm_star)

num_KF_state = 6;  %[dq1 dq2 dq3 b1 b2 b3]

q_est = x_hat_0(1:4); q_est = q_est(:);
b_est = x_hat_0(5:7); b_est = b_est(:);

w_hat = w_measure(:) - b_est;
w_hat_mag = sqrt(sum(w_hat.^2));
w_hatx = [   0           -w_hat(3)         w_hat(2);
             w_hat(3)     0               -w_hat(1);
            -w_hat(2)     w_hat(1)         0];

% propagate   
if w_hat_mag > 1e-12
   dtheta_k = w_hat_mag*dt_KF;
   cos_th = cos(dtheta_k/2);
   sin_th_over_w = sin(dtheta_k/2)/w_hat_mag;
   
   q_Phi = [  cos_th*eye(3)-sin_th_over_w*w_hatx      sin_th_over_w*w_hat;
            -sin_th_over_w*w_hat'                   cos_th];
        
   cos_th = cos(w_hat_mag*dt_KF);
   sin_th = sin(w_hat_mag*dt_KF);
   Phi_1 = eye(3) - w_hatx*sin_th/w_hat_mag + (w_hatx*w_hatx)*((1-cos_th)/w_hat_mag^2);
   Phi_2 = -eye(3)*dt_KF + w_hatx*((1-cos_th)/w_hat_mag^2) - (w_hatx*w_hatx)*((w_hat_mag*dt_KF-sin_th)/w_hat_mag^3);
else
   q_Phi = [  eye(3)-(dt_KF/2)*w_hatx          (dt_KF/2)*w_hat;
            -(dt_KF/2)*w_hat'                  1];
        
   Phi_1 = eye(3) - w_hatx*dt_KF;
   Phi_2 = -eye(3)*dt_KF;
end

q_est_minus = q_Phi*q_est; q_est_minus = q_est_minus/norm(q_est_minus);
b_est_minus = b_est;
dcm_BR_minus = q2dcm(q_est_minus);

Q = [(sgm_v^2*dt_KF+(dt_KF^3/3)*sgm_u^2)*eye(3)                 -(dt_KF^2/2)*sgm_u^2*eye(3)-(dt_KF^3/6)*sgm_u^2*w_hatx;
     -(dt_KF^2/2)*sgm_u^2*eye(3)-(dt_KF^3/6)*sgm_u^2*w_hatx     sgm_u^2*dt_KF*eye(3)];

Phi = [Phi_1 Phi_2; zeros(3) eye(3)];
P1 = Phi*P0*Phi' + Q;

% update
num_star = size(rB_star_measure,2);
rB_star_hat = dcm_BR_minus*rR_star_all;
H_k = zeros(3*num_star,6);
R = sgm_star^2*eye(num_star*3);
for xdx = 1:num_star
    vec = rB_star_hat(:,xdx);
    vec_x = [ 0          -vec(3)      vec(2); 
              vec(3)      0          -vec(1); 
             -vec(2)      vec(1)      0];
    st_idx = 3*(xdx-1)+1;
    H_k(st_idx:st_idx+2,:) = [vec_x zeros(3)];
end

K_k = P1*H_k'/(H_k*P1*H_k' + R);
P1 = (eye(num_KF_state)-K_k*H_k)*P1;
delta_x = K_k*(rB_star_measure(:)-rB_star_hat(:));

% quaternion & bias update
dq_13 = 2*delta_x(1:3); dq_13 = dq_13(:);

q = q_est_minus(:);
qx = [0 -q(3) q(2);
      q(3) 0 -q(1);
      -q(2) q(1) 0];
quat_update_matrix = [q(4)*eye(3)+qx; -q(1:3)'];

q_hat_plus = q_est_minus + quat_update_matrix*dq_13; q_hat_plus = q_hat_plus/norm(q_hat_plus);
b_hat_plus = b_est_minus(:) + delta_x(4:6);

x_hat_1 = [q_hat_plus(:); b_hat_plus(:)];

end