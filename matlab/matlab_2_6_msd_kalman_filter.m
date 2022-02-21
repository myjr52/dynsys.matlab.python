clear;

%% simulation setting
m_mass = 1.0; %[kg]
k_spring = 0.5; %[N/m]
c_damper = 0.01; %[N/(m/s)]
msd_const = [m_mass k_spring c_damper];

init_pos = 0.0; %[m]
init_vel = 0.0; % [m/s]

init_time = 0; %[s]
final_time = 59; %[s]

Delta_t = 0.01; %[s]

time_interval = [init_time final_time];

%% process noise 
num_w = floor((final_time-init_time)/Delta_t)+1;
sigma_beta = sqrt(0.5); 
sigma_w =sigma_beta/sqrt(Delta_t);
wk_noise = sigma_w*(randn(num_w,1));

%% sensor
sigma_v = 0.75;%[0.75 1.5];
sample_step = 20;
current_step = 0;
Delta_tk_KF = Delta_t*sample_step;
zk_all = 1;%[0 0];

%% kalman filter
A_KF = [1 Delta_tk_KF; -(k_spring*Delta_tk_KF)/m_mass 1-(c_damper*Delta_tk_KF)/m_mass ];
H_KF = [1 0];%eye(2);

% kalman filter initialize
Q_KF = [0 0; 0 (Delta_tk_KF*sigma_w)^2];
R_KF = sigma_v^2;%[sigma_v(1)^2 0; 0 sigma_v(2)^2];
x_hat_plus = [0 0]';
P_plus_KF = 0.1*eye(2);
x_hat_all = x_hat_plus(:)';

%% Data 
x0 = [init_pos init_vel];
t0 = init_time;
tf = t0 + Delta_t;

tout_all = zeros(num_w,1);
xout_all = zeros(num_w,2);

tout_all(1) = t0;
xout_all(1,:) = x0;
xtrue_sample = x0;
P_all = P_plus_KF(:)';

%% simulation
for idx=2:num_w
    
    wk = wk_noise(idx);
    
    % RK45
    [tout,xout] = ode45( ...
        @(time,state)msd_noisy(time,state,wk,msd_const),...
        [t0 tf],x0);
    
    tout_all(idx) = tout(end);
    xout_all(idx,:) = xout(end,:);
    
    x0 = xout(end,:);
    
    % Kalman Filter
    if current_step >= sample_step
        
        % Prediction
        x_hat_minus = A_KF*x_hat_plus;
        P_minus_KF = A_KF*P_plus_KF*A_KF' + Q_KF;
        
        % Sensor measurement
        z_k = x0(1) + sigma_v*randn(1);%x0(:) + sigma_v(:).*randn(size(x0(:)));
        current_step = 1;
        zk_all = [ zk_all; z_k(:)'];
        
        % Update
        K_KF = P_minus_KF*H_KF'/(H_KF*P_plus_KF*H_KF'+R_KF);
        x_hat_plus = x_hat_minus + K_KF*(z_k - H_KF*x_hat_minus);
        P_plus_KF = (eye(2)-K_KF*H_KF)*P_minus_KF;
        
        % collect KF data
        x_hat_all = [x_hat_all; x_hat_plus(:)'];
        P_all = [P_all; P_plus_KF(:)'];
        
        % sample the true state to compare with the estimated
        xtrue_sample = [xtrue_sample; x0(:)'];
        
    else
        current_step = current_step + 1;
    end
    
    % time interval update
    t0 = tf;
    tf = t0 + Delta_t;
    
end

%% Error & 3-sigma bound calculation
x_error = xtrue_sample-x_hat_all;
position_3sigma = 3*sqrt(P_all(:,1));
velocity_3sigma = 3*sqrt(P_all(:,4));

%% plot
tout_measure = 0:Delta_tk_KF:(length(zk_all)-1)*Delta_tk_KF;
    
figure; clf;
subplot(211);
plot(tout_measure,zk_all(:,1),'k.');
hold on;
plot(tout_all,xout_all(:,1),'b--');
plot(tout_measure,x_hat_all(:,1),'r-');
legend('Measurement','True','Estimated')
axis([init_time final_time -10 10]);
set(gca,'FontSize',12);
ylabel('position [m]');
xlabel('time [s]');
subplot(212);
plot(tout_all,xout_all(:,2),'b--');
hold on;
plot(tout_measure,x_hat_all(:,2),'r-');
legend('True','Estimated');
axis([init_time final_time -10 10]);
set(gca,'FontSize',12);
ylabel('velocity [m/s]');
xlabel('time [s]');

figure; clf;
subplot(211);
plot(tout_measure,x_error(:,1))
hold on
plot(tout_measure,position_3sigma,'r--')
plot(tout_measure,-position_3sigma,'r--')
set(gca,'FontSize',12);
ylabel('[m]');
legend('Position Estimation Error','3\sigma bounds');
axis([init_time final_time -2 2]);
subplot(212);
plot(tout_measure,x_error(:,2))
hold on
plot(tout_measure,velocity_3sigma,'r--')
plot(tout_measure,-velocity_3sigma,'r--')
set(gca,'FontSize',12);
ylabel('[m/s]');
xlabel('time [s]');
axis([init_time final_time -10 10]);
legend('Velocity Estimation Error','3\sigma bounds');

%% system
function dxdt = msd_noisy(time,state,wk,msd_const)
    x1 = state(1);
    x2 = state(2);
    m = msd_const(1);
    k = msd_const(2);
    c = msd_const(3);
    
    dxdt = zeros(2,1);
    dxdt(1) = x2;
    dxdt(2) = -(k/m)*x1 - (c/m)*x2 + wk;
end