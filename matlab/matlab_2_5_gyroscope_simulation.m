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

%% Set initial values & change non-SI units into the SI Units
dt = 0.05; % [seconds]
time_init = 0;
time_final = 120;
time = time_init:dt:time_final;
N_sample = length(time);

% standard deviation of the bias, sigma_beta_xyz
sigma_beta_xyz = [0.05 0.04 0.06]; % [(degrees/s)/sqrt(s)]
sigma_beta_xyz = sigma_beta_xyz*(pi/180); % [(rad/s)/sqrt(s)]
sigma_eta_xyz = sigma_beta_xyz/sqrt(dt);

% standard devitation of the white noise, sigma_v
sigma_v = 0.01; %[degrees/s]
sigma_v = sigma_v*(pi/180); %[rad/s]

% initial beta(t)
beta = (2*rand(3,1)-1)*0.03; % +/- 0.03[degrees/s]
beta = beta*(pi/180); % [radians/s]

% prepare the data store
w_all = zeros(N_sample,3);
w_measure_all = zeros(N_sample,3);

%% main simulation loop
for idx=1:N_sample
   
    time_c = time(idx);
    w_true(1,1) = 0.1*sin(2*pi*0.005*time_c); % [rad/s]
    w_true(2,1) = 0.05*cos(2*pi*0.01*time_c + 0.2); %[rad/s]
    w_true(3,1) = 0.02; %[rad/s]
    
    % beta(t)
    eta_u = sigma_eta_xyz(:).*randn(3,1);
    dbeta = eta_u*dt;
    beta = beta + dbeta;
    
    % eta_v(t)
    eta_v = sigma_v*randn(3,1);
    
    % w_tilde
    w_measurement = w_true + beta + eta_v;
    
    % store history
    w_all(idx,:) = w_true(:)';
    w_measure_all(idx,:) = w_measurement(:)';
    
end
    
% plot in degrees/s
figure;
plot(time,w_all*(180/pi));
hold on;
plot(time,w_measure_all*(180/pi),'--');
set(gca,'FontSize',14);
ylabel('$[^\circ/s]$','Interpreter','latex');
xlabel('time [s]','Interpreter','latex');
legend('$\omega_x$','$\omega_y$','$\omega_z$', ...
    '$\tilde{\omega}_x$','$\tilde{\omega}_y$','$\tilde{\omega}_z$', ...
    'Interpreter','latex','Location','SouthWest');
