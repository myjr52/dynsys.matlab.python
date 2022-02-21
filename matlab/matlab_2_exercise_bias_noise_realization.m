clear;

%% Set initial values & change non-SI units into the SI Units

% number of stochastic process trial
N_realize = 10;

dt = 0.1; % [seconds]
time_init = 0;
time_final = 120;
time = time_init:dt:time_final;
N_sample = length(time);

% standard deviation of the bias, sigma_beta_xyz
sigma_beta_xyz = [0.01 0.01 0.02]; % [degrees/sqrt(s)]
sigma_beta_xyz = sigma_beta_xyz*(pi/180); % [rad/sqrt(s)]
sigma_eta_xyz = sigma_beta_xyz/sqrt(dt);

% store all beta history
beta_all = zeros(N_sample,3);

%% multiple realization of bias noise
figure(1);
for idx = 1:N_realize
    % initial beta(t)
    beta = (2*rand(1,3)-1)*0.03; % +/- 0.03[degrees/s]
    
    % from here all units are in SI except the plot parts
    beta = beta*(pi/180); % [radians/s]
    
    % main simulation loops
    for jdx=1:N_sample
        beta_all(jdx,:) = beta;
        
        eta_u = sigma_eta_xyz.*randn(1,3);
        dbeta = eta_u*dt;
        beta = beta + dbeta;
        
        beta_all(jdx,:) = beta(:)';
    end
    
    % plot all realization of beta in degrees/s
    figure(1);
    
    subplot(311);
    plot(time,beta_all(:,1)*180/pi,'r-'); 
    hold on;
    set(gca,'FontSize',14);
    ylabel('$\beta_x(t) [^\circ/s]$','Interpreter','latex');
    axis([time_init time_final -0.5 0.5]);
    
    subplot(312);
    plot(time,beta_all(:,2)*180/pi,'b-'); 
    hold on;
    set(gca,'FontSize',14);
    ylabel('$\beta_y(t) [^\circ/s]$','Interpreter','latex');
    axis([time_init time_final -0.5 0.5]);
    
    subplot(313);
    plot(time,beta_all(:,3)*180/pi,'g-'); 
    hold on;
    set(gca,'FontSize',14);
    ylabel('$\beta_z(t) [^\circ/s]$','Interpreter','latex');
    xlabel('time [s]','Interpreter','latex');
    axis([time_init time_final -0.5 0.5]);
end