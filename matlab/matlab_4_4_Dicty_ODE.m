clear;

init_time = 0; % [min]
final_time = 1800.0; % [min]
time_interval = [init_time final_time];

ki_para_org = [2.0; 0.9; 2.5; 1.5; 0.6; 0.8; 1.0; 1.3; 0.3; 0.8; 0.7; 4.9; 23.0; 4.5];


% start with a random initial condition and simiulate long enough to
% reach the stable oscillation
x0 = rand(7,1);
[~,xout] = ode45(@(time,state) Dicty_cAMP(time,state,ki_para_org), time_interval, x0);

% robustness
%delta_worst = [-1 -1 1 1 -1 1 1 -1 1 1 -1 1 -1 1]';
%delta_worst = sign(2*rand(14,1)-1);
delta_worst = [1    -1     1     1    -1     1    -1     1    -1     1    -1     1    -1     1]';
p_delta = 20;
ki_para=ki_para_org.*(1+(p_delta/100)*delta_worst);

% start with the initial condtion on the oscillation trajectory
x0 = xout(end,:);
time_interval = [0 60]; % [min]
[tout,xout] = ode45(@(time,state) Dicty_cAMP(time,state,ki_para), time_interval, x0);

figure;
plot(tout,xout(:,5),'-');
hold on;
plot(tout,xout(:,7),'-.');
set(gca,'FontSize',14);
xlabel('time [min]');
ylabel('[\muM]');
legend('i-cAMP','CAR1');


%% Dicty cAMP network
function dxdt = Dicty_cAMP(~,state,ki_para)
    ACA   = state(1);
    PKA   = state(2);
    ERK2  = state(3);
    REGA  = state(4);
    icAMP = state(5);
    ecAMP = state(6);
    CAR1  = state(7);
    
    k1 = ki_para(1);
    k2 = ki_para(2);
    k3 = ki_para(3);
    k4 = ki_para(4);
    k5 = ki_para(5);
    k6 = ki_para(6);
    k7 = ki_para(7);
    k8 = ki_para(8);
    k9 = ki_para(9);
    k10 = ki_para(10);
    k11 = ki_para(11);
    k12 = ki_para(12);
    k13 = ki_para(13);
    k14 = ki_para(14);
    
    dACA_dt = k1*CAR1 - k2*ACA*PKA;
    dPKA_dt = k3*icAMP - k4*PKA;
    dERK2_dt = k5*CAR1 - k6*PKA*ERK2;
    dREGA_dt = k7 - k8*ERK2*REGA;
    dicAMP_dt = k9*ACA - k10*REGA*icAMP;
    decAMP_dt = k11*ACA - k12*ecAMP;
    dCAR1_dt = k13*ecAMP - k14*CAR1;
    
    dxdt = [dACA_dt; dPKA_dt; dERK2_dt; dREGA_dt; dicAMP_dt; decAMP_dt; dCAR1_dt];
end