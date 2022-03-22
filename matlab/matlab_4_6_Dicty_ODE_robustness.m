clear;

FitnessFunction = @(delta)Dicty_x1_square_integral(delta);
lb = -ones(1,14);
ub = ones(1,14);
opt_opts = optimoptions('ga','Display','iter');

[delta_worst,fval] = ga(FitnessFunction,14,[],[],[],[],lb,ub,[],opt_opts);

%% Cost function to be minimized for robustness analysis
function J_cost = Dicty_x1_square_integral(delta)

    ki_para_org = [2.0; 0.9; 2.5; 1.5; 0.6; 0.8; 1.0; 1.3; 0.3; 0.8; 0.7; 4.9; 23.0; 4.5];
    p_delta = 2; % [percents]
    ki_para=ki_para_org.*(1+(p_delta/100)*delta(:));

    x0 = rand(7,1);
    dt = 0.1;  % [min]
    t0 = 600;  % [min]
    tf = 1200; % [min]
    
    time_interval = 0:dt:tf; % [min]

    [~,xout] = ode45(@(time,state) Dicty_cAMP(time,state,ki_para), time_interval, x0);

    N_t0 = floor(t0/dt);
    
    ACA = xout(N_t0:end,1);
    PKA = xout(N_t0:end,2);
    CAR1 = xout(N_t0:end,7);
    J_cost = sum((ki_para(1)*CAR1 - ki_para(2)*(ACA.*PKA)).^2)*dt*0.5;
    
end


%% Dicty cAMP network
function dxdt = Dicty_cAMP(time,state,ki_para)
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