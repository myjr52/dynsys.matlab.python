
# parameters
kP = 50
kI = 5e-6
gamma_G = 8e-4
ks2 = 4e-4

#----------------- Change Here------------------------------
Kdeg = 1e-3   #1e-3    # PI-Deg gain
X3 = 1 
KF = 3  # 0.62; 3      # Pd scaling

eta = 1e-4             # integrator self-annihilation  

Pd = 1                     # Desired P
time_for_step_change = 8   # [hours]
#-----------------------------------------------------------

kon = 5e-5
kcat = 1.6*2
kdg = 8e-8
ks1 = 3

E_total = 0.2

# simulation time values
time_current = 0       # initial time 
time_final   = 3600*16 # final time [min]
tspan = [time_current, time_final]


%% simulation
para = [kP kI gamma_G ks2 kon kcat kdg ks1 Pd KF Kdeg X3 eta time_for_step_change];
ode_option = odeset('RelTol',1e-3,'AbsTol',1e-6);
state_t0    = 0.1*ones(1,8); 
state_t0(3) = E_total - state_t0(1);
state_t0(4) = 1e-3; 
state_t0(5) = 1e-3;
sol = ode15s(@(time, state)ES_PI_Half_Subtraction(time, state, para),tspan, ...
    state_t0, ode_option);

%----------------- Change Here------------------------------
figure; 
clf;
subplot(211);
time_hr = sol.x/3600; % [hour]
P_history = sol.y(2,:); 

%plot(time_hr,Pd*ones(size(time_hr)),'r--');
plot([0:time_for_step_change time_for_step_change:time_final/3600], ...
     [ones(1,9)*Pd ones(1,9)*Pd/2],'r--');

hold on;
plot(time_hr,P_history,'b');
set(gca,'FontSize',14);
ylabel('[P(t)] [a.u.]');
legend('desired [P]', 'achieved [P]');
axis([0 time_final/3600 0 1.2]);
subplot(212);
plot(sol.x/3600,sol.y(6,:));
set(gca,'FontSize',14);
xlabel('time [hour]');
ylabel('[S] [a.u.]');
axis([0 time_final/3600, 0 120]);

%% E-S PI Control Half Subtraction
function dxdt = ES_PI_Half_Subtraction(time,state,ki_para)
    E   = state(1);
    P   = state(2);
    ES  = state(3);
    X1  = state(4);
    X2  = state(5);
    S   = state(6);
    DP  = state(7);
    Xs  = state(8);
    
    kP = ki_para(1);
    kI = ki_para(2);
    gamma_G = ki_para(3);
    ks2 = ki_para(4);
    kon = ki_para(5);
    kcat = ki_para(6);
    kdg = ki_para(7);
    ks1 = ki_para(8);
    Pd = ki_para(9);
    KF = ki_para(10);
    Kdeg = ki_para(11);
    X3 = ki_para(12);
    eta = ki_para(13);
    time_for_step_change = ki_para(14);
    
    if time > 3600*time_for_step_change
        Pd = Pd/2;
    end
    
    Pd = Pd*KF; % Pd scaling
    
    dE_dt = -kon*E*S + kcat*ES;
    dP_dt = kcat*ES - kdg*P - Kdeg*X3*P;
    dES_dt = kon*E*S - kcat*ES;
    dX1_dt = kI*DP - eta*X1;
    dX2_dt = -gamma_G*X2 + gamma_G*kP*DP;
    dS_dt = ks2*X1 + ks2*X2 - ks2*S;
    dDP_dt = ks1*Pd - ks1*DP*Xs - ks1*DP;
    dXs_dt = -ks1*DP*Xs + ks1*P;
    
    dxdt = [dE_dt; dP_dt; dES_dt; dX1_dt; dX2_dt; dS_dt; dDP_dt; dXs_dt];
    
    
end