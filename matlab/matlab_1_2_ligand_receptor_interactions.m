clear;

init_receptor = 0.01; % [nM]
init_ligand = 0.0415; %[nM]
init_complex = 0.0; %[kg]

init_time = 0; % [min]
final_time = 180.0; % [min]
time_interval = [init_time final_time];

kon = 0.0972; % [1/(min nM)] 
koff = 0.24; % [1/min]
kt = 0.02; %[1/min]
ke = 0.15; % [1/min]

ft = 0.0; % [nM/min]
QR = 0.0166; % [nM/min]
R_max = 0.415; %[nM]

sim_para = [kon koff kt ke ft QR R_max];

x0 = [init_receptor init_ligand init_complex];
[tout,xout] = ode45(@(time,state) RLC_kinetics(time,state,sim_para), time_interval, x0);

figure(1); clf;
subplot(311);
plot(tout,xout(:,1))
ylabel('Receptor [nM]');
xlabel('time [min]');
axis([time_interval 0 0.5]);
subplot(312);
plot(tout,xout(:,2))
ylabel('Ligand [nM]');
xlabel('time [min]');
axis([time_interval 0 0.05]);
subplot(313);
plot(tout,xout(:,3))
ylabel('Complex [nM]');
xlabel('time [min]');
axis([time_interval 0 0.004]);

%% 1000 simulations
Num_Sim = 1000;

R_min_all = zeros(1,Num_Sim);
L_min_all = zeros(1,Num_Sim);
C_min_all = zeros(1,Num_Sim);

for idx=1:Num_Sim
    
    idx
    
    init_receptor = rand(1)*0.2; % [nM]
    init_ligand = rand(1)*0.05; %[nM]
    init_complex = rand(1)*0.01; %[kg]
    x0 = [init_receptor init_ligand init_complex];
    [tout,xout] = ode45(@(time,state) RLC_kinetics(time,state,sim_para), time_interval, x0);
    
    xmin = min(xout);
    R_min_all(idx) = xmin(1);
    L_min_all(idx) = xmin(2);
    C_min_all(idx) = xmin(3);
end



%% RLC Kinetics
function dxdt = RLC_kinetics(time,state, sim_para)
    R = state(1);
    L = state(2);
    C = state(3);

    kon = sim_para(1); 
    koff = sim_para(2); 
    kt = sim_para(3);
    ke = sim_para(4);
    ft = sim_para(5);
    QR = sim_para(6);
    R_max = sim_para(7);
    
    if R > R_max         
        QR = 0;
    end
    
    dxdt = zeros(3,1);
    dxdt(1) = -kon*R*L + koff*C - kt*R + QR;
    dxdt(2) = -kon*R*L + koff*C + ft;
    dxdt(3) = kon*R*L - koff*C - ke*C;
end