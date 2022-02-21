clear

syms kon kcat kdg Kdeg kI eta gamma_G kP ks2 ks1 Pd X3 real;
syms E P ES X1 X2 S DP Xs real;


dE_dt = -kon*E*S + kcat*ES;
dP_dt = kcat*ES - kdg*P - Kdeg*X3*P;
dES_dt = kon*E*S - kcat*ES;
dX1_dt = kI*DP - eta*X1;
dX2_dt = -gamma_G*X2 + gamma_G*kP*DP;
dS_dt = ks2*X1 + ks2*X2 - ks2*S;
dDP_dt = ks1*Pd - ks1*DP*Xs - ks1*DP;
dXs_dt = -ks1*DP*Xs + ks1*P;

fx = [dE_dt; dP_dt; dES_dt; dX1_dt; dX2_dt; dS_dt; dDP_dt; dXs_dt];
state = [E; P; ES; X1; X2; S; DP; Xs];

dfdx = jacobian(fx,state);

%% Steady-state
Ess = 0.1998;    
Pss = 0.4999;    
ESss = 0.0002;    
X1ss = 0.0558;   
X2ss = 50.0415;   
Sss = 50.0723;    
DPss = 1.0001;    
Xsss = 0.4999;

dfdx_at_ss = subs(dfdx,{E, P, ES, X1, X2, S, DP, Xs},{Ess, Pss, ESss, X1ss, X2ss, Sss, DPss, Xsss});

%% nomial stability with the nominal parameters
kP = 50;
kI = 5e-6;
gamma_G = 8e-4;
ks2 = 4e-4;
Kdeg = 1e-3;
X3 = 1; 
KF = 3; 
eta = 1e-4;            
Pd = 1;
kon = 5e-5;
kcat = 1.6*2;
kdg = 8e-8;
ks1 = 3;


dfdx_nominal = subs(dfdx_at_ss, ...
    {sym('kP'), sym('kI'),sym('gamma_G'),sym('ks2'),sym('Kdeg'), ...
    sym('X3'),sym('KF'),sym('eta'),sym('Pd'),sym('kon'),sym('kcat'),sym('kdg'),sym('ks1')}, ...
    {kP,kI,gamma_G,ks2,Kdeg,X3,KF,eta,Pd,kon,kcat,kdg,ks1});

dfdx_nominal_val = eval(dfdx_nominal);