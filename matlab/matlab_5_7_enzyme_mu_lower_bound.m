clear

syms kon kcat kdg Kdeg kI eta gamma_G kP ks2 ks1 Pd X3 real;
syms E P ES X1 X2 S DP Xs real;
syms d1 d2 d3 d4 d5 d6 d7 d8 d9 d10 d11 d12 d13 d14 d15 d16 real;


dE_dt = -(kon+d1)*E*S + (kcat+d2)*ES;
dP_dt = (kcat+d2)*ES - (kdg+d3)*P - (Kdeg+d4)*(X3+d5)*P;
dES_dt = (kon+d1)*E*S - (kcat+d2)*ES;
dX1_dt = (kI+d6)*DP - (eta+d7)*X1;
dX2_dt = -(gamma_G+d8)*X2 + (gamma_G+d9)*(kP+d10)*DP;
dS_dt = (ks2+d11)*X1 + (ks2+d12)*X2 - (ks2+d13)*S;
dDP_dt = (ks1+d14)*Pd - (ks1+d15)*DP*Xs - (ks1+d16)*DP;
dXs_dt = -(ks1+d15)*DP*Xs + (ks1+d16)*P;

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


%% mu-analysis
Ns = 5000;
eps = 1e-6;

num_state = 8;
num_delta = 16;
A0 =  eval(subs(dfdx_nominal,{d1 d2 d3 d4 d5 d6 d7 d8 d9 d10 d11 d12 d13 d14 d15 d16},{zeros(1,16)}));

num_omega = 10;
omega_all = [0 logspace(-3,-1,num_omega)];
num_omega = num_omega + 1;

mu_lb = zeros(1,num_omega);

%% lower bound using geomatric approach
for wdx=1:num_omega
    omega = omega_all(wdx);
    Mjw = inv(1j*omega*eye(num_state)-A0);

    d_lb = 1e-6;
    d_ub = 10;
    d_ulb = d_ub - d_lb;
    
    if omega==0
        size_check = 2;
    else
        size_check = 4;
    end

    while d_ulb > eps

        d = (d_lb+d_ub)/2;

        sign_all = [];

        for idx=1:Ns
            delta_vec = rand(1,num_delta)*d-d/2;
            rand_face = randi(num_delta,1);
            delta_vec(rand_face) = d/2;

            Delta = eval(subs(dfdx_nominal,{d1 d2 d3 d4 d5 d6 d7 d8 d9 d10 d11 d12 d13 d14 d15 d16}, ...
                {delta_vec})) -A0;

            I_MD = det(eye(num_state)-Mjw*Delta);
            fR = sign(real(I_MD));
            fI = sign(imag(I_MD));

            sign_all = unique([sign_all; fR fI],'row');

        end
        
        if size(sign_all,1) == size_check
            d_ub = d;
        else
            d_lb = d;
        end
        
        d_ulb = d_ub - d_lb;
        
%         omega
%         sign_all
%         [d_lb d_ub]

    end

    mu_lb(wdx) = 2/d_ub;
end