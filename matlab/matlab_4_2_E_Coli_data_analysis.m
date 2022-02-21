clear

%% E. coli Active Enzyme Experiment data

% [Ref] Charles Yanofsky and Virginia Horn. Role of regulatory features of the trp
% operon of Escherichia coli in mediating a response to a nutritional shift.
% Journal of Bacteriology, 176(20):6245–6254, October 1994

time_A = [0 20 38 59 89 119 149];
Enzy_A = [25 657 617 618 577 577 567];

time_B = [0 29 60 89 179];
Enzy_B = [0 1370 1362 1291 913];

time_C = [0 29 58 88 118 178];
Enzy_C = [0 754 888 763 704 683];

time_Text = [0 1200];

% experiment A
load exp_A_0_99_delta_best;
Santillan_Model_Fit_Cost(delta_best_A, time_Text, time_A, Enzy_A, true)

% experiment B
load exp_B_0_99_delta_best;
Santillan_Model_Fit_Cost(delta_best_B, time_Text, time_B, Enzy_B, true)

% experiment C
load exp_C_0_99_delta_best;
Santillan_Model_Fit_Cost(delta_best_C, time_Text, time_C, Enzy_C, true)

% noise strength
delta_best_all = [delta_best_A; delta_best_B; delta_best_C];
delta_rank=(std(delta_best_all).^2)./abs(mean(delta_best_all));
delta_rank(isnan(delta_rank))=0;
[delta_rank_val,delta_rank_idx]=sort(delta_rank);

%% Cost function for the model fitting
function J_cost = Santillan_Model_Fit_Cost(delta, tspan, time_exp, Act_Enzy_exp, plot_sw)

    if nargin < 5
        plot_sw = false;
    end
    
    try
        num_state = 12;
        model_para = Santillans_Tryptophan_Model_constants(delta);
        
        % Initially the culture in the medium with presence of the external tryptophan
        T_ext = 400*(model_para(13)/model_para(14)); % 400 times of T(t) steady-state
        ode_option = odeset('RelTol',1e-3,'AbsTol',1e-6,'Events',@negativeConcentration);
        state_t0 = zeros(1,num_state);
        sol = ode45(@(time, state)Santillan_E_coli_Tryptophan(time, state, ...
            model_para, T_ext),tspan, state_t0, ode_option);
        OF_MF_E_T_IC = mean(sol.y(:,end-50:end),2); % it reaches to the steady-state
        
        sol2=sol;
        
        % No external tryptophan medium shift experiment
        T_ext = 0;
        state_t0 =  OF_MF_E_T_IC(:); % the steady state becomes the initial condition
        tspan_sim = [0 time_exp(end)];
        sol = ode45(@(time, state)Santillan_E_coli_Tryptophan(time, state, ...
            model_para, T_ext), tspan_sim, state_t0, ode_option);
        
        % evaluate the Enzyme and the Tryptophan at the given measurent time
        state_at_time_exp = deval(sol,time_exp);
        E_at_time_exp = state_at_time_exp(3,:);
        T_at_time_exp = state_at_time_exp(4,:);
        
        % calculate the active enzyme using the model
        n_H = model_para(6);
        K_i = model_para(13)/model_para(14);
        EA_model = (K_i^n_H./(K_i^n_H + T_at_time_exp.^n_H)).*E_at_time_exp;
        
        % normalize the active enzyme
        y_bar = EA_model/EA_model(end);
        y_tilde = Act_Enzy_exp/Act_Enzy_exp(end);
        
        if plot_sw
            figure;
            EA_full = (K_i^n_H./(K_i^n_H + sol.y(4,:).^n_H)).*sol.y(3,:);
            plot(sol.x,EA_full/EA_full(end),'b-');
            hold on;
            plot(time_exp,y_tilde,'rx');
            set(gca,'FontSize',12);
            xlabel('time [min]');
            ylabel('Normalized Active Enzyme');
            axis([0 150 0 2]);
            grid on;
        end
        
        % calculate the cost
        J_cost = sum((y_bar-y_tilde).^2);
        
    catch
        J_cost = 1e3;
    end
    
end

%% Santillan's model delayed differential equation
function dxdt = Santillan_E_coli_Tryptophan(time, state_all, parameters, T_ext)
    
    state_org = state_all;
    state_all(state_all<0) = 0.0;

    %------------------------------------------------
    % Uncertain parameters
    %------------------------------------------------
    tau_p           = parameters(1);        
    tau_m           = parameters(2);      
    tau_rho         = parameters(3);     
    tau_e           = parameters(4);
    R               = parameters(5);             
    n_H             = parameters(6);          
    b               = parameters(7);            
    e               = parameters(8);             
    f               = parameters(9);             
    O               = parameters(10);
    k_mr            = parameters(11);       
    k_pr            = parameters(12);                  
    k_mi            = parameters(13);      
    k_pi            = parameters(14);     
    k_mt            = parameters(15);        
    k_pt            = parameters(16);  
    c               = parameters(17);            
    d               = parameters(18);           
    gama            = parameters(19);    
    T_consume_rate  = parameters(20);  
    P               = parameters(21);              
    rho             = parameters(22);  
    mu              = parameters(23);

    %----------------------------------
    % Dependent variables
    %----------------------------------
    K_i         = k_mi/k_pi;
    K_t         = k_mt/k_pt;
    K_r         = k_mr/k_pr;

    k_rho       = 1/(tau_rho*rho);
    k_p         = 1/(tau_p*P);
    kdD         = rho*k_rho/30;

    %-----------------------------------
    % Steady-state
    %-----------------------------------
    T_SS = K_i;
    K_g  = T_SS/10/2;  % < (steady-state Tryptophan concentration = 4.1)/10[/2(?)]
    g_SS = T_consume_rate*(K_i + K_g)/K_i;
    G_SS = g_SS*K_i/(K_i+K_g);
    
    R_A_SS =  T_SS/(T_SS+K_t)*R;
    O_F_SS = (K_r*mu*O)/(K_r*k_p*(1-exp(-mu*tau_p))+mu*(K_r+R_A_SS));
    M_F_SS = k_p*P*O_F_SS*exp(-mu*tau_m)*(1-b*(1-exp(-K_i/c))) ...
            /(k_rho*rho*(1-exp(-mu*tau_rho))+kdD+mu);
    E_SS = (k_rho*rho*M_F_SS*exp(-mu*tau_e))/(2*(gama+mu));
   
    K = 2*(G_SS + mu*K_i)/E_SS;
    
    % state
    O_F = state_all(1); 
    M_F = state_all(2); 
    E   = state_all(3);
    T   = state_all(4);

    % delayed state
    state_tau_p     = state_all(5:6);
    state_tau_m     = state_all(7:8);
    state_tau_rho   = state_all(9:10);
    state_tau_e     = state_all(11:12);
    
    A_tau_p     = [0 1; -12/tau_p^2     -6/tau_p]; 
    A_tau_m     = [0 1; -12/tau_m^2     -6/tau_m];
    A_tau_rho   = [0 1; -12/tau_rho^2   -6/tau_rho];
    A_tau_e     = [0 1; -12/tau_e^2     -6/tau_e];
    B_tau       = [0; 1];
    C_tau_p     = [0 -12/tau_p];
    C_tau_m     = [0 -12/tau_m];
    C_tau_rho   = [0 -12/tau_rho];
    C_tau_e     = [0 -12/tau_e];
    D_tau       = 1;
    
    % dxdt = Ax + Bu
    dO_F_tau_p      = A_tau_p*state_tau_p(:)        + B_tau*O_F;
    dO_F_tau_m      = A_tau_m*state_tau_m(:)        + B_tau*O_F;
    dM_F_tau_rho    = A_tau_rho*state_tau_rho(:)    + B_tau*M_F;
    dM_F_tau_e      = A_tau_e*state_tau_e(:)        + B_tau*M_F;
    
    % y = Cx + Du
    O_F_tau_p   = C_tau_p*state_tau_p(:)        + D_tau*O_F;
    O_F_tau_m   = C_tau_m*state_tau_m(:)        + D_tau*O_F;
    M_F_tau_rho = C_tau_rho*state_tau_rho(:)    + D_tau*M_F;
    M_F_tau_e   = C_tau_e*state_tau_e(:)        + D_tau*M_F;
    
    d_delay_dt = [dO_F_tau_p(:); dO_F_tau_m(:); dM_F_tau_rho(:); dM_F_tau_e(:)];
    
    % auxilary variables
    A_T = b*(1-exp(-T/c));
    E_A = K_i^n_H/(K_i^n_H + T^n_H)*E;
    R_A = T/(T+K_t)*R;
    G   = g_SS*T/(T+K_g);
    F   = d*T_ext/(e + T_ext*(1+T/f));

    % kinetics
    dOF_dt = K_r/(K_r + R_A)*(mu*O - k_p*P*(O_F - O_F_tau_p*exp(-mu*tau_p))) - mu*O_F;
    dMF_dt = k_p*P*O_F_tau_m*exp(-mu*tau_m)*(1-A_T) ... 
        - k_rho*rho*(M_F - M_F_tau_rho*exp(-mu*tau_rho)) - (kdD + mu)*M_F;
    dE_dt = 0.5*k_rho*rho*M_F_tau_e*exp(-mu*tau_e) - (gama + mu)*E;
    dT_dt = K*E_A - G + F - mu*T;
    
    if state_org(1) < 0 && dOF_dt < 0
        dOF_dt = 0;
    end
    if state_org(2) < 0 && dMF_dt < 0
        dMF_dt = 0;
    end
    if state_org(3) < 0 && dE_dt < 0
        dE_dt = 0;
    end    
    if state_org(4) < 0 && dT_dt < 0
        dT_dt = 0;
    end
    
    dOF_MF_E_T_dt = [dOF_dt dMF_dt dE_dt dT_dt]';
        
    % return all state
    dxdt = [dOF_MF_E_T_dt; d_delay_dt];
    
end



%% uncertain parameters & initial conditions
function [perturbed_model_para] = Santillans_Tryptophan_Model_constants(delta)

    %------------------------------------------------
    % Uncertain ranges without experimental evidences
    %------------------------------------------------
    Santillan_tau_p   = 0.1*(1 + delta(1));                 % 1
    Santillan_tau_m   = 0.1*(1 + delta(2));                 % 2
    Santillan_tau_rho = 0.05*(1 + delta(3));                % 3
    Santillan_tau_e   = 0.66*(1 + delta(4));                % 4
    
    Santillan_R = 0.8*(1 + delta(5));                       % 5

    Santillan_n_H = 2 + delta(6);                           % 6
                                                            % nominal = 1.2
                                                            % delta_nominal = -0.8

    Santillan_b =  0.65 + 0.35*delta(7);                    % 7  [0.3, 1.0]
                                                            % nominal = 0.85
                                                            % delta_nominal = 0.5714

    Santillan_e = 0.9*(1 + delta(8));                       % 8
    Santillan_f = 380*(1 + delta(9));                       % 9
    
    Santillan_O = 3.32e-3*(1 + delta(10));                  % 10
    
    Santillan_k_mr = 1.2e-2*(1 + delta(11));                % 11 value in [Ref] & its supplementary is different
    Santillan_k_pr = 4.6*(1 + delta(12));                   % 12 value in [Ref] & its supplementary is different
                                                            % but the ratio, kmr/kpr is the same  

    Santillan_k_mi = 7.2e-2*(1 + delta(13));                % 13
    Santillan_k_pi = 1.76e-2*(1 + delta(14));               % 14

    Santillan_k_mt = 2.1e4*(1 + delta(15));                 % 15
    Santillan_k_pt = 348*(1 + delta(16));                   % 16
    
    Santillan_c = 4e-2*(1 + delta(17));                     % 17
    Santillan_d = 23.5*(1 + delta(18));                     % 18
    
    Santillan_gama = 0.01*(1 + delta(19));                  % 19
                                                            % nominal value 0 
                                                            % delta nominal = -1 
                                                              
    %----------------------------------
    % Uncertain ranges from experiments
    %----------------------------------
    Santillan_T_consume_rate = 21.5 + 7.5*delta(20);        % 20
                                                            % range 14 ~ 29
                                                            % nominal 22.7 -> 0.16

    Santillan_P = 2.785 + 0.675*delta(21);   
                                                            % 21
               % range 2.11 - 3.46 micro-Molar, 
               % nominal 2.6 -> -0.2741
               % 1250 molecule per cell, cell average volume 6.0e-16 - 9.8e-16
               % liters, average volumn = (6.0 + 9.8)/2*1e-16 = 7.9e-16 liters
               % 1250 molecule = 1250/6.022e23 = 2.0757e-21 mole
               % 2.0757e-21/7.9e-16 = 2.62e-6 Molar = 2.62 micro-Molar

    Santillan_rho = 3.12 + 0.75*delta(22); 
                                                            % 22
               % range 2.37 - 3.87 micro-Molar, 
               % nominal 2.9 -> -0.2933
               % 1400 molecule per cell, cell average volume 6.0e-16 - 9.8e-16
               % liters, average volumn = (6.0 + 9.8)/2*1e-16 = 7.9e-16 liters
               % 1400 molecule = 1400/6.022e23 = 2.3248e-21 mole
               % 2.3248e-21/7.9e-16 = 2.94e-6 Molar = 2.94 micro-Molar

    Santillan_mu = 0.0259 + 0.0159*delta(23); 
                                                            % 23
               % range 0.01 ~ 0.0418 [min^-1],
               % nominal 0.01 -> -1
               % actual range from 0.6 h^-1 ~ 2.5 h^-1
    
    %% return values
    num_para = 23;
    perturbed_model_para = zeros(1,num_para);
    perturbed_model_para(1) = Santillan_tau_p; 
    perturbed_model_para(2) = Santillan_tau_m; 
    perturbed_model_para(3) = Santillan_tau_rho; 
    perturbed_model_para(4) = Santillan_tau_e;
    perturbed_model_para(5) = Santillan_R;
    perturbed_model_para(6) = Santillan_n_H;
    perturbed_model_para(7) = Santillan_b;
    perturbed_model_para(8) = Santillan_e;
    perturbed_model_para(9) = Santillan_f;
    perturbed_model_para(10) = Santillan_O;
    perturbed_model_para(11) = Santillan_k_mr;
    perturbed_model_para(12) = Santillan_k_pr;
    perturbed_model_para(13) = Santillan_k_mi;
    perturbed_model_para(14) = Santillan_k_pi;
    perturbed_model_para(15) = Santillan_k_mt;
    perturbed_model_para(16) = Santillan_k_pt;
    perturbed_model_para(17) = Santillan_c;
    perturbed_model_para(18) = Santillan_d;
    perturbed_model_para(19) = Santillan_gama;
    perturbed_model_para(20) = Santillan_T_consume_rate;
    perturbed_model_para(21) = Santillan_P;
    perturbed_model_para(22) = Santillan_rho;
    perturbed_model_para(23) = Santillan_mu;

end

function [value, isterminal, direction] = negativeConcentration(~,state)
    tol = 1e-1;
    OF_MF_E_T = state(1:4);
    delay_output = state(5:2:11);
    all_positive_state = [OF_MF_E_T(:)' delay_output(:)'];
    value = any(all_positive_state<tol)-1;
    isterminal = 1;
    direction = 0;
end
