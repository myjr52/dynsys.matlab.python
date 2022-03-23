#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2022 Jongrae.K

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy import optimize

# Santillan's model delayed differential equation
def Santillan_E_coli_Tryptophan(time, state_all, parameters, T_ext):
    
    state_org = state_all
    state_all[state_all<0] = 0.0
    
    #------------------------------------------------
    # Uncertain parameters
    #------------------------------------------------
    tau_p           = parameters[0]        
    tau_m           = parameters[1]      
    tau_rho         = parameters[2]     
    tau_e           = parameters[3]
    R               = parameters[4]             
    n_H             = parameters[5]          
    b               = parameters[6]            
    e               = parameters[7]             
    f               = parameters[8]             
    O               = parameters[9]
    k_mr            = parameters[10]       
    k_pr            = parameters[11]                  
    k_mi            = parameters[12]      
    k_pi            = parameters[13]     
    k_mt            = parameters[14]        
    k_pt            = parameters[15]  
    c               = parameters[16]            
    d               = parameters[17]           
    gama            = parameters[18]    
    T_consume_rate  = parameters[19]  
    P               = parameters[20]              
    rho             = parameters[21]  
    mu              = parameters[22]

    #----------------------------------
    # Dependent variables
    #----------------------------------
    K_i         = k_mi/k_pi
    K_t         = k_mt/k_pt
    K_r         = k_mr/k_pr

    k_rho       = 1/(tau_rho*rho)
    k_p         = 1/(tau_p*P)
    kdD         = rho*k_rho/30

    #-----------------------------------
    # Steady-state
    #-----------------------------------
    T_SS = K_i
    K_g  = T_SS/10/2   # < (steady-state Tryptophan concentration = 4.1)/10[/2(?)]
    g_SS = T_consume_rate*(K_i + K_g)/K_i
    G_SS = g_SS*K_i/(K_i+K_g)
    
    R_A_SS =  T_SS/(T_SS+K_t)*R
    O_F_SS = (K_r*mu*O)/(K_r*k_p*(1-np.exp(-mu*tau_p))+mu*(K_r+R_A_SS))
    M_F_SS = k_p*P*O_F_SS*np.exp(-mu*tau_m)*(1-b*(1-np.exp(-K_i/c))) \
        /(k_rho*rho*(1-np.exp(-mu*tau_rho))+kdD+mu)
    E_SS = (k_rho*rho*M_F_SS*np.exp(-mu*tau_e))/(2*(gama+mu))
   
    K = 2*(G_SS + mu*K_i)/E_SS
    
    # state
    O_F = state_all[0] 
    M_F = state_all[1] 
    E   = state_all[2]
    T   = state_all[3]

    # delayed state
    state_tau_p     = state_all[4:6];   state_tau_p.resize((2,1))
    state_tau_m     = state_all[6:8];   state_tau_m.resize((2,1))
    state_tau_rho   = state_all[8:10];  state_tau_rho.resize((2,1))
    state_tau_e     = state_all[10::];  state_tau_e.resize((2,1))
    
    A_tau_p     = np.array([[0,1], [-12/tau_p**2,     -6/tau_p]]) 
    A_tau_m     = np.array([[0,1], [-12/tau_m**2,     -6/tau_m]])
    A_tau_rho   = np.array([[0,1], [-12/tau_rho**2,   -6/tau_rho]])
    A_tau_e     = np.array([[0,1], [-12/tau_e**2,     -6/tau_e]])
    B_tau       = np.array([[0], [1]])
    C_tau_p     = np.array([[0, -12/tau_p]])
    C_tau_m     = np.array([[0, -12/tau_m]])
    C_tau_rho   = np.array([[0, -12/tau_rho]])
    C_tau_e     = np.array([[0, -12/tau_e]])
    D_tau       = np.array([[1]])
    
    # dxdt = Ax + Bu
    dO_F_tau_p      = A_tau_p@state_tau_p        + B_tau@np.array([[O_F]])
    dO_F_tau_m      = A_tau_m@state_tau_m        + B_tau@np.array([[O_F]])
    dM_F_tau_rho    = A_tau_rho@state_tau_rho    + B_tau@np.array([[M_F]])
    dM_F_tau_e      = A_tau_e@state_tau_e        + B_tau@np.array([[M_F]])
    
    # y = Cx + Du
    O_F_tau_p   = C_tau_p@state_tau_p        + D_tau@np.array([[O_F]]) 
    O_F_tau_m   = C_tau_m@state_tau_m        + D_tau@np.array([[O_F]]) 
    M_F_tau_rho = C_tau_rho@state_tau_rho    + D_tau@np.array([[M_F]]) 
    M_F_tau_e   = C_tau_e@state_tau_e        + D_tau@np.array([[M_F]]) 

    # make 1x1 array to scalar
    O_F_tau_p=O_F_tau_p[0][0]
    O_F_tau_m=O_F_tau_m[0][0]
    M_F_tau_rho=M_F_tau_rho[0][0]
    M_F_tau_e=M_F_tau_e[0][0]
    
    d_delay_dt = np.vstack((dO_F_tau_p,dO_F_tau_m,dM_F_tau_rho,dM_F_tau_e))
    d_delay_dt = d_delay_dt.squeeze()
    
    # auxilary variables
    A_T = b*(1-np.exp(-T/c))
    E_A = K_i**n_H/(K_i**n_H + T**n_H)*E
    R_A = T/(T+K_t)*R
    G   = g_SS*T/(T+K_g)
    F   = d*T_ext/(e + T_ext*(1+T/f))

    # kinetics
    dOF_dt = K_r/(K_r + R_A)*(mu*O - k_p*P*(O_F - O_F_tau_p*np.exp(-mu*tau_p))) - mu*O_F
    dMF_dt = k_p*P*O_F_tau_m*np.exp(-mu*tau_m)*(1-A_T) \
        - k_rho*rho*(M_F - M_F_tau_rho*np.exp(-mu*tau_rho)) - (kdD + mu)*M_F
    dE_dt = 0.5*k_rho*rho*M_F_tau_e*np.exp(-mu*tau_e) - (gama + mu)*E;
    dT_dt = K*E_A - G + F - mu*T;
    
    if state_org[0] < 0 and dOF_dt < 0:
        dOF_dt = 0
    if state_org[1] < 0 and dMF_dt < 0:
        dMF_dt = 0;
    if state_org[2] < 0 and dE_dt < 0:
        dE_dt = 0
    if state_org[3] < 0 and dT_dt < 0:
        dT_dt = 0
    
    dOF_MF_E_T_dt = np.array([dOF_dt, dMF_dt, dE_dt, dT_dt])
        
    # return all state
    dxdt = np.hstack((dOF_MF_E_T_dt,d_delay_dt))
    
    return dxdt

# uncertain parameters & initial conditions
#
# [Ref] Moisés Santillán and Michael C. Mackey. Dynamic reguiation of the tryptophan 
# operon: A modeling study and comparison with experimental data.
# Proceedings of the National Academy of Sciences, 98(4):1364–1369, February 2001.
#
def Santillans_Tryptophan_Model_constants(delta):

    #------------------------------------------------
    # Uncertain ranges without experimental evidences
    #------------------------------------------------
    Santillan_tau_p   = 0.1*(1 + delta[0])                  # 1
    Santillan_tau_m   = 0.1*(1 + delta[1])                  # 2
    Santillan_tau_rho = 0.05*(1 + delta[2])                 # 3
    Santillan_tau_e   = 0.66*(1 + delta[3])                 # 4
    
    Santillan_R = 0.8*(1 + delta[4])                        # 5

    Santillan_n_H = 2 + delta[5]                            # 6
                                                            # nominal = 1.2
                                                            # delta_nominal = -0.8

    Santillan_b = 0.65 + 0.35*delta[6]                      # 7  [0.3, 1.0]
                                                            # nominal = 0.85
                                                            # delta_nominal = 0.5714

    Santillan_e = 0.9*(1 + delta[7])                        # 8
    Santillan_f = 380*(1 + delta[8])                        # 9
    
    Santillan_O = 3.32e-3*(1 + delta[9])                    # 10
    
    Santillan_k_mr = 1.2e-2*(1 + delta[10])                 # 11 value in [Ref] & its supplementary is different
    Santillan_k_pr = 4.6*(1 + delta[11])                    # 12 value in [Ref] & its supplementary is different
                                                            # but the ratio, kmr/kpr is the same  

    Santillan_k_mi = 7.2e-2*(1 + delta[12])                 # 13
    Santillan_k_pi = 1.76e-2*(1 + delta[13])                # 14

    Santillan_k_mt = 2.1e4*(1 + delta[14])                  # 15
    Santillan_k_pt = 348*(1 + delta[15])                    # 16
    
    Santillan_c = 4e-2*(1 + delta[16])                      # 17
    Santillan_d = 23.5*(1 + delta[17])                      # 18
    
    Santillan_gama = 0.01*(1 + delta[18])                   # 19
                                                            # nominal value 0 
                                                            # delta nominal = -1 
                                                              
    #----------------------------------
    # Uncertain ranges from experiments
    #----------------------------------
    Santillan_T_consume_rate = 21.5 + 7.5*delta[19]         # 20
                                                            # range 14 ~ 29
                                                            # nominal 22.7 -> 0.16

    Santillan_P = 2.785 + 0.675*delta[20]   
                                                            # 21
               # range 2.11 - 3.46 micro-Molar, 
               # nominal 2.6 -> -0.2741
               # 1250 molecule per cell, cell average volume 6.0e-16 - 9.8e-16
               # liters, average volumn = (6.0 + 9.8)/2*1e-16 = 7.9e-16 liters
               # 1250 molecule = 1250/6.022e23 = 2.0757e-21 mole
               # 2.0757e-21/7.9e-16 = 2.62e-6 Molar = 2.62 micro-Molar

    Santillan_rho = 3.12 + 0.75*delta[21] 
                                                            # 21
               # range 2.37 - 3.87 micro-Molar, 
               # nominal 2.9 -> -0.2933
               # 1400 molecule per cell, cell average volume 6.0e-16 - 9.8e-16
               # liters, average volumn = (6.0 + 9.8)/2*1e-16 = 7.9e-16 liters
               # 1400 molecule = 1400/6.022e23 = 2.3248e-21 mole
               # 2.3248e-21/7.9e-16 = 2.94e-6 Molar = 2.94 micro-Molar

    Santillan_mu = 0.0259 + 0.0159*delta[22] 
                                                            # 22
               # range 0.01 ~ 0.0417 [min^-1],
               # nominal 0.01 -> -1
               # actual range from 0.6 h^-1 ~ 2.5 h^-1
    
    # return values
    num_para = 23
    perturbed_model_para = np.zeros(num_para)
    perturbed_model_para[0] = Santillan_tau_p 
    perturbed_model_para[1] = Santillan_tau_m 
    perturbed_model_para[2] = Santillan_tau_rho 
    perturbed_model_para[3] = Santillan_tau_e
    perturbed_model_para[4] = Santillan_R
    perturbed_model_para[5] = Santillan_n_H
    perturbed_model_para[6] = Santillan_b
    perturbed_model_para[7] = Santillan_e
    perturbed_model_para[8] = Santillan_f
    perturbed_model_para[9] = Santillan_O
    perturbed_model_para[10] = Santillan_k_mr
    perturbed_model_para[11] = Santillan_k_pr
    perturbed_model_para[12] = Santillan_k_mi
    perturbed_model_para[13] = Santillan_k_pi
    perturbed_model_para[14] = Santillan_k_mt
    perturbed_model_para[15] = Santillan_k_pt
    perturbed_model_para[16] = Santillan_c
    perturbed_model_para[17] = Santillan_d
    perturbed_model_para[18] = Santillan_gama
    perturbed_model_para[19] = Santillan_T_consume_rate
    perturbed_model_para[20] = Santillan_P
    perturbed_model_para[21] = Santillan_rho
    perturbed_model_para[22] = Santillan_mu

    return perturbed_model_para

# check negative states to stop the integrator
def negativeConcentration(time,state,parameters, T_ext):
    tol = -1e-1;
    OF_MF_E_T = state[0:4]
    delay_output = state[4::2]
    all_positive_state = np.hstack((OF_MF_E_T,delay_output))
    value = 1-float(any(all_positive_state<tol))
    return value

# Cost function for the model fitting
def Santillan_Model_Fit_Cost(delta, tspan, time_exp, Act_Enzy_exp, plot_sw):

    try: 
        num_state = 12;
        model_para = Santillans_Tryptophan_Model_constants(delta);
        
        negativeConcentration.terminal = True
        negativeConcentration.direction = 0
        
        # Initially the culture in the medium with presence of the external tryptophan
        T_ext = 400*(model_para[12]/model_para[13]); # 400 times of T(t) steady-state
        time_eval = np.linspace(tspan[0],tspan[1],1000)
        state_t0 = np.zeros(num_state)
        
        sol = solve_ivp(Santillan_E_coli_Tryptophan, tspan,
                    state_t0, events=negativeConcentration, args=(model_para, T_ext), 
                    t_eval=time_eval, rtol=1e-6, atol=1e-9)
        OF_MF_E_T_IC = np.mean(sol.y[:,-50:-1],axis=1) # it reaches to the steady-state
           
        # No external tryptophan medium shift experiment
        T_ext = 0
        state_t0=OF_MF_E_T_IC # the steady state becomes the initial condition
        sol = solve_ivp(Santillan_E_coli_Tryptophan, (tspan[0], time_exp[-1]),  
                    state_t0, args=(model_para, T_ext), 
                    t_eval=time_exp, rtol=1e-3, atol=1e-6)
        
        # evaluate the Enzyme and the Tryptophan at the given measurent time
        state_at_time_exp = sol.y[0:4,:]
        E_at_time_exp = state_at_time_exp[2,:]
        T_at_time_exp = state_at_time_exp[3,:]
        
        # calculate the active enzyme using the model
        n_H = model_para[5]
        K_i = model_para[12]/model_para[13]
        EA_model = (K_i**n_H/(K_i**n_H + T_at_time_exp**n_H))*E_at_time_exp    
        
        # normalize the active enzyme
        y_bar = EA_model/EA_model[-1]
        y_tilde = Act_Enzy_exp/Act_Enzy_exp[-1]
        
        if plot_sw:
            time_eval = np.linspace(tspan[0],time_exp[-1],500)
            sol = solve_ivp(Santillan_E_coli_Tryptophan, (tspan[0], time_exp[-1]), 
                            state_t0, events=negativeConcentration, args=(model_para, T_ext), 
                            t_eval=time_eval, rtol=1e-3, atol=1e-6)
            EA_full = (K_i**n_H/(K_i**n_H + sol.y[3,:]**n_H))*sol.y[2,:]
            
            plt.plot(sol.t,EA_full/EA_full[-1],'b-')
            plt.plot(time_exp,y_tilde,'rx')
            plt.axis([0, 180, 0, 2])
            plt.xlabel('time [min]',fontsize=20)
            plt.ylabel('Normalized Active Enzyme',fontsize=20)
            plt.xticks([0,50,100,150,180],fontsize=20)
            plt.yticks([0,0.5,1.0,1.5,2.0],fontsize=20)
            plt.grid()
            #plt.savefig('python_E_Coli_Exp_C.pdf',bbox_inches='tight')
            
        # calculate the cost
        J_cost = np.sum((y_bar-y_tilde)**2)       
        
    except:
        J_cost = 1e3
    
    return J_cost

#-----------------------------------------
# E. coli Active Enzyme Experiment data
#-----------------------------------------

# [Ref] Charles Yanofsky and Virginia Horn. Role of regulatory features of the trp
# operon of Escherichia coli in mediating a response to a nutritional shift.
# Journal of Bacteriology, 176(20):6245–6254, October 1994

experiment_num = 2 # 1(A), 2(B), 3(C)

time_A = np.array([0, 20, 38, 59, 89, 119, 149])
Enzy_A = np.array([25, 657, 617, 618, 577, 577, 567])

time_B = np.array([0, 29, 60, 89, 179])
Enzy_B = np.array([0, 1370, 1362, 1291, 913])

time_C = np.array([0, 29, 58, 88, 118, 178])
Enzy_C = np.array([0, 754, 888, 763, 704, 683])

# choose experiment
if experiment_num==1:
        time_exp = time_A 
        Enzy_exp = Enzy_A
elif experiment_num==2:
        time_exp = time_B 
        Enzy_exp = Enzy_B
elif experiment_num==3: 
        time_exp = time_C 
        Enzy_exp = Enzy_C     

#-----------------------------------------
# Main Model Fitting Optimization
#-----------------------------------------
delta_dim = 23;

# time span for obtaining the steady-state
time_span = np.array([0, 1200]) # [minutes] 

state_all = np.random.randn(12)
delta = 0.99*(2*np.random.rand(23)-1)

Act_Enzy_exp = Enzy_exp
plot_sw = False
bounds = [(-0.99,0.99)]*delta_dim

# result = optimize.differential_evolution(Santillan_Model_Fit_Cost, bounds, 
#                                           args=(time_span, time_exp, Act_Enzy_exp, plot_sw), 
#                                           updating='deferred', disp=True, popsize=200, maxiter=100, workers=4)

# with open('exp_B_0_99bnds.npy', 'wb') as f:
#       np.save(f, result)     

# plot results
if experiment_num==1:
    ########################AAAAAAAAAAAAAAAAAAAAAAAAAAAA
    opt_result = np.load('exp_A_0_99bnds.npy',allow_pickle=True)
    delta_best = opt_result.item(0).x
    delta_best = np.array([ 0.4139, -0.6501, -0.5079,  0.6557,  0.8028, -0.6979, -0.0143,
       -0.7478, -0.6088, -0.4282, -0.9513,  0.9848,  0.4909, -0.6512,
        0.1611, -0.4746, -0.0649,  0.4922, -0.5255, -0.4848, -0.6853,
        0.2661,  0.1039])
    # delta_best = np.array([ 0.41391058, -0.65007115, -0.50789937,  0.65569299,  0.80284194,
    #     -0.69788653, -0.01426885, -0.7477905 , -0.60884596, -0.42819439,
    #     -0.95128225,  0.98480601,  0.49091824, -0.6512251 ,  0.16107301,
    #     -0.47463868, -0.06487504,  0.49217244, -0.52551436, -0.48475391,
    #     -0.68534604,  0.26606286,  0.10390048])
elif experiment_num==2:
    ########################BBBBBBBBBBBBBBBBBBBBBBBBBBBB
    opt_result = np.load('exp_B_0_99bnds.npy',allow_pickle=True)
    delta_best = opt_result.item(0).x
    # delta_best = np.array([-0.8165, -0.9378, -0.279 ,  0.9639,  0.99  , -0.9403,  0.7599,
    #     0.3749,  0.9659,  0.0057, -0.9891,  0.5645,  0.3566, -0.6389,
    #     -0.4876, -0.4666, -0.3607,  0.8619, -0.6752, -0.0358,  0.2272,
    #     0.417 ,  0.3176])
    delta_best = np.array([-0.8164609821840342   , -0.9377737042517797   ,
        -0.2790376682215831   ,  0.9639037189678477   ,
        0.9899716285333048   , -0.9403394600209112   ,
        0.7598769991256485   ,  0.3748592645823837   ,
        0.9659493912171094   ,  0.0057193188983498434,
        -0.9890715372649793   ,  0.5645345097583515   ,
        0.35655407244896503  , -0.6388522834084576   ,
        -0.4876497370209529   , -0.4666376552759634   ,
        -0.3606741179703174   ,  0.8618508223375989   ,
        -0.6752330337194928   , -0.03575084441063797  ,
        0.22717349615134524  ,  0.4169772469116487   ,
        0.3175712292813526])
elif experiment_num==3:
    ######################CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    opt_result = np.load('exp_C_0_99bnds.npy',allow_pickle=True)
    delta_best = opt_result.item(0).x
    delta_best = np.array([-9.5300e-01, -1.9603e-01,  8.5486e-01,  7.6837e-01, -8.4957e-04,
       -9.6020e-01, -8.8099e-02, -7.1178e-01,  9.4534e-01,  3.3016e-02,
       -9.7867e-01,  6.2270e-02,  2.5430e-01, -7.5352e-01,  6.1063e-01,
       -8.7082e-01,  5.9102e-01,  8.5359e-01, -1.0110e-02, -8.9946e-01,
        8.2926e-01,  2.6782e-01,  2.1704e-01])
    # delta_best = np.array([-9.52997307e-01, -1.96025254e-01,  8.54855055e-01,  7.68370196e-01,
    #     -8.49567218e-04, -9.60202292e-01, -8.80987283e-02, -7.11775804e-01,
    #     9.45341324e-01,  3.30157377e-02, -9.78673243e-01,  6.22704640e-02,
    #     2.54295149e-01, -7.53516288e-01,  6.10627689e-01, -8.70823889e-01,
    #     5.91017184e-01,  8.53594505e-01, -1.01095725e-02, -8.99463857e-01,
    #     8.29264858e-01,  2.67819560e-01,  2.17044466e-01])
    

J_opt=Santillan_Model_Fit_Cost(delta_best,time_span,time_exp,Enzy_exp,True)
print(J_opt)
print(Santillans_Tryptophan_Model_constants(delta_best))
