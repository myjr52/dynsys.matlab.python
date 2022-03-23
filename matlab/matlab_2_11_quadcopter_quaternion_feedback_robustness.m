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

%% simulation parameters
init_time = 0; % [s]
final_time = 120.0; % [s]
time_interval = [init_time final_time];

J_inertia = [0.005 -0.001  0.004;
            -0.001  0.006 -0.002;
             0.004 -0.002  0.004]; % vehicle moment of inertia [kg m^2]
J_inv = inv(J_inertia);
 
C_T = 8.8e-7; % motor thruster coefficient [N/(rad/s)^2]
C_D = 11.3e-8; % motor drag coefficient [Nm/(rad/s)^2]w_motor_fblr_squared_desired(w_motor_fblr_squared_desired<0) = 0;
L_arm = 0.127; % length from the centre of quadcopter to the motor [m]

q0 = (1/sqrt(4))*[1 1 -1 1]'; % initial quaternion
w0 = [0.1 -0.2 0.1]'; % initial angular velocity

r0 = [0 0 -50]';
v0 = [0 0 0]';

state_0 = [q0; w0; r0; v0]; % states including q0 and omega0

Num_MC = 6;
ts_all = inf(1,Num_MC);
dJ_norm_all = inf(1,Num_MC);

%% simulation
ode_options = odeset('RelTol',1e-6,'AbsTol',1e-9, 'MaxStep', 0.01);

parfor g_MC_idx = 1:Num_MC
    not_find_dJ = true;

    while not_find_dJ

        dJ = diag(0.002*randn(3,1));
        J_inertia_perturbed = J_inertia + dJ;
        pd_cond = min(eig(J_inertia_perturbed))>0;
        j3_cond = J_inertia_perturbed(1,1)+J_inertia_perturbed(2,2) > J_inertia_perturbed(3,3);
        j2_cond = J_inertia_perturbed(1,1)+J_inertia_perturbed(3,3) > J_inertia_perturbed(2,2);
        j1_cond = J_inertia_perturbed(2,2)+J_inertia_perturbed(3,3) > J_inertia_perturbed(1,1);

        if pd_cond && j1_cond && j2_cond && j3_cond
            not_find_dJ = false;
        end

    end
    dJ_norm_all(g_MC_idx) = norm(dJ);
    
    J_inv_perturbed = inv(J_inertia_perturbed);

    quadcopter_uav = {J_inertia_perturbed, J_inv, C_T, C_D, L_arm};

    [tout,state_out] = ode45(@(time,state) dqdt_dwdt_drvdt(time,state,quadcopter_uav), ...
        time_interval, state_0, ode_options);

    qout = state_out(:,1:4);
    wout = state_out(:,5:7);
    rout = state_out(:,8:10);
    vout = state_out(:,11:13);

    q13 = qout(:,1:3);
    q13_norm = sqrt(sum(q13.^2,2));
    q13_ts = int32(q13_norm>0.01);
    q13_ts = cumsum(q13_ts);
    q13_ts = tout(q13_ts==q13_ts(end));
    q13_ts = q13_ts(1);
    ts_all(g_MC_idx) = q13_ts;

    fprintf('#%1d %6.5f %4.2f\n',g_MC_idx,norm(dJ),q13_ts);
end

% end monte-carlo simulation

%% plot
figure;
subplot(121)
plot(dJ_norm_all*1e3, ts_all,'.');
set(gca,'FontSize',12);
ylabel('$t_s$ [s]','interpreter','latex');
xlabel('$||\Delta J||\times10^{-3}$ [kg m$^2$]','interpreter','latex');
axis([0 8 35 90])
subplot(222)
histogram(dJ_norm_all*1e3,'Normalization','probability');
set(gca,'FontSize',12);
xlabel('$||\Delta J||\times10^{-3}$ [kg m$^2$]','interpreter','latex');
ylabel('pdf');
axis([0 8 0 0.1])
subplot(224)
histogram(ts_all,'Normalization','probability');
set(gca,'FontSize',12);
xlabel('$t_s$ [s]','interpreter','latex');
ylabel('pdf');
axis([38 90 0 0.15])

%% functions
function dstate_dt = dqdt_dwdt_drvdt(time,state,quadcopter_uav)

    m_quadcopter = 0.49; %[kg]
    grv_acce = 9.81; %[m/s^2]
    
    q_current = state(1:4);
    q_current = q_current(:)/norm(q_current);
    
    w_current = state(5:7);
    w_current = w_current(:);
    
    rv_current = state(8:13);
    rv_current = rv_current(:);
    
    J_inertia = quadcopter_uav{1};
    inv_J = quadcopter_uav{2};
    C_T = quadcopter_uav{3};
    C_D = quadcopter_uav{4};
    L_arm = quadcopter_uav{5};
    
    q_13 = q_current(1:3); q_13 =q_13(:); q4 = q_current(4);
    q13x = [0 -q_13(3) q_13(2);
            q_13(3) 0 -q_13(1);
            -q_13(2) q_13(1) 0];
    C_BR = (q4^2-q_13'*q_13)*eye(3) + 2*(q_13*q_13')-2*q4*q13x;
    
    %--------------------------------
    % Begin: this part is controller
    %--------------------------------
    w_motor_fblr_desired = quaternion_feedback_and_altitude_control(q_current, ...
        w_current, rv_current, J_inertia, C_T, C_D, L_arm, C_BR, m_quadcopter, grv_acce);
    %--------------------------------
    % End: this part is controller
    %--------------------------------
    
    % Motor Force & Torque
    FM_w_Motor = propeller_motor_actuator(C_T, C_D, L_arm, w_motor_fblr_desired);
    
    F_motor = [0; 0; FM_w_Motor(1)];
    M_torque = FM_w_Motor(2:4);
    
    motor_force_in_R = C_BR'*F_motor;
    
    % Kinematics & Dynamics
    dqdt = dqdt_attitude_kinematics(q_current,w_current);
    dwdt = dwdt_attitude_dynamics(w_current, J_inertia, inv_J, M_torque);
    drvdt = drdt_linear_dynamics(rv_current, m_quadcopter, grv_acce, motor_force_in_R);
    
    dstate_dt = [dqdt(:); dwdt(:); drvdt(:)];

end

function dqdt = dqdt_attitude_kinematics(q_true,w_true)
    q_true = q_true(:); 
    w_true = w_true(:);
    
    wx = [  0           -w_true(3)  w_true(2);
            w_true(3)   0           -w_true(1);
            -w_true(2)  w_true(1)   0];
        
    Omega = [   -wx         w_true;
                -w_true'    0];
            
    dqdt = 0.5*Omega*q_true;
end

function dwdt = dwdt_attitude_dynamics(w_true, J_inertia, inv_J_inertia, M_torque)

    w_true = w_true(:);
    Jw = J_inertia*w_true;
    
    Jw_dot = -cross(w_true,Jw) + M_torque(:);
    
    dwdt = inv_J_inertia*Jw_dot;
end

function drvdt = drdt_linear_dynamics(rv_true, mass, grv_const, motor_force_in_R)
    r = rv_true(1:3); r=r(:);
    v = rv_true(4:6); v=v(:);
    
    drdt = v(:);
    dvdt = [0;0;grv_const] + motor_force_in_R(:)/mass;
    
    drvdt = [drdt; dvdt];
end

function FM_w_Motor = propeller_motor_actuator(C_T,C_D,L_arm,w_command)

    % assume perfect motor angular velocity control
    w_motor = w_command(:);

    F_fblr = C_T*(w_motor.^2);
    tau_fblr = C_D*(w_motor.^2);
    
    F_motor = -sum(F_fblr);
    M_motor = [ L_arm*(F_fblr(3)-F_fblr(4));
                L_arm*(F_fblr(1)-F_fblr(2));
                sum(tau_fblr(1:2))-sum(tau_fblr(3:4))];
            
    FM_w_Motor = [F_motor; M_motor; w_motor];
end

function w_motor_fblr_squared_desired = propeller_motor_FM2w_conversion(F_M_desired, C_T, C_D, L_arm)
    F_M_desired = F_M_desired(:);
    
    inv_C_T = 1/C_T;
    inv_C_D = 1/C_D;
    inv_2_L_C_T = 2/(L_arm*C_T);
    
    Conv_Mat = 0.25*[   -inv_C_T     0               inv_2_L_C_T    inv_C_D;
                        -inv_C_T     0               -inv_2_L_C_T     inv_C_D;
                        -inv_C_T     inv_2_L_C_T     0               -inv_C_D;
                        -inv_C_T     -inv_2_L_C_T    0               -inv_C_D];
    
    w_motor_fblr_squared_desired = Conv_Mat*F_M_desired;
    w_motor_fblr_squared_desired(w_motor_fblr_squared_desired<0) = 0;
end

function w_motor_fblr_desired = quaternion_feedback_and_altitude_control(q_current, ...
    w_current, rv_current, J_inertia, ...
    C_T, C_D, L_arm, C_BR, mass_quadcopter, grv_acce)
    
    zR_desired = -30; %[m]
    zdotR_desired = 0; %[m/s]
    K_qf = 0.01*eye(3);
    C_qf = 0.001*eye(3);
    k1 = 0.1;
    k2 = 0.5;
    
    J_inertia_model = [0.005 -0.001  0.004;
                    -0.001  0.006 -0.002;
                    0.004 -0.002  0.004]; % vehicle moment of inertia [kg m^2]
    
    J_inertia = J_inertia_model;
                
    q_13 = q_current(1:3); q_13 =q_13(:);
    w = w_current(:);
    
    Fmg_R = grv_acce*mass_quadcopter; %[N]
    Falt_R = k1*(zR_desired-rv_current(3))+k2*(zdotR_desired-rv_current(6));
    F_desired_R = [0;0;-Fmg_R+Falt_R];
    F_desired_B = C_BR*F_desired_R;
    
    u_qf = -K_qf*q_13 - C_qf*w - cross(w,J_inertia*w);
    M_Desired = u_qf;
    
    F_M_desired = [F_desired_B(3); M_Desired];
                
    w_motor_fblr_squared_desired = propeller_motor_FM2w_conversion(F_M_desired, ...
        C_T, C_D, L_arm);
   
    w_motor_fblr_desired = sqrt(w_motor_fblr_squared_desired);

end
