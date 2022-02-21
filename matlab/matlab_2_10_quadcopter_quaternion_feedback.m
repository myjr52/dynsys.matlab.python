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

quadcopter_uav = {J_inertia, J_inv, C_T, C_D, L_arm};

q0 = (1/sqrt(4))*[1 1 -1 1]'; % initial quaternion
w0 = [0.1 -0.2 0.1]'; % initial angular velocity

r0 = [0 0 -50]';
v0 = [0 0 0]';

state_0 = [q0; w0; r0; v0]; % states including q0 and omega0

% use global variables only for saving values
global global_motor_time_FM_w_all;
global_motor_time_FM_w_all = [];

% minimum time interval for saving values to the global
dt_save = 0.005; %[s]

%% simulation
ode_options = odeset('RelTol',1e-6,'AbsTol',1e-9, 'MaxStep', 0.01);
[tout,state_out] = ode45(@(time,state) dqdt_dwdt_drvdt(time,state,quadcopter_uav,dt_save), ...
    time_interval, state_0, ode_options);

qout = state_out(:,1:4);
wout = state_out(:,5:7);
rout = state_out(:,8:10);
vout = state_out(:,11:13);

time_Motor = global_motor_time_FM_w_all(:,1); %[s]
Force_Motor = global_motor_time_FM_w_all(:,2); %[N]
Torque_Motor = global_motor_time_FM_w_all(:,3:5); %[Nm]
Omega_Motor = global_motor_time_FM_w_all(:,6:9)*(60/(2*pi)); %[rpm]

% clear all global variables
clearvars -global

%% plot
figure;
subplot(211)
plot(tout,qout(:,1),'b-',tout,qout(:,2),'r--',tout,qout(:,3),'g-.',tout,qout(:,4),'m:')
set(gca,'FontSize',14);
ylabel('quaternion');
xlabel('time [s]');
legend('q_1','q_2','q_3','q_4');
axis([init_time final_time -1.1 1.1]);

subplot(212)
plot(tout,wout(:,1),'r-',tout,wout(:,2),'b--',tout,wout(:,3),'m-.')
set(gca,'FontSize',14);
ylabel('\omega [rad/s]');
xlabel('time [s]');
legend('\omega_1','\omega_2','\omega_3');
axis([init_time final_time -3 3]);

figure;
subplot(211);
plot(time_Motor,Force_Motor)
set(gca,'FontSize',14);
ylabel('Force [N]');
xlabel('time [s]');
axis([init_time final_time -10 0]);
subplot(212)
plot(time_Motor,Torque_Motor(:,1),'r-',time_Motor,Torque_Motor(:,2),'b--', ...
    time_Motor,Torque_Motor(:,3),'m-.')
set(gca,'FontSize',14);
ylabel('Torque [Nm]');
xlabel('time [s]');
legend('M_1','M_2','M_3');
axis([init_time final_time -0.02 0.02]);

figure;
plot(time_Motor,Omega_Motor(:,1),'b-',time_Motor,Omega_Motor(:,2),'r--', ...
    time_Motor,Omega_Motor(:,3),'g-.',time_Motor,Omega_Motor(:,4),'m:')
set(gca,'FontSize',14);
ylabel('Motor Speed [rpm]');
xlabel('time [s]');
legend('\omega_{m1}','\omega_{m2}','\omega_{m3}','\omega_{m4}','Location','southeast');

figure;
subplot(211)
plot(tout,-rout(:,3))
set(gca,'FontSize',14);
ylabel('Altitude [m]');
axis([init_time final_time 0 60]);
subplot(212)
plot(tout,vout(:,1),'r-',tout,vout(:,2),'b--',tout,vout(:,3),'m-.')
set(gca,'FontSize',14);
ylabel('Velocity [m/s]');
xlabel('time [s]');
legend('v_1','v_2','v_3');
axis([init_time final_time -20 20]);

%% functions
function dstate_dt = dqdt_dwdt_drvdt(time,state,quadcopter_uav,dt_save)

    m_quadcopter = 0.49; %[kg]
    grv_acce = 9.81; %[m/s^2]

    global global_motor_time_FM_w_all;
    
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
    
    if time < 1e-200
        global_motor_time_FM_w_all = [time FM_w_Motor(:)'];
    elseif time > global_motor_time_FM_w_all(end,1)+dt_save
        global_motor_time_FM_w_all = [global_motor_time_FM_w_all; time FM_w_Motor(:)'];
    end
    
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