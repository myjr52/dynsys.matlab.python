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
final_time = 10.0; % [s]
time_interval = [init_time final_time];

J_inertia = [0.005 -0.001  0.004;
            -0.001  0.006 -0.002;
             0.004 -0.002  0.004]; % vehicle moment of inertia [kg m^2]
J_inv = inv(J_inertia);
 
C_T = 8.8e-7; % motor thruster coefficient [N/(rad/s)^2]
C_D = 11.3e-8; % motor drag coefficient [Nm/(rad/s)^2]
L_arm = 0.127; % length from the centre of quadcopter to the motor [m]

quadcopter_uav = {J_inertia, J_inv, C_T, C_D, L_arm};

q0 = [0 0 0 1]'; % initial quaternion
w0 = [0 0 0]'; % initial angular velocity

state_0 = [q0; w0]; % states including q0 and omega0

% use global variables only for saving values
global global_motor_time_FM_all;
global_motor_time_FM_all = [];

% minimum time interval for saving values to the global
dt_save = 0.05; %[s]

%% simulation
ode_options = odeset('RelTol',1e-6,'AbsTol',1e-9, 'MaxStep', 0.01);
[tout,state_out] = ode45(@(time,state) dqdt_dwdt(time,state,quadcopter_uav,dt_save), ...
    time_interval, state_0, ode_options);

qout = state_out(:,1:4);
wout = state_out(:,5:7);

time_Motor = global_motor_time_FM_all(:,1);
Force_Motor = global_motor_time_FM_all(:,2);
Torque_Motor = global_motor_time_FM_all(:,3:5);

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

subplot(212)
plot(tout,wout(:,1),'r-',tout,wout(:,2),'b--',tout,wout(:,3),'m-.')
set(gca,'FontSize',14);
ylabel('\omega [rad/s]');
xlabel('time [s]');
legend('\omega_1','\omega_2','\omega_3');

figure;
subplot(211);
plot(time_Motor,Force_Motor)
set(gca,'FontSize',14);
ylabel('Force [N]');
xlabel('time [s]');
axis([init_time final_time -15 0]);
subplot(212);
plot(time_Motor,Torque_Motor(:,1),'r-',time_Motor,Torque_Motor(:,2),'b--', ...
    time_Motor,Torque_Motor(:,3),'m-.')
set(gca,'FontSize',14);
ylabel('Torque [Nm]');
xlabel('time [s]');
legend('M_1','M_2','M_3');

%% functions
function dstate_dt = dqdt_dwdt(time,state,quadcopter_uav,dt_save)

    global global_motor_time_FM_all;
    
    q_current = state(1:4);
    q_current = q_current(:)/norm(q_current);
    
    w_current = state(5:7);
    w_current = w_current(:);
    
    J_inertia = quadcopter_uav{1};
    inv_J = quadcopter_uav{2};
    C_T = quadcopter_uav{3};
    C_D = quadcopter_uav{4};
    L_arm = quadcopter_uav{5};
    
    %--------------------------------
    % Begin: this part is controller
    %--------------------------------
    M_Desired = [    0.00001+0.0005*sin(2*time);
                    -0.00002+0.0001*cos(0.75*time);
                    -0.0001]; % [Nm]
                
    mg = 10; %[N]
    F_M_desired = [-mg; M_Desired];
                
    w_motor_fblr_squared_desired = propeller_motor_FM2w_conversion(F_M_desired, ...
        C_T, C_D, L_arm);
    w_motor_fblr_desired = sqrt(w_motor_fblr_squared_desired);
    %--------------------------------
    % End: this part is controller
    %--------------------------------
    
    % Motor Force & Torque
    FM_Motor = propeller_motor_actuator(C_T, C_D, L_arm, w_motor_fblr_desired);
    M_torque = FM_Motor(2:4);
    
    if time < 1e-200
        global_motor_time_FM_all = [time FM_Motor(:)'];
    elseif time > global_motor_time_FM_all(end,1)+dt_save
        global_motor_time_FM_all = [global_motor_time_FM_all; time FM_Motor(:)'];
    end
    
    % Kinematics & Dynamics
    dqdt = dqdt_attitude_kinematics(q_current,w_current);
    dwdt = dwdt_attitude_dynamics(w_current, J_inertia, inv_J, M_torque);
    
    dstate_dt = [dqdt(:); dwdt(:)];

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

function FM_Motor = propeller_motor_actuator(C_T,C_D,L_arm,w_command)

    % assume perfect motor angular velocity control
    w_motor = w_command(:);

    F_fblr = C_T*(w_motor.^2);
    tau_fblr = C_D*(w_motor.^2);
    
    F_motor = -sum(F_fblr);
    M_motor = [ L_arm*(F_fblr(3)-F_fblr(4));
                L_arm*(F_fblr(1)-F_fblr(2));
                sum(tau_fblr(1:2))-sum(tau_fblr(3:4))];
            
    FM_Motor = [F_motor; M_motor];
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
