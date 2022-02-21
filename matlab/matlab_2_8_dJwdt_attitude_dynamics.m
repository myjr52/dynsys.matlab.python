clear;

%% simulation parameters
init_time = 0; % [s]
final_time = 10.0; % [s]
time_interval = [init_time final_time];

J_inertia = [0.005 -0.001  0.004;
            -0.001  0.006 -0.002;
             0.004 -0.002  0.004]; % vehicle moment of inertia [kg m^2]
J_inv = inv(J_inertia);
 
J_inv_J_inertia = [J_inertia; J_inv];

q0 = [0 0 0 1]'; % initial quaternion
w0 = [0 0 0]'; % initial angular velocity

state_0 = [q0; w0]; % states including q0 and omega0

%% simulation
ode_options = odeset('RelTol',1e-6,'AbsTol',1e-9, 'MaxStep', 0.01);
[tout,state_out] = ode45(@(time,state) dqdt_dwdt(time,state,J_inv_J_inertia), ...
    time_interval, state_0, ode_options);

qout = state_out(:,1:4);
wout = state_out(:,5:7);

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

%% functions
function dstate_dt = dqdt_dwdt(time,state,J_inv_J_inertia)

    q_current = state(1:4);
    q_current = q_current(:)/norm(q_current);
    
    w_current = state(5:7);
    w_current = w_current(:);
    
    J_inertia = J_inv_J_inertia(1:3,:);
    inv_J = J_inv_J_inertia(4:6,:);
    
    M_torque = [    0.00001+0.0005*sin(2*time);
                    -0.00002+0.0001*cos(0.1*time);
                    -0.0001]; % [Nm]
    
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