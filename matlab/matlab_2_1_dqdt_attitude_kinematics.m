clear;

init_time = 0; % [s]
final_time = 60.0; % [s]
time_interval = [init_time final_time];

q0 = [0 0 0 1]';
[tout,qout] = ode45(@(time,state) dqdt_attitude_kinematics(time,state), ...
    time_interval, q0);

figure;
plot(tout,qout(:,1),'b-',tout,qout(:,2),'r--',tout,qout(:,3),'g-.',tout,qout(:,4),'m:')
ylabel('quaternion');
xlabel('time [s]');
legend('q1','q2','q3','q4');
set(gca,'FontSize',14);

figure(2);
plot(tout,sum(qout.^2,2)-1);
hold on;


function dqdt = dqdt_attitude_kinematics(time,state)
    q_true = state(:);
    
    w_true(1) = 0.1*sin(2*pi*0.005*time); % [rad/s]
    w_true(2) = 0.05*cos(2*pi*0.01*time + 0.2); %[rad/s]
    w_true(3) = 0.02; %[rad/s]
    w_true = w_true(:);
    
    wx = [  0           -w_true(3)  w_true(2);
            w_true(3)   0           -w_true(1);
            -w_true(2)  w_true(1)   0];
        
    Omega = [   -wx         w_true;
                -w_true'    0];
            
    dqdt = 0.5*Omega*q_true;
end