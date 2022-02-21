clear;

grv_const = 9.81; % [m/s^2]
init_pos = 0.0; %[m]
init_vel = 0.5; % [m/s]
init_mass = 5.0; %[kg]

init_time = 0; % [s]
final_time = 5.0; % [s]
time_interval = [init_time final_time];

x0 = [init_pos init_vel init_mass];
[tout,xout] = ode45(@(time,state) free_falling_obj(state,grv_const,time), time_interval, x0);

figure(1);
plot(tout,xout(:,1))
ylabel('position [m]');
xlabel('time [s]');

figure(2);
plot(tout,xout(:,2))
ylabel('velocity [m/s]');
xlabel('time [s]');

figure(3);
plot(tout,xout(:,3))
ylabel('m(t) [kg]');
xlabel('time [s]');

function dxdt = free_falling_obj(state,grv_const,time)
    x1 = state(1);
    x2 = state(2);
    x3 = state(3);

    dxdt = zeros(3,1);
    dxdt(1) = x2;
    dxdt(2) = grv_const + (x3-2)*(x2/x3);
    dxdt(3) = -x3 + 2;
end