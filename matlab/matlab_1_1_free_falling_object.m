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
