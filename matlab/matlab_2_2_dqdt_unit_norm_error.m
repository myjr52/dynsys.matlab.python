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

init_time = 0; % [s]
final_time = 6000.0; % [s]
N_time = final_time*10;
time_interval = linspace(init_time,final_time,N_time);

q0 = [0 0 0 1]';
RelTol_cases = [1e-3 1e-6 1e-9];
N_sim = length(RelTol_cases);

err_idx = zeros(N_sim,N_time);
for idx=1:N_sim
    ode_options = odeset('RelTol',RelTol_cases(idx),'AbsTol',1e-3*RelTol_cases(idx));
    [tout,qout] = ode45(@(time,state) dqdt_attitude_kinematics(time,state), time_interval, q0, ode_options);
    err_idx(idx,:) = abs(sum(qout.^2,2)-1);
end

figure;
plot(tout,log(err_idx));
xlabel('time [s]');
ylabel('log|{\bf q}^T {\bf q} - 1|');
set(gca,'FontSize',14);
legendStrings = "RelTol = " + string(RelTol_cases);
legend(legendStrings,'Location','southeast');

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
