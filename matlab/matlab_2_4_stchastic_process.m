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

% numer of time samplng & number of stochastic process trial
N_sample = 100;
N_realize = 500;

% time 
dt = 0.1; % [seconds]
time_init = 0;
time_final = dt*N_sample;
time = linspace(time_init,time_final,N_sample);

% declare memory space for x_rand_all to include all trials
x_rand_all = zeros(N_realize,N_sample);

% time varying mean and sqrt(variance) at the time instance
mu_all = linspace(-2,2,N_sample);
sigma_all = linspace(0.1,1.5,N_sample);

% for a fixed time instance, generate the random numbers
% with the mean and the variance at the fixed time
for idx=1:N_sample
    mu_t = mu_all(idx);
    sigma_t = sigma_all(idx);
    
    x_rand = mu_t+sigma_t*randn(N_realize,1);
    x_rand_all(:,idx) = x_rand;
end

% plot all trials with respect to the time

% Warning: this part is only executed with the small N_trial, e.g., 5
% the plot takes really long and causing the computer crashed with the
% large N_trial, e.g., 1000

if N_realize < 10
    figure;
    plot(time,x_rand_all,'k-');
    set(gca,'FontSize',14);
    xlabel('time [s]');
    ylabel('x(t)');
end

% approximate mean and variance from the realisation
% and compare with the true
mu_approx = mean(x_rand_all);
sigma2_approx = var(x_rand_all);
figure;
subplot(211);
plot(time,mu_all); 
hold on;
plot(time,mu_approx,'r--');
set(gca,'FontSize',14);
ylabel('$\mu(t)$','Interpreter','latex');
legend('True','Estimated','Location','southeast');
subplot(212);
plot(time,sigma_all.^2);
hold on;
plot(time,sigma2_approx,'r--');
set(gca,'FontSize',14);
ylabel('$[\sigma(t)]^2$','Interpreter','latex');
xlabel('time [s]');
legend('True','Estimated','Location','southeast');

% esimate the pdf for each instance using N-trials at each instance
N_bin = 100;
x_bin = linspace(-5,5,N_bin);
dx=mean(diff(x_bin));
px_all = zeros(N_bin-1,N_sample);
for jdx=1:N_sample
    x_rand = x_rand_all(:,jdx);
    N_occur = histcounts(x_rand,x_bin);
    px_at_t = N_occur/(dx*N_realize);
    px_all(:,jdx) = px_at_t(:);
end

% plot the estimated pdf
figure;
surf(time,x_bin(1:end-1)+0.5*dx,px_all);
set(gca,'FontSize',14);
xlabel('time [s]','Interpreter','latex');
ylabel('sampling space $x$','Interpreter','latex');
zlabel('$\hat{p}(x)$','Interpreter','latex');
