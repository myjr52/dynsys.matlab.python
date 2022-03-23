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

clear

kP = 20;
kI = 2.5e-4;
gamma_G = 8e-4;
ks2 = 4e-4;

A_PI = [0 0 0; 0 -gamma_G 0; ks2 ks2 -ks2];
B_PI = [kI; gamma_G*kP; 0];
C_PI = [0 0 1];
D_PI = 0;

sys_PI = ss(A_PI,B_PI,C_PI,D_PI);

A_true = 0;
B_true = 1;
C_true = kI;
D_true = kP;

sys_true_PI = ss(A_true,B_true,C_true,D_true);

% bode plots
freq = logspace(-7,-3,1000); % [rad/time]
[mm1,pp1]=bode(sys_PI,freq);
[mm2,pp2]=bode(sys_true_PI,freq);

figure; clf;
subplot(211);
semilogx(freq,20*log10(squeeze(mm1)));
hold on;
semilogx(freq,20*log10(squeeze(mm2)),'r--');
set(gca,'FontSize',14);
ylabel('Magnitude [dB]');
legend('Approximation','True');
subplot(212);
semilogx(freq,squeeze(pp1));
hold on;
semilogx(freq,squeeze(pp2),'r--');
set(gca,'FontSize',14);
ylabel('Phase [\circ]');
xlabel('Frequency [rad/time]');


% step response and impulse response
time_sim = linspace(0,30000,300000);
[ys1,~,xs1]=step(sys_PI,time_sim);
[ys2,~]=step(sys_true_PI,time_sim);
[yp1,~,xp2]=impulse(sys_PI,time_sim);
[yp2,~]=impulse(sys_true_PI,time_sim);

figure; clf;
subplot(211);
plot(time_sim/60, ys1);
hold on;
plot(time_sim/60, ys2,'r--');
set(gca,'FontSize',14);
ylabel('[a.u.]');
xlabel('time [minutes]')
title('Step Response');
legend('approximated PI','true PI');
subplot(212);
plot(time_sim/60, yp1);
hold on;
plot(time_sim/60, yp2,'r--');
set(gca,'FontSize',14);
legend('approximated PI','true PI');
ylabel('[a.u.]');
xlabel('time [minutes]')
title('Impulse Response');

