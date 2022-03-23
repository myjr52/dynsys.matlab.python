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

alpha = 0.266;

A = -2;
B = [0.1 0.5 0.1/alpha];
C = [1; 1; 0];
D = [0 0 0; 0 0 0; 0 alpha 0];

N_omega = 300;
omega = logspace(-2,3,N_omega);
mu_ub = zeros(N_omega,1);

for idx=1:N_omega
    jw = omega(idx)*sqrt(-1);
    Mjw = C*inv(jw-A)*B+D;
    [U, S, V] = svd(Mjw);

    mu_ub(idx) = max(diag(S));
    
end

semilogx(omega,mu_ub);
axis([1e-2 1e3 0 1.1]);
ylabel('$\bar{\sigma}\, [M(j\omega)]$',fontsize=14, interpreter='latex');
xlabel('$\omega$ [rad/time]',fontsize=14,interpreter='latex');
