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

% generate random quaternion
eig_ax = 2*(rand(3,1)-0.5);
eig_ax = eig_ax/norm(eig_ax); % random axis normalized

eig_ang = pi*2*(rand(1)-0.5); % angle between +/- pi [rad]

% random quaternion
q13 = eig_ax*sin(eig_ang/2);
q4 = cos(eig_ang/2);

% random quaternion to dcm and the dcm to quaternion
quat = [q13(:); q4];
dcm = q2dcm(quat);

r1R = [-0.6794 -0.3237 -0.6586]'; r1R = r1R/norm(r1R);
r2R = [-0.7296  0.5858  0.3528]'; r2R = r2R/norm(r2R);
r3R = [-0.2718  0.6690 -0.6918]'; r3R = r3R/norm(r3R);
r4R = [-0.2062 -0.3986  0.8936]'; r4R = r4R/norm(r4R); 
r5R = [ 0.6858 -0.7274 -0.0238]'; r5R = r5R/norm(r5R);

r1B = dcm*r1R;
r2B = dcm*r2R;
r3B = dcm*r3R;
r4B = dcm*r4R;
r5B = dcm*r5R;

rB_all = [r1B r2B r3B r4B r5B];
rR_all = [r1R r2R r3R r4R r5R];
a_i = ones(size(rB_all,2),1);

q_est = QUEST(rB_all,rR_all,a_i);
dcm_est = q2dcm(q_est);

fprintf('quaternion error = %4.2e\n',norm(quat-q_est))
fprintf('dcm error = %4.2e\n',norm(dcm-dcm_est))

function dcm = q2dcm(quat) 
% be careful matlab aerospace toolbox has quat2dcm function

q13 =quat(1:3); q13 = q13(:);
q4 = quat(4);

q13x = [ 0          -q13(3)       q13(2);
         q13(3)      0           -q13(1);
        -q13(2)      q13(1)       0];

dcm = (q4^2-q13'*q13)*eye(3) + 2*q13*q13' - 2*q4*q13x;
end


function q_est = QUEST(w,v,a_w)

num_obs = size(v,2);

% construct B & z
B = zeros(3);
z = zeros(3,1);
for idx=1:num_obs
    B = B+a_w(idx)*w(:,idx)*v(:,idx)';
    z = z+a_w(idx)*cross(w(:,idx),v(:,idx));
end

S = B + B';
sgm = B(1,1)+B(2,2)+B(3,3);

% trace of adj(S) for S symmetric
kappa =     (S(2,2)*S(3,3)-S(2,3)^2) ...
    +   (S(1,1)*S(3,3)-S(1,3)^2) ...
    +   (S(1,1)*S(2,2)-S(1,2)^2);

delta = det(S);

a = sgm^2 - kappa;
b = sgm^2 + z'*z;
c = delta  +z'*S*z;
d = z'*S*S*z;

%--------------------------------------
% Newton-Raphson Method to find a root
%--------------------------------------
lamda_c = 10; % initial guess is 10
tol_NR = true;
max_num_itr = 500;
cur_itr = 1;

while tol_NR
    lamda_p =  lamda_c;
    
    f_c = lamda_c^4 - (a+b)*lamda_c^2 - c*lamda_c + (a*b+c*sgm-d);
    dfdx = 4*lamda_c^3 - 2*(a+b)*lamda_c - c;
    lamda_c = lamda_c - f_c/dfdx;
    
    if abs(lamda_c - lamda_p) < 1e-6
        tol_NR = false;
    end
    
    cur_itr = cur_itr + 1;
    
    if cur_itr > max_num_itr
        break
    end
    
end
lamda = lamda_c;

Y = (((sgm+lamda)*eye(3)-S))\z;

q_est = 1/sqrt(1+sum(Y.*Y))*[Y; 1];


end
