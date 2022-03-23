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

% assume random quaternion and the dcm are the satellite attitude
% with respect to the reference
quat_BR = [q13(:); q4];
dcm_BR = q2dcm(quat_BR);

% dcm from body frame to sensor frame
dcm_SB = [  1   0   0
            0   -1  0
            0   0   -1];

% generate random stars in reference frame
total_num_star = 200;        
rR_all = 2*rand(3,total_num_star)-1;
rR_all = rR_all./kron(ones(3,1),sqrt(sum(rR_all.^2))); % normalize it

rB_all = dcm_BR*rR_all;
rS_all = dcm_SB*rB_all;

% star sensor direction in the sensor frame
star_sensor_in_S = [0 0 1]';
star_sensor_fov = 12*pi/180; % 12 degrees in [radian]

% angles between the starts and the star sensor
th_all = acos(rS_all'*star_sensor_in_S);
star_seen = th_all < star_sensor_fov;
rB_seen = rB_all(:,star_seen);

if size(rB_seen,2) > 1
    rR_seen = rR_all(:,star_seen);
    
    if size(rB_seen,2) == 2
        rB_3 = cross(rB_seen(:,1),rB_seen(:,2));
        rB_seen = [rB_seen rB_3(:)];
        
        
        rR_3 = cross(rR_seen(:,1),rR_seen(:,2));
        rR_seen = [rR_seen rR_3(:)];
    end
    
    a_i = ones(size(rB_seen,2),1);
    q_est = QUEST(rB_seen,rR_seen,a_i);
    dcm_est = q2dcm(q_est);
    
    fprintf('quaternion error = %4.2e\n',norm(quat_BR-q_est));
    fprintf('dcm error = %4.2e\n',norm(dcm_BR-dcm_est));
    
else
    fprintf('No enough number of stars are seen\n')
end

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
