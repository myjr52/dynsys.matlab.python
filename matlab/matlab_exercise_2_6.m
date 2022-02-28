clear;

% generate random quaternion
eig_ax = 2*(rand(3,1)-0.5);
eig_ax = eig_ax/norm(eig_ax); % random axis normalized

eig_ang = pi*2*(rand(1)-0.5); % angle between +/- pi [rad]

% random quaternion
q13 = eig_ax*sin(eig_ang/2);
q4 = cos(eig_ang/2);

% random quaternion to dcm and the dcm to quaternion
dcm = q2dcm([q13(:); q4]);
quat = dcm2q(dcm);

% compare quaternion
norm([q13; q4] - quat(:))

% dcm from the converted quaternion
dcm2 = q2dcm(quat);

% compare dcm
norm(dcm-dcm2)

function dcm = q2dcm(quat) 
% be careful matlab aerospace toolbox has quat2dcm function

q13 =quat(1:3); q13 = q13(:);
q4 = quat(4);

q13x = [ 0          -q13(3)       q13(2);
         q13(3)      0           -q13(1);
        -q13(2)      q13(1)       0];

dcm = (q4^2-q13'*q13)*eye(3) + 2*q13*q13' - 2*q4*q13x;
end


function quat = dcm2q(dcm) 
% be careful matlab aerospace toolbox has dcm2quat function

quat = zeros(4,1);

a1 = (1 + dcm(1,1) - dcm(2,2) - dcm(3,3))/4;
a2 = (1 + dcm(2,2) - dcm(1,1) - dcm(3,3))/4;
a3 = (1 + dcm(3,3) - dcm(1,1) - dcm(2,2))/4;
a4 = (1 + dcm(1,1) + dcm(2,2) + dcm(3,3))/4;

[a_max, a_idx] = max([a1 a2 a3 a4]);
quat(a_idx) = sqrt(a_max);

switch a_idx
    case 1
        quat(2) = (dcm(1,2)+dcm(2,1))/(4*quat(1));
        quat(3) = (dcm(1,3)+dcm(3,1))/(4*quat(1));
        quat(4) = (dcm(2,3)-dcm(3,2))/(4*quat(1));
    case 2
        quat(1) = (dcm(1,2)+dcm(2,1))/(4*quat(2));
        quat(3) = (dcm(2,3)+dcm(3,2))/(4*quat(2));
        quat(4) = (dcm(3,1)-dcm(1,3))/(4*quat(2));
    case 3
        quat(1) = (dcm(1,3)+dcm(3,1))/(4*quat(3));
        quat(2) = (dcm(2,3)+dcm(3,2))/(4*quat(3));
        quat(4) = (dcm(1,2)-dcm(2,1))/(4*quat(3));
    case 4
        quat(1) = (dcm(2,3)-dcm(3,2))/(4*quat(4));
        quat(2) = (dcm(3,1)-dcm(1,3))/(4*quat(4));
        quat(3) = (dcm(1,2)-dcm(2,1))/(4*quat(4));
end

if quat(4) < 0
    quat = - quat;
end
        
end