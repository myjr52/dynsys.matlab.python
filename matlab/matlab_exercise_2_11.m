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

syms q1 q2 q3 q4 real;
syms w1 w2 w3 real;


q = [q1 q2 q3 q4]';
qinv = [-q1 -q2 -q3 q4]';

w = [w1 w2 w3]';

wx = [0 -w3 w2;
      w3 0 -w1;
      -w2 w1 0];
Omega = [-wx w; 
         -w' 0];
     
dqdt = (1/2)*Omega*q;

% obtain dqinv/dt: Equation (4) in the solution
qinv_dot = -1*simplify(expand(quat_x_quat(qinv,quat_x_quat(dqdt, qinv))));

syms p1 p2 p3 p4 real
p = [p1 p2 p3 p4]';

% prove dqdt (x) p = (1/2) Omega (q (x) p) : Equation (6) in the solution
simplify(expand(quat_x_quat(dqdt,p)-(1/2)*Omega*quat_x_quat(q,p)))

% prove p (x) (q (x) r) = (p (x) q) (x) r: Equation (7) in the solution
syms r1 r2 r3 r4 real
r = [r1 r2 r3 r4]';
simplify(expand(quat_x_quat(p,quat_x_quat(q,r)) - quat_x_quat(quat_x_quat(p,q),r)))

function qxq = quat_x_quat(qp, q)


qp_v = qp(1:3); qp_v = qp_v(:); qp4 = qp(4);
q_v = q(1:3); q_v = q_v(:); q4 = q(4);

qxq_v = qp4*q_v  + qp_v*q4 - cross(qp_v,q_v);
qxq_4 = qp4*q4 - qp_v'*q_v;

qxq = [qxq_v(:); qxq_4];

end
