clear;

%----------------------------------------------------------------------
% define symbols
%----------------------------------------------------------------------
% define time interval
syms Dt real;

% define symbols for aircraft's and target's control inputs
syms ux0 uy0 ux1 uy1 real;
syms wx0 wy0 wx1 wy1 real;

% define symbols for the initial conditions
syms xa0 ya0 vxa0 vya0 real;
syms xt0 yt0 real; 

% target characteristics 
syms th w_max real;

%----------------------------------------------------------------------
% Dynamics
%----------------------------------------------------------------------
% aircraft & target dynamics
Fa = eye(4) + [zeros(2) Dt*eye(2); zeros(2,4)];
Ga = [zeros(2); Dt*eye(2)];
Ca = eye(2,4);

Ft = eye(2);
Gt = Dt*eye(2);
Ct = eye(2);

% control inputs
u_vec_0 = [ux0 uy0]';
w_vec_0 = [wx0 wy0]';
u_vec_1 = [ux1 uy1]';
w_vec_1 = [wx1 wy1]';

% initial conditions
xa_vec_0 = [xa0 ya0 vxa0 vya0]';
xt_vec_0 = [xt0 yt0]';

% state propagation
xa_k_plus_1 = Fa*xa_vec_0    + Ga*u_vec_0;
xa_k_plus_2 = Fa*xa_k_plus_1 + Ga*u_vec_1;
y_k_plus_1 = Ca*xa_k_plus_1;
y_k_plus_2 = Ca*xa_k_plus_2;

xt_k_plus_1 = Ft*xt_vec_0    + Gt*w_vec_0;
xt_k_plus_2 = Ft*xt_k_plus_1 + Gt*w_vec_1;
z_k_plus_1 = Ct*xt_k_plus_1;
z_k_plus_2 = Ct*xt_k_plus_2;

%----------------------------------------------------------------------
% calculate the cost function in the original form
%----------------------------------------------------------------------
dyz_1=(y_k_plus_1-z_k_plus_1);
dyz_2=(y_k_plus_2-z_k_plus_2);
J_over_dt_2 = simplify(expand(dyz_1'*dyz_1+dyz_2'*dyz_2));

pretty(collect(J_over_dt_2,[ux0 uy0 ux1 uy1]));

ux_poly = coeffs(J_over_dt_2, ux0);
alpha = ux_poly(2)/ux_poly(3);

uy_poly = coeffs(ux_poly(1), uy0);
beta = uy_poly(2)/uy_poly(3);

gama = uy_poly(1);

poly_recover=(alpha*(Dt^4)*ux0 + Dt^4*ux0^2 + beta*(Dt^4)*uy0 + Dt^4*uy0^2) + gama;

% check alpha, beta, gama are correct: the following must return zero
zero_check = eval(expand(poly_recover-J_over_dt_2));
fprintf('Is this zero? %4.2f \n',zero_check);

