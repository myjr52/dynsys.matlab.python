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