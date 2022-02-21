clear;

N_alpha   = 100;
alpha_all = linspace(0.05,0.4,N_alpha);
max_mu    = zeros(N_alpha,1);

for gdx = 1:N_alpha
    
    alpha = alpha_all(gdx);

    A = -2;
    B = [0.1 0.5 0.1/alpha];
    C = [1; 1; 0];
    D = [0 0 0; 0 0 0; 0 alpha 0];
    
    N_omega = 300;
    omega = logspace(-2,3,N_omega);
    mu_ub = zeros(N_omega,1);
    
    for idx = 1:N_omega
        jw = omega(idx)*sqrt(-1);
        Mjw = C*inv(jw-A)*B+D;
        [U,S,V] = svd(Mjw);
    
        mu_ub(idx) = max(diag(S));
    end
    
    max_mu(gdx) = max(mu_ub);
end

plot(alpha_all, max_mu);
axis([alpha_all(1),alpha_all(end),0.4,1.6]);
ylabel('$\max \{\bar{\sigma}\, [M(j\omega)]\}$',fontsize=14, interpreter='latex');
xlabel('$\alpha$',fontsize=14,interpreter='latex');