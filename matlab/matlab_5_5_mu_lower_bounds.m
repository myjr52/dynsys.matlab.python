clear;

Ns = 5000;
eps = 1e-6;

d_lb = 1e-3;
d_ub = 10;

d_ulb = d_ub - d_lb;

omega = 0;
Mjw = 1/(omega*1j + 2);

num_delta = 2;

while d_ulb > eps
    
    d = (d_lb+d_ub)/2;
    
    delta_1 = rand(Ns,1)*d-d/2;
    delta_2 = rand(Ns,1)*d-d/2;
    
    rand_face = randi(num_delta,Ns,1);
    delta_1(rand_face==1) = d/2;
    delta_2(rand_face==2) = d/2;
   
    Delta = A_delta(delta_1,delta_2)-A_delta(0,0);
    
    if length(unique(sign(real(1-Mjw*Delta)))) == 2
        d_ub = d;
    else
        d_lb = d;
    end
        
    d_ulb = d_ub - d_lb;
    
end

mu_lb = 2/d_ub;

function A_delta_val = A_delta(delta_1, delta_2)
    
    A_delta_val = -2+delta_1+sin(delta_1.*delta_2);

end