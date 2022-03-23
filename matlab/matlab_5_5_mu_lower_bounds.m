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
