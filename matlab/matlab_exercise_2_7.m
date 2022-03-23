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

r1R = [-0.6794 -0.3237 -0.6586]';
r2R = [-0.7296  0.5858  0.3528]';
r3R = [-0.2718  0.6690 -0.6918]';
r4R = [-0.2062 -0.3986 0.8936]';
r5R = [0.6858 -0.7274 -0.0238]';

A1 = blkdiag(r1R',r1R',r1R');
A2 = blkdiag(r2R',r2R',r2R');
A3 = blkdiag(r3R',r3R',r3R');
A4 = blkdiag(r4R',r4R',r4R');
A5 = blkdiag(r5R',r5R',r5R');

A = [A1;A2;A3;A4;A5];


r1B = [-0.2147 -0.7985 0.562]';
r2B = [-0.7658 0.4424 0.4667]';
r3B = [-0.8575 -0.4610 -0.228]';
r4B = [0.4442 0.6863 0.5758]';
r5B = [0.9407 -0.1845 -0.2847]';

rB = [r1B; r2B; r3B; r4B; r5B];

C_BR=(reshape((A'*A)\(A'*rB),3,3))'
