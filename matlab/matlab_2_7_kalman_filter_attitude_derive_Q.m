clear;

% symbols for omega, noise and the variances
syms w1 w2 w3 Dt nv1 nv2 nv3 nu1 nu2 nu3 real;
syms sgm2_u sgm2_v real; % these are variance, i.e., sigma-squared

wx=[ 0 -w3 w2; w3 0 -w1; -w2 w1 0];
nv=[nv1;nv2;nv3];
nu=[nu1;nu2;nu3];
wc=[nv;nu];

% F & G Matrices
F = [-wx -eye(3); zeros(3,6)];
G = [-eye(3) zeros(3); zeros(3) eye(3)];

% e^{Ft}
Phi = eye(6) + F*Dt + (1/2)*(F^2)*Dt^2 + (1/6)*(F^3)*Dt^3 + (1/24)*(F^4)*Dt^4;

% wd before integral
wd = Phi*wc;

% E(wd wd^T)
cov_wd = simplify(expand(wd*wd'));
Q_cov = sym(zeros(6));

eqn2=sgm2_u==nu1^2;
eqn3=sgm2_u==nu2^2;
eqn4=sgm2_u==nu3^2;
eqn5=sgm2_v==nv1^2;
eqn6=sgm2_v==nv2^2;
eqn7=sgm2_v==nv3^2;

syms q11 q12 q13 q21 q22 q23 q31 q32 q33 real;

eqn_1=q11==expand(cov_wd(1,1));
PPT_11 = subs(eqn_1,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_11 = subs(rhs(PPT_11),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_2=q12==expand(cov_wd(1,2));
PPT_12 = subs(eqn_2,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_12 = subs(rhs(PPT_12),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_3=q13==expand(cov_wd(1,3));
PPT_13 = subs(eqn_3,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_13 = subs(rhs(PPT_13),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_4=q21==expand(cov_wd(2,1));
PPT_21 = subs(eqn_4,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_21 = subs(rhs(PPT_21),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_5=q22==expand(cov_wd(2,2));
PPT_22 = subs(eqn_5,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_22 = subs(rhs(PPT_22),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_6=q23==expand(cov_wd(2,3));
PPT_23 = subs(eqn_6,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_23 = subs(rhs(PPT_23),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0}); 

eqn_7=q31==expand(cov_wd(3,1));
PPT_31 = subs(eqn_7,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_31 = subs(rhs(PPT_31),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_8=q32==expand(cov_wd(3,2));
PPT_32 = subs(eqn_8,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_32 = subs(rhs(PPT_32),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_9=q33==expand(cov_wd(3,3));
PPT_33 = subs(eqn_9,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_33 = subs(rhs(PPT_33),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0}); 

%-------------------------------------
% integral from t_k to t_k + Delta t
%-------------------------------------
eqn_1=q11==expand(int(PPT_11,Dt,[0 Dt]));
eqn_2=q12==expand(int(PPT_12,Dt,[0 Dt]));
eqn_3=q13==expand(int(PPT_13,Dt,[0 Dt]));
eqn_4=q21==expand(int(PPT_21,Dt,[0 Dt]));
eqn_5=q22==expand(int(PPT_22,Dt,[0 Dt]));
eqn_6=q23==expand(int(PPT_23,Dt,[0 Dt]));
eqn_7=q31==expand(int(PPT_31,Dt,[0 Dt]));
eqn_8=q32==expand(int(PPT_32,Dt,[0 Dt]));
eqn_9=q33==expand(int(PPT_33,Dt,[0 Dt]));

% ignore higher-order terms
Q_cov(1,1) = subs(rhs(eqn_1),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0}); 
Q_cov(1,2) = subs(rhs(eqn_2),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(1,3) = subs(rhs(eqn_3),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(2,1) = subs(rhs(eqn_4),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(2,2) = subs(rhs(eqn_5),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(2,3) = subs(rhs(eqn_6),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(3,1) = subs(rhs(eqn_7),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(3,2) = subs(rhs(eqn_8),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(3,3) = subs(rhs(eqn_9),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});

%------------------------------
syms q14 q15 q16 q24 q25 q26 q34 q35 q36 real;
eqn_1=q14==expand(cov_wd(1,4));
PPT_14 = subs(eqn_1,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_14 = subs(rhs(PPT_14),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_2=q15==expand(cov_wd(1,5));
PPT_15 = subs(eqn_2,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_15 = subs(rhs(PPT_15),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_3=q16==expand(cov_wd(1,6));
PPT_16 = subs(eqn_3,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_16 = subs(rhs(PPT_16),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_4=q24==expand(cov_wd(2,4));
PPT_24 = subs(eqn_4,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_24 = subs(rhs(PPT_24),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_5=q25==expand(cov_wd(2,5));
PPT_25 = subs(eqn_5,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_25 = subs(rhs(PPT_25),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_6=q26==expand(cov_wd(2,6));
PPT_26 = subs(eqn_6,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_26 = subs(rhs(PPT_26),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0}); 

eqn_7=q34==expand(cov_wd(3,4));
PPT_34 = subs(eqn_7,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_34 = subs(rhs(PPT_34),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_8=q35==expand(cov_wd(3,5));
PPT_35 = subs(eqn_8,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_35 = subs(rhs(PPT_35),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_9=q36==expand(cov_wd(3,6));
PPT_36 = subs(eqn_9,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_36 = subs(rhs(PPT_36),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0}); 


%-------------------------------------
% integral from t_k to t_k + Delta t
%-------------------------------------
eqn_1=q14==expand(int(PPT_14,Dt,[0 Dt]));
eqn_2=q15==expand(int(PPT_15,Dt,[0 Dt]));
eqn_3=q16==expand(int(PPT_16,Dt,[0 Dt]));
eqn_4=q24==expand(int(PPT_24,Dt,[0 Dt]));
eqn_5=q25==expand(int(PPT_25,Dt,[0 Dt]));
eqn_6=q26==expand(int(PPT_26,Dt,[0 Dt]));
eqn_7=q34==expand(int(PPT_34,Dt,[0 Dt]));
eqn_8=q35==expand(int(PPT_35,Dt,[0 Dt]));
eqn_9=q36==expand(int(PPT_36,Dt,[0 Dt]));

% ignore higher-order terms
Q_cov(1,4) = subs(rhs(eqn_1),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0}); 
Q_cov(1,5) = subs(rhs(eqn_2),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(1,6) = subs(rhs(eqn_3),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(2,4) = subs(rhs(eqn_4),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(2,5) = subs(rhs(eqn_5),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(2,6) = subs(rhs(eqn_6),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(3,4) = subs(rhs(eqn_7),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(3,5) = subs(rhs(eqn_8),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(3,6) = subs(rhs(eqn_9),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});

Q_cov(4,1) = Q_cov(1,4);
Q_cov(5,1) = Q_cov(1,5);
Q_cov(6,1) = Q_cov(1,6);
Q_cov(4,2) = Q_cov(2,4);
Q_cov(5,2) = Q_cov(2,5);
Q_cov(6,2) = Q_cov(2,6);
Q_cov(4,3) = Q_cov(3,4);
Q_cov(5,3) = Q_cov(3,5);
Q_cov(6,3) = Q_cov(3,6);


%------------------------------
syms q44 q45 q46 q54 q55 q56 q64 q65 q66 real;
eqn_1=q44==expand(cov_wd(4,4));
PPT_44 = subs(eqn_1,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_44 = subs(rhs(PPT_44),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_2=q45==expand(cov_wd(4,5));
PPT_45 = subs(eqn_2,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_45 = subs(rhs(PPT_45),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_3=q46==expand(cov_wd(4,6));
PPT_46 = subs(eqn_3,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_46 = subs(rhs(PPT_46),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_4=q54==expand(cov_wd(5,4));
PPT_54 = subs(eqn_4,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_54 = subs(rhs(PPT_54),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_5=q55==expand(cov_wd(5,5));
PPT_55 = subs(eqn_5,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_55 = subs(rhs(PPT_55),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_6=q56==expand(cov_wd(5,6));
PPT_56 = subs(eqn_6,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_56 = subs(rhs(PPT_56),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0}); 

eqn_7=q64==expand(cov_wd(6,4));
PPT_64 = subs(eqn_7,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_64 = subs(rhs(PPT_64),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_8=q65==expand(cov_wd(6,5));
PPT_65 = subs(eqn_8,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_65 = subs(rhs(PPT_65),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0});

eqn_9=q66==expand(cov_wd(6,6));
PPT_66 = subs(eqn_9,{rhs(eqn2),rhs(eqn3),rhs(eqn4),rhs(eqn5),rhs(eqn6),rhs(eqn7)},{lhs(eqn2),lhs(eqn3),lhs(eqn4),lhs(eqn5),lhs(eqn6),lhs(eqn7)});
PPT_66 = subs(rhs(PPT_66),{nv1,nv2,nv3,nu1,nu2,nu3},{0,0,0,0,0,0}); 

%-------------------------------------
% integral from t_k to t_k + Delta t
%-------------------------------------
eqn_1=q44==expand(int(PPT_44,Dt,[0 Dt]));
eqn_2=q45==expand(int(PPT_45,Dt,[0 Dt]));
eqn_3=q46==expand(int(PPT_46,Dt,[0 Dt]));
eqn_4=q54==expand(int(PPT_54,Dt,[0 Dt]));
eqn_5=q55==expand(int(PPT_55,Dt,[0 Dt]));
eqn_6=q56==expand(int(PPT_56,Dt,[0 Dt]));
eqn_7=q64==expand(int(PPT_64,Dt,[0 Dt]));
eqn_8=q65==expand(int(PPT_65,Dt,[0 Dt]));
eqn_9=q66==expand(int(PPT_66,Dt,[0 Dt]));

% ignore higher-order terms
Q_cov(4,4) = subs(rhs(eqn_1),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0}); 
Q_cov(4,5) = subs(rhs(eqn_2),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(4,6) = subs(rhs(eqn_3),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(5,4) = subs(rhs(eqn_4),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(5,5) = subs(rhs(eqn_5),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(5,6) = subs(rhs(eqn_6),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(6,4) = subs(rhs(eqn_7),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(6,5) = subs(rhs(eqn_8),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});
Q_cov(6,6) = subs(rhs(eqn_9),{Dt^9,Dt^8,Dt^7,Dt^6,Dt^5,Dt^4},{0,0,0,0,0,0});











