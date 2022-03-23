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

% simulation time values
time_current = 0;    % initial time 
time_final   = 120.0; % final time [min]
time_record  = time_current; % data record time
dt_record    = 0.1; % minimum time interval for data recording
max_num_data = floor((time_final-time_current)/dt_record+0.5);

% kinetic parameters for the Laub-Loomis Dicty cAMP oscillation 
% network model from k1 to k14
ki_para_org = [2.0; 0.9; 2.5; 1.5; 0.6; 0.8; 1.0; 1.3; 0.3; 0.8; 0.7; 4.9; 23.0; 4.5];
Cell_Vol = 3.672e-14; % [litre]
NA = 6.022e23;        % Avogadro's number
num_molecule_species = 7;

% robustness
%delta_worst = [-1 -1 1 1 -1 1 1 -1 1 1 -1 1 -1 1]';
delta_worst = [1    -1     1     1    -1     1    -1     1    -1     1    -1     1    -1     1]';
%delta_worst = sign(2*rand(14,1)-1);
p_delta = 20;
ki_para=ki_para_org.*(1+(p_delta/100)*delta_worst);

% initial number of molecules
ACA   = 35403; % [# of molecules]
PKA   = 32888; % [# of molecules]
ERK2  = 11838; % [# of molecules]
REGA  = 27348; % [# of molecules]
icAMP = 15489; % [# of molecules]
ecAMP = 4980;  % [# of molecules]
CAR1  = 25423; % [# of molecules]

ACA   = 66535; % [# of molecules]
PKA   = 24282; % [# of molecules]
ERK2  = 996; % [# of molecules]
REGA  = 25443; % [# of molecules]
icAMP = 14638; % [# of molecules]
ecAMP = 6365;  % [# of molecules]
CAR1  = 21697; % [# of molecules]

% storing data
species_all = zeros(max_num_data, num_molecule_species+1);
species_all(1,:) = [time_current ACA PKA ERK2 REGA icAMP ecAMP CAR1];
data_idx = 1;

while data_idx < max_num_data
   
    propensity_a(1) = ki_para(1)*CAR1;
    propensity_a(2) = ki_para(2)*ACA*PKA/(NA*Cell_Vol*1e-6);
    propensity_a(3) = ki_para(3)*icAMP;
    propensity_a(4) = ki_para(4)*PKA;
    propensity_a(5) = ki_para(5)*CAR1;
    propensity_a(6) = ki_para(6)*PKA*ERK2/(NA*Cell_Vol*1e-6);
    propensity_a(7) = ki_para(7)*(NA*Cell_Vol*1e-6);
    propensity_a(8) = ki_para(8)*ERK2*REGA/(NA*Cell_Vol*1e-6);
    propensity_a(9) = ki_para(9)*ACA;
    propensity_a(10) = ki_para(10)*REGA*icAMP/(NA*Cell_Vol*1e-6);
    propensity_a(11) = ki_para(11)*ACA;
    propensity_a(12) = ki_para(12)*ecAMP;
    propensity_a(13) = ki_para(13)*ecAMP;
    propensity_a(14) = ki_para(14)*CAR1;
    
    % determine the reaction time tau
    sum_propensity_a = sum(propensity_a);
    tau = exprnd(1/sum_propensity_a);

    % determine the reaction
    normalized_propensity_a = propensity_a/sum_propensity_a;
    cumsum_propensity_a = cumsum(normalized_propensity_a);
    which_reaction = rand(1);
    reaction_idx = cumsum((cumsum_propensity_a-which_reaction)<0);
    reaction = reaction_idx(end)+1;

    % update number of molecules
    switch reaction
        case 1
            ACA = ACA + 1;
        case 2
            ACA = ACA - 1;
        case 3
            PKA = PKA + 1;
        case 4
            PKA = PKA - 1;
        case 5
            ERK2 = ERK2 + 1;
        case 6
            ERK2 = ERK2 - 1;
        case 7
            REGA = REGA + 1;
        case 8
            REGA = REGA - 1;
        case 9
            icAMP = icAMP + 1;
        case 10
            icAMP = icAMP - 1;
        case 11
            ecAMP = ecAMP + 1;
        case 12
            ecAMP = ecAMP - 1;
        case 13
            CAR1 = CAR1 + 1;
        case 14
            CAR1 = CAR1 - 1;
        otherwise
            error('Wrong reaction number!');
    end
    
    time_current = time_current + tau;
    
    if time_record < time_current
       data_idx = data_idx + 1;
       species_all(data_idx,:) = [time_current ACA PKA ERK2 REGA icAMP ecAMP CAR1];
       time_record = time_record + dt_record;
       disp(time_record);
    end
    
end

figure;
plot(species_all(:,1),species_all(:,6),'-');
hold on;
plot(species_all(:,1),species_all(:,8),'-.');
set(gca,'FontSize',14);
xlabel('time [min]');
ylabel('[# of molecules]');
legend('i-cAMP','CAR1');
