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

% number of cells
num_cell = 20;

% simulation time values
time_current = 0;    % initial time 
time_final   = 300; % final time [min]
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
delta_worst = [-1 -1 1 1 -1 1 1 -1 1 1 -1 1 -1 1]';
p_delta = 2;
ki_para=ki_para_org.*(1+(p_delta/100)*delta_worst);

num_reaction = length(ki_para);

% nominal initial number of molecules
ACA_n   = 66535; % [# of molecules]
PKA_n   = 24282; % [# of molecules]
ERK2_n  = 996;   % [# of molecules]
REGA_n  = 25443; % [# of molecules]
icAMP_n = 14638; % [# of molecules]
ecAMP_n = 6365;  % [# of molecules]
CAR1_n  = 21697; % [# of molecules]

% initial conditions for each cell
ACA     = ACA_n*ones(1,num_cell) + randi(ACA_n,1,num_cell);
PKA     = PKA_n*ones(1,num_cell) + randi(PKA_n,1,num_cell);
ERK2    = ERK2_n*ones(1,num_cell) + randi(ERK2_n,1,num_cell);
REGA    = REGA_n*ones(1,num_cell) + randi(REGA_n,1,num_cell);
icAMP   = icAMP_n*ones(1,num_cell) + randi(icAMP_n,1,num_cell);
CAR1    = CAR1_n*ones(1,num_cell) + randi(CAR1_n,1,num_cell);

% total external cAMP
ecAMP_total = num_cell*ecAMP_n + randi(ecAMP_n,1);
% the total external cAMP distributed equally to each cell
ecAMP   = (ecAMP_total/num_cell)*ones(1,num_cell);

% storing data: only store the total external cAMP
species_all = zeros(max_num_data, 2);
species_all(1,:) = [time_current ecAMP_total];
data_idx = 1;

propensity_a = zeros(num_reaction,num_cell);

while data_idx < max_num_data
   
    % each cell
    for idx_cell=1:num_cell 
        propensity_a(1,idx_cell) = ki_para(1)*CAR1(idx_cell);
        propensity_a(2,idx_cell) = ki_para(2)*ACA(idx_cell)*PKA(idx_cell)/(NA*Cell_Vol*1e-6);
        propensity_a(3,idx_cell) = ki_para(3)*icAMP(idx_cell);
        propensity_a(4,idx_cell) = ki_para(4)*PKA(idx_cell);
        propensity_a(5,idx_cell) = ki_para(5)*CAR1(idx_cell);
        propensity_a(6,idx_cell) = ki_para(6)*PKA(idx_cell)*ERK2(idx_cell)/(NA*Cell_Vol*1e-6);
        propensity_a(7,idx_cell) = ki_para(7)*(NA*Cell_Vol*1e-6);
        propensity_a(8,idx_cell) = ki_para(8)*ERK2(idx_cell)*REGA(idx_cell)/(NA*Cell_Vol*1e-6);
        propensity_a(9,idx_cell) = ki_para(9)*ACA(idx_cell);
        propensity_a(10,idx_cell) = ki_para(10)*REGA(idx_cell)*icAMP(idx_cell)/(NA*Cell_Vol*1e-6);
        propensity_a(11,idx_cell) = ki_para(11)*ACA(idx_cell);
        propensity_a(12,idx_cell) = ki_para(12)*ecAMP(idx_cell);
        propensity_a(13,idx_cell) = ki_para(13)*ecAMP(idx_cell);
        propensity_a(14,idx_cell) = ki_para(14)*CAR1(idx_cell);
    end
    
    % determine the reaction time tau
    sum_propensity_a = sum(propensity_a(:));
    tau = exprnd(1/sum_propensity_a);

    % determine the reaction
    normalized_propensity_a = propensity_a(:)/sum_propensity_a;
    cumsum_propensity_a = cumsum(normalized_propensity_a);
    which_reaction = rand(1);
    reaction_idx = cumsum((cumsum_propensity_a-which_reaction)<0);
    active_reaction = reaction_idx(end)+1;
    reaction = rem(active_reaction,num_reaction);
    if reaction==0
       reaction = num_reaction;
    end
    reaction_cell_num = ceil(active_reaction/num_reaction);

    % update number of molecules
    switch reaction
        case 1
            ACA(reaction_cell_num) = ACA(reaction_cell_num) + 1;
        case 2
            ACA(reaction_cell_num) = ACA(reaction_cell_num) - 1;
        case 3
            PKA(reaction_cell_num) = PKA(reaction_cell_num) + 1;
        case 4
            PKA(reaction_cell_num) = PKA(reaction_cell_num) - 1;
        case 5
            ERK2(reaction_cell_num) = ERK2(reaction_cell_num) + 1;
        case 6
            ERK2(reaction_cell_num) = ERK2(reaction_cell_num) - 1;
        case 7
            REGA(reaction_cell_num) = REGA(reaction_cell_num) + 1;
        case 8
            REGA(reaction_cell_num) = REGA(reaction_cell_num) - 1;
        case 9
            icAMP(reaction_cell_num) = icAMP(reaction_cell_num) + 1;
        case 10
            icAMP(reaction_cell_num) = icAMP(reaction_cell_num) - 1;
        case 11
            %ecAMP(reaction_cell_num) = ecAMP(reaction_cell_num) + 1;
            ecAMP_total = ecAMP_total + 1;
        case 12
            %ecAMP(reaction_cell_num) = ecAMP(reaction_cell_num) - 1;
            ecAMP_total = ecAMP_total - 1;
        case 13
            CAR1(reaction_cell_num) = CAR1(reaction_cell_num) + 1;
        case 14
            CAR1(reaction_cell_num) = CAR1(reaction_cell_num) - 1;
        otherwise
            error('Wrong reaction number!');
    end
    
    % distribute the total ecAMP equally to all cells, where allows
    % non-integer numbers
    ecAMP = floor(ecAMP_total/num_cell+0.5)*ones(1,num_cell);
    
    time_current = time_current + tau;
    
    if time_record < time_current
       data_idx = data_idx + 1;
       species_all(data_idx,:) = [time_current ecAMP_total];
       time_record = time_record + dt_record;
       disp(time_record);
    end
    
end

figure;
plot(species_all(:,1),species_all(:,2),'-');
set(gca,'FontSize',14);
xlabel('time [min]');
ylabel('[# of molecules]');
title('total external cAMP');
