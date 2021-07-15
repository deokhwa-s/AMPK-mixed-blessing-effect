%Plot 3D bifurcation diagrams varying V_IR and IC of AMPK
%change IC of DEPTOR if needed
clc ; clear 
tic % start timer

%Set initial parameters and range of bifurcation analysis
tt = 0:1000; %Timespan
N = 1001; %Number of evalution points
Abundance = linspace(0,1000, N);%Re-define range of Abundance
DEPTOR_IC = 1000; %Modify initial condition of DEPTOR
y0 = [50;0;100;0;0;100;0;250;0;200;0;0;0;DEPTOR_IC;0;250;0;0;250;0]; %Initial condition

%Import model parameters
param = importdata('modelParameters.txt');

%Set points for V_IR in the range 1e-6 to 0.2 for simulation
param_1 = linspace(1e-6, 0.01, 20);
param_2 = linspace(0.01, 0.1, 20);
param_3 = linspace(0.1, 0.2, 20);
param_range = [param_1(1:end-1),param_2(1:end-1), param_3];

%Initialize matrices to store resulting steady state concentrations 
output_max_val = zeros(N, length(param_range), 2);
output_min_val = zeros(N, length(param_range), 2);
output_sum = zeros(N, length(param_range), 2);
output_product = zeros(N, length(param_range), 2);
output_div = zeros(N, length(param_range), 2);
output_pAMPK = zeros(N, length(param_range), 2);

%Solve the system ODEs at each V_IR and IC of [AMPK]
for j = 1:length(param_range)
    param(1) = param_range(j) * 3600; %Change V_IR convert units from per sec -> per h
    parfor i = 1:length(Abundance)
        y0_copy = y0; %Copy IC for parallel computing
        y0_copy(16) = Abundance(i); %Change IC of AMPK
        [t, answer] = ode23s(@(t, x) dR2(t, x, param), tt, y0_copy); %Solve the system of ODEs  

        %Find the max/min, sum, product, quotient of pmTORC1 and pmTORC2
        output_max_val(i, j, :) = max(answer(end-200:end, [9, 11])); %Max of steady state
        output_min_val(i, j, :) = min(answer(end-200:end, [9, 11])); %Min of steady state
		
        %Max/min of steady state for pmTORC1 + pmTORC2
        output_sum(i, j, :) = [max(answer(end-200:end, 9) + answer(end-200:end, 11)); 
                               min(answer(end-200:end, 9) + answer(end-200:end, 11))];
		
        %Max/min of steady state for pmTORC1 * pmTORC2                   
        output_product(i, j, :) = [max(answer(end-200:end, 9) .* answer(end-200:end, 11));		
                                  min(answer(end-200:end, 9) .* answer(end-200:end, 11))];		
		
        %Max/min of steady state for pmTORC1 / pmTORC2                   
        output_div(i, j, :) = [max(answer(end-200:end, 9) ./ answer(end-200:end, 11));		
                              min(answer(end-200:end, 9) ./ answer(end-200:end, 11))];		
        
        %Max/min of steady state for pAMPK                  
        output_pAMPK(i, j, :) = max(answer(end-200:end, 17)); %Max of steady state
        output_pAMPK(i, j, :) = min(answer(end-200:end, 17)); %Min of steady state
        
        fprintf('i = %d j = %d\n',i, j); % Print current indices of loop
    end
end


elapsedTime = toc % end timer
save(append('bifurcation_AMPK_V_IR_3D_DEPTOR_', num2str(DEPTOR_IC)))%Save simulation results