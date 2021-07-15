%Plot 3D bifurcation diagrams varying V_pmTORC2 and IC of AMPK
clc; clear;
tic % start timer

%Import model parameters
param = importdata('modelParameters.txt');

%Set parameters for simulation
Abundance = 0:1000; %Range of initial condition of AMPK to simulate
N = length(Abundance); %Number of points in AMPK
tt = 0:1000; %Timespan
y0 = [50;0;100;0;0;100;0;250;0;200;0;0;0;350;0;250;0;0;250;0]; %Initial condition


%Set points for V_pmTORC2 in the range 1e-6 to 0.02 for simulation
param_1 = linspace(1e-6, 0.001, 20);
param_2 = linspace(0.001, 0.01, 20);
param_3 = linspace(0.01, 0.02, 20);
param_range = [param_1(1:end-1),param_2(1:end-1), param_3];


%Initialize matrices to store results
output_max_val = zeros(N, length(param_range), 2);
output_min_val = zeros(N, length(param_range), 2);
output_sum = zeros(N, length(param_range), 2);
output_product = zeros(N, length(param_range), 2);
output_div = zeros(N, length(param_range), 2);
output_pAMPK = zeros(N, length(param_range), 2);

%Solve the ODE for each parameter and initial condition of AMPK
for i = 1:length(Abundance)
   y0(16) = Abundance(i); %Change IC of AMPK
   parfor j = 1:length(param_range)
       param_copy = param;
       param_copy(26) = param_range(j) * 3600; %Change V_pmTORC2
       
       %Solve the system of ODEs
       [t, answer] = ode23s(@(t, x) dR2(t, x, param_copy), tt, y0); 
       
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
save bifurcation_V_pmTORC2 %Store simulation results