%Plot 3D bifurcation diagrams
%varying AMPK and DEPTOR levels
clc; clear;

%Import model parameters
param = importdata('modelParameters.txt');

%Set parameters for simulation
param(1) = 0.01 * 3600; %Increase V_IR to avoid oscillation
Abundance = 0:5:1000; %Range of initial condition of AMPK and DEPTOR to simulate
tt = 0:1000; %Timespan
y0 = [50;0;100;0;0;100;0;250;0;200;0;0;0;350;0;250;0;0;250;0]; %Initial condition
N = length(Abundance);

%Initialize matrices to store resulting steady state concentrations 
output_max_val = zeros(N, N, 2);
output_min_val = zeros(N, N, 2);
output_sum = zeros(N, N, 2);
output_product = zeros(N, N, 2);
output_div = zeros(N, N, 2);


%Solve the system ODEs at each IC of AMPK and DEPTOR
for j = 1:N
    y0(14) = Abundance(j); %Change IC of DEPTOR
    parfor i = 1:N
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
        
        fprintf('i = %d j = %d\n',i, j); % Print current indices of loop
    end
end

save(append('bifurcation_AMPK_DEPTOR_V_IR_', sprintf('%0.0e',param(1)/3600))) %Store simulation results