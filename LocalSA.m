%Local Sensitivity Analysis for the mTORC1/AMPK/ULK1 Pathway
clc; clear;            
tic %start timing

%Import model parameters
param = importdata('modelParameters.txt');

%Define parameters names for constructing the heatmap
parameter_names = {'V_IR'; 'Km_IR'; 'V_pIR'; 'Km_pIR'; 'K_IRS_by_pIR'; 
    'Km_IRS_by_pIR'; 'V_pIRS'; 'Km_pIRS'; 'K_AKT_by_pIRS'; 'Km_AKT_by_pIRS'; 
    'K_AKT_by_pmTORC2'; 'Km_AKT_by_pmTORC2'; 'V_pAKT'; 'Km_pAKT'; 'K_mTORC1_by_pAKT'; 
    'Km_mTORC1_by_pAKT'; 'K_pmTORC1'; 'K_pmTORC1_by_pAMPK'; 'Km_pmTORC1_by_pAMPK'; 
    'K_pmTORC1_by_pULK1'; 'Km_pmTORC1_by_pULK1'; 'K_mTORC2_by_pIRS'; 'Km_mTORC2_by_pIRS'; 
    'K_mTORC2_by_pAMPK'; 'Km_mTORC2_by_pAMPK'; 'V_pmTORC2'; 'Km_pmTORC2'; 
    'K_DEPTOR_by_pmTORC1'; 'Km_DEPTOR_by_pmTORC1'; 'K_DEPTOR_by_pmTORC2'; 
    'Km_DEPTOR_by_pmTORC2'; 'V_pDEPTOR'; 'Km_pDEPTOR'; 'K_mTORC1_DEPTOR_form'; 
    'K_mTORC1_DEPTOR_diss'; 'K_mTORC2_DEPTOR_form'; 'K_mTORC2_DEPTOR_diss'; 
    'K_IRS_to_iIRS'; 'Km_IRS_to_iIRS'; 'V_iIRS'; 'Km_iIRS'; 'K_AMPK'; 'K_AMPK_by_SIRT1'; 
    'Km_AMPK'; 'K_pAMPK'; 'K_pAMPK_by_pULK1'; 'K_pAMPK_by_pmTORC1'; 'Km_pAMPK'; 'K_SIRT1'; 
    'K_SIRT1_by_pAMPK'; 'Km_SIRT1'; 'K_SIRT1_diss'; 'K_ULK1'; 'K_ULK1_by_pAMPK'; 'Km_ULK1';
    'K_pULK1'; 'K_pULK1_by_pmTORC1'; 'Km_pULK1'}; 

parameter_names = strrep(parameter_names, '_', '\_'); % Replace '_' by '\_'

%Define variable names  
variable_names = {
    'IR'; 'pIR'; 'IRS'; 'pIRS'; 'iIRS'; 'AKT'; 
    'pAKT'; 'mTORC1'; 'pmTORC1'; 'mTORC2'; 'pmTORC2'; 
    'mTORC1-DEPTOR'; 'mTORC2-DEPTOR'; 'DEPTOR'; 'pDEPTOR'; 
    'AMPK'; 'pAMPK'; 'SIRT1'; 'ULK1'; 'pULK1'};  

%Define the initial condition and time span for solving the system of ODEs
tt = 0:1000; %Timespan
y0 = [50;0;100;0;0;100;0;250;0;200;0;0;0;350;0;250;0;0;250;0]; %Initial condition

%Calculate the local sensitivty for each variable
Sensitivity = zeros(length(parameter_names), length(variable_names)); %Define matrix to store the sensitivities

%Solve the system of ODEs using the initial parameter set
[t, answer] = ode23s(@(t, x) dR2(t, x, param), tt, y0);
S_p0_full = answer(end, :); % S(p0)

%Loop through each variable and parameter and find the sensitivty by
%approximating the derivatives and store the results in the sensitivity matrix
for i = 1:length(variable_names)
    parfor j = 1:length(parameter_names)    
        
        %Set the vector param_modified with the modified parameter
        param_modified = param; %copy original parameter set
        p = param_modified(j); %Find the parameter of interest
        dp = 0.005 * p; %Set dp as 0.005*p
        param_modified(j) = p + dp; %Perturb parameter of interest
        
        %Get S(p0)
        S_p0_copy = S_p0_full; %Copy for parallel computing 
        S_p0 = S_p0_copy(i); %Find S(p0)
       
        %Solve the system of ODEs with the modified parameters,and find the
        %terminal point of the solution for the variable of interest
        [t, answer] = ode23s(@(t, x) dR2(t, x, param_modified), tt, y0);
        S_p0_dp = answer(end, i); % Get S(p0+dp)
        
        % Calculate the relative local sensitivity
        % Sensitivity = p/S * dS/dp
        % where dS/dp approximated by dS/dp = [S(p0+dp)-S(p0)]/dp
        Sensitivity(j, i) = p/S_p0*(S_p0_dp - S_p0)/dp; 
        
        %Print current i and j
        fprintf('i = %d, j = %d\n', i,j)
    end
end
 

%Construct the heat map
figure(1)
sensitivity_heatmap = heatmap(variable_names, parameter_names, Sensitivity, 'Colormap', jet);
sensitivity_heatmap.Title = 'Local Sensitivity Analysis';
% sensitivity_heatmap.GridVisible = 'off'; %Turn grid off
sensitivity_heatmap.FontSize = 8;

toc %stop timing
save localSAresults %Store simulation results

