%Global Sensitivity Analysis for the mTORC1/AMPK/ULK1 Pathway
clc; clear;            
tic %start timing

%Set number of random samples
N = 100000; % Number of parameter sets

%Parameter ranges for Km, Vmax, kc and kd
Km_range = [1 1000];
Vmax_range = [1e-3 10] * 3600; %Convert units from per sec -> per h
Kc_range = [1e-6 1] * 3600; %Convert units from per sec -> per h
Kd_range = [1e-6 10] * 3600; %Convert units from per sec -> per h

%Random sampling of parameter within the parameter range using
%Latin hypercube sampling
V_IR                 = (Vmax_range(end)-Vmax_range(1))*lhsdesign(N,1) + Vmax_range(1);
Km_IR                = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
V_pIR                = (Vmax_range(end)-Vmax_range(1))*lhsdesign(N,1) + Vmax_range(1);
Km_pIR               = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
K_IRS_by_pIR         = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
Km_IRS_by_pIR        = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
V_pIRS               = (Vmax_range(end)-Vmax_range(1))*lhsdesign(N,1) + Vmax_range(1);
Km_pIRS              = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
K_AKT_by_pIRS        = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
Km_AKT_by_pIRS       = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
K_AKT_by_pmTORC2     = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
Km_AKT_by_pmTORC2    = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
V_pAKT               = (Vmax_range(end)-Vmax_range(1))*lhsdesign(N,1) + Vmax_range(1);
Km_pAKT              = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
K_mTORC1_by_pAKT     = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
Km_mTORC1_by_pAKT    = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
K_pmTORC1            = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
K_pmTORC1_by_pAMPK   = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
Km_pmTORC1_by_pAMPK  = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
K_pmTORC1_by_pULK1   = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
Km_pmTORC1_by_pULK1  = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
K_mTORC2_by_pIRS     = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
Km_mTORC2_by_pIRS    = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
K_mTORC2_by_pAMPK    = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
Km_mTORC2_by_pAMPK   = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
V_pmTORC2            = (Vmax_range(end)-Vmax_range(1))*lhsdesign(N,1) + Vmax_range(1);
Km_pmTORC2           = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
K_DEPTOR_by_pmTORC1  = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
Km_DEPTOR_by_pmTORC1 = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
K_DEPTOR_by_pmTORC2  = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
Km_DEPTOR_by_pmTORC2 = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
V_pDEPTOR            = (Vmax_range(end)-Vmax_range(1))*lhsdesign(N,1) + Vmax_range(1);
Km_pDEPTOR           = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
K_mTORC1_DEPTOR_form = (Kd_range(end)-Kd_range(1))*lhsdesign(N,1) + Kd_range(1);
K_mTORC1_DEPTOR_diss = (Kd_range(end)-Kd_range(1))*lhsdesign(N,1) + Kd_range(1);
K_mTORC2_DEPTOR_form = (Kd_range(end)-Kd_range(1))*lhsdesign(N,1) + Kd_range(1);
K_mTORC2_DEPTOR_diss = (Kd_range(end)-Kd_range(1))*lhsdesign(N,1) + Kd_range(1);
K_IRS_to_iIRS        = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
Km_IRS_to_iIRS       = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
V_iIRS               = (Vmax_range(end)-Vmax_range(1))*lhsdesign(N,1) + Vmax_range(1);
Km_iIRS              = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
K_AMPK               = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1); 
K_AMPK_by_SIRT1      = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
Km_AMPK              = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
K_pAMPK              = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
K_pAMPK_by_pULK1     = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
K_pAMPK_by_pmTORC1   = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
Km_pAMPK             = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
K_SIRT1              = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
K_SIRT1_by_pAMPK     = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
Km_SIRT1             = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1);
K_SIRT1_diss         = (Kd_range(end)-Kd_range(1))*lhsdesign(N,1) + Kd_range(1);
K_ULK1               = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
K_ULK1_by_pAMPK      = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
Km_ULK1              = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1); 
K_pULK1              = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
K_pULK1_by_pmTORC1   = (Kc_range(end)-Kc_range(1))*lhsdesign(N,1) + Kc_range(1);
Km_pULK1             = (Km_range(end)-Km_range(1))*lhsdesign(N,1) + Km_range(1); 




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
    'IR'; 'PIR'; 'IRS'; 'pIRS'; 'iIRS'; 'AKT'; 
    'pAKT'; 'mTORC1'; 'pmTORC1'; 'mTORC2'; 'pmTORC2'; 
    'mTORC1-Deptor'; 'mTORC2-Deptor'; 'Deptor'; 'pDeptor'; 
    'AMPK'; 'pAMPK'; 'SIRT1'; 'ULK1'; 'pULK1'};  

%Construct matrix for sampled parameter sets
%Each row of the array param corresponds to a sampled set of parameters
param = [V_IR, Km_IR, V_pIR, Km_pIR, K_IRS_by_pIR, Km_IRS_by_pIR, V_pIRS,...
    Km_pIRS, K_AKT_by_pIRS, Km_AKT_by_pIRS, K_AKT_by_pmTORC2, Km_AKT_by_pmTORC2,...
    V_pAKT, Km_pAKT, K_mTORC1_by_pAKT, Km_mTORC1_by_pAKT, K_pmTORC1, K_pmTORC1_by_pAMPK,...
    Km_pmTORC1_by_pAMPK, K_pmTORC1_by_pULK1, Km_pmTORC1_by_pULK1, K_mTORC2_by_pIRS,...
    Km_mTORC2_by_pIRS, K_mTORC2_by_pAMPK, Km_mTORC2_by_pAMPK, V_pmTORC2, Km_pmTORC2, ...
    K_DEPTOR_by_pmTORC1, Km_DEPTOR_by_pmTORC1, K_DEPTOR_by_pmTORC2, Km_DEPTOR_by_pmTORC2, ...
    V_pDEPTOR, Km_pDEPTOR, K_mTORC1_DEPTOR_form, K_mTORC1_DEPTOR_diss, K_mTORC2_DEPTOR_form, ...
    K_mTORC2_DEPTOR_diss, K_IRS_to_iIRS, Km_IRS_to_iIRS, V_iIRS, Km_iIRS, K_AMPK, ...
    K_AMPK_by_SIRT1, Km_AMPK, K_pAMPK, K_pAMPK_by_pULK1, K_pAMPK_by_pmTORC1, Km_pAMPK, ...
    K_SIRT1, K_SIRT1_by_pAMPK, Km_SIRT1, K_SIRT1_diss, K_ULK1, K_ULK1_by_pAMPK, Km_ULK1, ...
    K_pULK1, K_pULK1_by_pmTORC1, Km_pULK1];

%Find number of parameters/variables
len_param = length(parameter_names); %Number of parameters
len_var = length(variable_names); %Number of variables

%Solve the system of ODEs for each sampled parameter set
tt = 0:3000; %Timespan
y0 = [50;0;100;0;0;100;0;250;0;200;0;0;0;350;0;250;0;0;250;0]; %Initial condition
results = zeros(N, len_var); %Pre-allocate matrix to store results

%Solve the system of ODEs and store the terminal concentration vector in time 
parfor i = 1:N   
    [t, answer] = ode23s(@(t, x) dR2(t, x, param(i,:)), tt, y0);
    results(i,:) = answer(end,:); %Capture terminal point of the simulation
    fprintf('i = %d\n',i) %Print current iteration number
end


%Calculate mean and variance
results_mean = mean(results)';
results_var = var(results)';
results_table = table( variable_names, results_mean, results_var ); %Show table of mean / variance


%Plot the Probability Density Distributions for each Output
fig = figure(1);
hold on
for i = 1:len_var
   subplot(5,4,i)
   histogram(results(:,i), 250, 'Normalization', 'pdf');
   title(variable_names{i})
end
hold off

%Set plot properties
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Probability');
xlabel(han,'Abundance');
sgtitle('Probability Density Functions of the Output', 'fontsize', 10)


%Partial Rank Correlation Coefficient
%Find the PRCC for each output 
rho_final = zeros(len_param, len_var); %Pre-allocate matrix 
pvalues_final = zeros(len_param, len_var); %Pre-allocate matrix

%Compute the PRCC for each variable with respect to the model parameters
parfor i = 1:length(variable_names)
    x = [param, results(:,i)]; %Construct matrix including the sampled parameter set 
                               %and the simulation results of each output
    [rho, pvalues] = partialcorr(x,'Type','Spearman'); %Calculate the PRCC
    rho_final(:,i) = rho(1:len_param, len_param+1); %Store the PRCC results 
    pvalues_final(:,i) = pvalues(1:len_param, len_param+1); %Store p-values
    fprintf('i = %d\n',i) %Print current iteration number
end


%Plot the PRCC heat map
figure(2)  
sensitivity_heatmap = heatmap(variable_names, parameter_names,rho_final,'colormap',jet);
% sensitivity_heatmap.GridVisible = 'off';
sensitivity_heatmap.FontSize = 8;
title('Global Sensitivity Analysis (PRCC)')


%Plot the p-value heat map
figure(3)
sensitivity_heatmap = heatmap(variable_names, parameter_names,pvalues_final,'colormap',jet);
title('p-value Heat Map')
sensitivity_heatmap.GridVisible = 'off';


%Plot PRCC heat map with only p < 0.05
rho_adjusted = rho_final; %copy PRCC matrix
for i = find(pvalues_final >= 0.05) %Set PRCC with p >= 0.05 to NaN
    rho_adjusted(i) = NaN;
end

figure(4)
sensitivity_heatmap = heatmap(variable_names, parameter_names,rho_adjusted,'colormap',jet);
title('PRCC Heat Map with p < 0.05')
sensitivity_heatmap.GridVisible = 'off';

%Plot 7D parallel plot 
figure(5)
idx = randi(N, 1, 150); %Randomly sample 150 simulation results out of N
parallelcoords(results(idx,[8,9,12,10,11,13, 15]) , 'labels', {'mTORC1', 'pmTORC1', 'mTORC1-DEPTOR', 'mTORC2', 'pmTORC2', 'mTORC2-DEPTOR', 'pDEPTOR'})
ylabel('Abundance')

toc %end timing
save globalSAresults %Store simulation results