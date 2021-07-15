%Plot 3D bifurcation diagrams for pmTORC1 and pmTORC2 varying V_pmTORC2
%and AMPK from saved matrices from simulation
clc; clear;

%Set initial parameter set as a column vector
param = importdata('modelParameters.txt');
y0 = [50;0;100;0;0;100;0;250;0;200;0;0;0;350;0;250;0;0;250;0]; %Initial condition

%Solve the ODE
tt = 0:1000; %Set timespan
[t, answer] = ode23s(@(t, x) dR2(t, x, param), tt, y0); %Solve the ODE

%Find baseline by computing the average of the max and min of the last 100
%points in time to normalize the plots
pmTORC1_baseline = (max(answer(end-100:end,9))+min(answer(end-100:end,9)))/2;
pmTORC2_baseline = (max(answer(end-100:end,11))+min(answer(end-100:end,11)))/2;


%Import saved data for plotting
%Check if 'bifurcation_V_pmTORC2.mat' is present. If not
%inform user to run bifurcationAnalysis2_3D.m
filename = 'bifurcation_V_pmTORC2.mat';
if isfile(filename)
    data = open(filename); %Import saved simulation data
else
    warning('Run ''bifurcationAnalysis2_3D.m'' before running this file')
end


%Import saved data for plotting
data = open('bifurcation_V_pmTORC2.mat'); %Import saved simulation data
Abundance = data.Abundance; %Define AMPK abundance
max_val = data.output_max_val; %Get matrix containing the maximum values
min_val = data.output_min_val; %Get matrix containing the minimum values
sum_val = data.output_sum; %Get matrix containing the sum
product_val = data.output_product; %Get matrix containing the product
div_val =  data.output_div; %Get matrix containing the quotient
div_val_norm =  div_val / pmTORC1_baseline * pmTORC2_baseline;
dim = size(max_val); %Find the dimensions of the matrix

%Define parameter range for V_pmTORC2
V_pmTORC2 = data.param_range;
lim = [0 0.02]; %Define limit of V_pmTORC2 to display


%Define mesh for plotting scatter plots
[X, Y] = meshgrid(V_pmTORC2, Abundance);
X =  reshape(X, [dim(1)*dim(2), 1]);
Y =  reshape(Y, [dim(1)*dim(2), 1]);

%Reshape matrices
max_val_new = reshape(max_val, [dim(1)*dim(2), 2]);
min_val_new = reshape(min_val, [dim(1)*dim(2), 2]);
sum_val_new = reshape(sum_val, [dim(1)*dim(2), 2]);
product_val_new = reshape(product_val, [dim(1)*dim(2), 2]);
div_val_new = reshape(div_val, [dim(1)*dim(2), 2]);
div_val_norm_new = reshape(div_val_norm, [dim(1)*dim(2), 2]);

%Create figure and define position
f = figure(1);
set(gcf, 'Position',  [100, 100, 1000, 400])

%pmTORC1 + pmTORC2
subplot(1,2,1)
surf(V_pmTORC2,Abundance,reshape(sum_val(:, 1:length(V_pmTORC2), 1), [dim(1), dim(2)]), 'LineStyle','none')
hold on 
surf(V_pmTORC2,Abundance,reshape(sum_val(:, 1:length(V_pmTORC2), 2), [dim(1), dim(2)]), 'LineStyle','none')
scatter3(X, Y, sum_val_new(:,1), 1, 'k.')
scatter3(X, Y, sum_val_new(:,2), 1, 'k.')
hold off
xlabel('V\_pmTORC2')
ylabel('AMPK')
zlabel('pmTORC1 + pmTORC2')
xlim(lim)
view([35, 7])

%pmTORC1 / pmTORC2 normalized
subplot(1,2,2)
surf(V_pmTORC2,Abundance,reshape(div_val_norm(:, 1:length(V_pmTORC2), 1), [dim(1), dim(2)]), 'LineStyle','none')
hold on 
surf(V_pmTORC2,Abundance,reshape(div_val_norm(:, 1:length(V_pmTORC2), 2), [dim(1), dim(2)]), 'LineStyle','none')
scatter3(X, Y, div_val_norm_new(:,1), 1, 'k.')
scatter3(X, Y, div_val_norm_new(:,2), 1, 'k.')
hold off
xlabel('V\_pmTORC2')
ylabel('AMPK')
zlabel('pmTORC1 / pmTORC2')
xlim(lim)
view([40, 25])

sgtitle('Bifurcation Plots of pmTORC1 and pmTORC2 varying V\_pmTORC2 and AMPK', 'fontsize', 10)