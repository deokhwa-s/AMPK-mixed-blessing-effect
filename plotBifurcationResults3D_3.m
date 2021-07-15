%Plot 3D bifurcation diagrams for pmTORC1 and pmTORC2 varying IC of DEPTOR
%and AMPK from saved matrices from simulation
clc; clear;
    
%Import saved data for plotting
%Check if 'bifurcation_AMPK_DEPTOR_V_IR_3e-03.mat' is present. If not
%inform user to run bifurcationAnalysis2_3D.m
filename = 'bifurcation_AMPK_DEPTOR_V_IR_3e-03.mat';
if isfile(filename)
    data = open(filename); %Import saved simulation data
else
    warning('Run ''bifurcationAnalysis3_3D.m'' before running this file')
end

Abundance = data.Abundance; %Define AMPK abundance
max_val = data.output_max_val; %Get matrix containing the maximum values
min_val = data.output_min_val; %Get matrix containing the minimum values
sum_val = data.output_sum; %Get matrix containing the sum
product_val = data.output_product; %Get matrix containing the product
div_val =  data.output_div / (data.y0(8) / data.y0(10)); %Get matrix containing the quotient and normalize
dim = size(max_val); %Find the dimensions of the matrix

%Define mesh for plotting scatter plots
[X, Y] = meshgrid(Abundance, Abundance);
X =  reshape(X, [dim(1)*dim(2), 1]);
Y =  reshape(Y, [dim(1)*dim(2), 1]);
div_val_new = reshape(div_val, [dim(1)*dim(2), 2]);

%pmTORC1 / pmTORC2
figure(1)
surf(Abundance,Abundance,reshape(div_val(:, :, 1), [dim(1), dim(2)]), 'LineStyle','none')
hold on 
surf(Abundance,Abundance,reshape(div_val(:, :, 2), [dim(1), dim(2)]), 'LineStyle','none')
scatter3(X, Y, div_val_new(:,1), 0.1, 'k.')
scatter3(X, Y, div_val_new(:,2), 0.1, 'k.')
hold off
xlabel('DEPTOR')
ylabel('AMPK')
zlabel('pmTORC1 / pmTORC2')

title(append('Bifurcation Plot of pmTORC1 / pmTORC2 at V\_IR = ', num2str(data.param(1)/3600)), 'fontsize', 10)