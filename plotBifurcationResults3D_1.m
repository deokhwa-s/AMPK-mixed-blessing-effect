%Plot 3D bifurcation diagrams for pmTORC1 and pmTORC2 varying V_IR
%and AMPK from saved matrices from simulation
clc; clear;

%Import saved data for plotting
%Check if 'bifurcation_AMPK_V_IR_3D_DEPTOR_350.mat' is present. If not
%inform user to run bifurcationAnalysis1_3D.m
filename = 'bifurcation_AMPK_V_IR_3D_DEPTOR_350.mat';
if isfile(filename)
    data = open(filename); %Import saved simulation data
else
    warning('Run ''bifurcationAnalysis1_3D.m'' before running this file')
end

Abundance = data.Abundance; %Define AMPK abundance
max_val = data.output_max_val; %Get matrix containing the maximum values
min_val = data.output_min_val; %Get matrix containing the minimum values
sum_val = data.output_sum; %Get matrix containing the sum
product_val = data.output_product; %Get matrix containing the product
div_val =  data.output_div; %Get matrix containing the quotient
dim = size(max_val); %Find the dimensions of the matrix

%Define parameter range for V_IR
V_IR = data.param_range;
lim = [0 0.2]; %Define limit of V_IR to display
view_angle = [-30, 45]; %Define view angle of plot

%Define mesh for plotting scatter plots
[X, Y] = meshgrid(V_IR, Abundance);
X =  reshape(X, [dim(1)*dim(2), 1]);
Y =  reshape(Y, [dim(1)*dim(2), 1]);

%Reshape matrices
max_val_new = reshape(max_val, [dim(1)*dim(2), 2]);
min_val_new = reshape(min_val, [dim(1)*dim(2), 2]);
sum_val_new = reshape(sum_val, [dim(1)*dim(2), 2]);
product_val_new = reshape(product_val, [dim(1)*dim(2), 2]);
div_val_new = reshape(div_val, [dim(1)*dim(2), 2]);


%Plot 3D surface plots
%pmTORC1 
figure(1)
subplot(3,2,1)
surf(V_IR,Abundance,reshape(max_val(:, 1:length(V_IR), 1), [dim(1), dim(2)]), 'LineStyle','none')
hold on 
surf(V_IR,Abundance,reshape(min_val(:, 1:length(V_IR), 1), [dim(1), dim(2)]), 'LineStyle','none')
scatter3(X, Y, max_val_new(:,1), 1, 'k.')
scatter3(X, Y, min_val_new(:,1), 1, 'k.')
hold off
xlabel('V\_IR')
ylabel('AMPK')
zlabel('pmTORC1')
xlim(lim)
view(view_angle)

%pmTORC2
subplot(3,2,2)
surf(V_IR,Abundance,reshape(max_val(:, 1:length(V_IR), 2), [dim(1), dim(2)]), 'LineStyle','none')
hold on 
surf(V_IR,Abundance,reshape(min_val(:, 1:length(V_IR), 2), [dim(1), dim(2)]), 'LineStyle','none')
scatter3(X, Y, max_val_new(:,2), 1, 'k.')
scatter3(X, Y, min_val_new(:,2), 1, 'k.')
hold off
xlabel('V\_IR')
ylabel('AMPK')
zlabel('pmTORC2')
xlim(lim)
view(view_angle)

%pmTORC1 + pmTORC2
subplot(3,2,3)
surf(V_IR,Abundance,reshape(sum_val(:, 1:length(V_IR), 1), [dim(1), dim(2)]), 'LineStyle','none')
hold on 
surf(V_IR,Abundance,reshape(sum_val(:, 1:length(V_IR), 2), [dim(1), dim(2)]), 'LineStyle','none')
scatter3(X, Y, sum_val_new(:,1), 1, 'k.')
scatter3(X, Y, sum_val_new(:,2), 1, 'k.')
hold off
xlabel('V\_IR')
ylabel('AMPK')
zlabel('pmTORC1 + pmTORC2')
xlim(lim)
view([120, 5])

%pmTORC1 * pmTORC2
subplot(3,2,4)
surf(V_IR,Abundance,reshape(product_val(:, 1:length(V_IR) ,1), [dim(1), dim(2)]), 'LineStyle','none')
hold on 
surf(V_IR,Abundance,reshape(product_val(:, 1:length(V_IR), 2), [dim(1), dim(2)]), 'LineStyle','none')
scatter3(X, Y, product_val_new(:,1), 1, 'k.')
scatter3(X, Y, product_val_new(:,2), 1, 'k.')
hold off
xlabel('V\_IR')
ylabel('AMPK')
% zlabel('pmTORC1 X pmTORC2')
zlabel('pmTORC1 X pmTORC2')
xlim(lim)
view(view_angle)

%pmTORC1 / pmTORC2
subplot(3,2,5)
surf(V_IR,Abundance,reshape(div_val(:, 1:length(V_IR), 1), [dim(1), dim(2)]), 'LineStyle','none')
hold on 
surf(V_IR,Abundance,reshape(div_val(:, 1:length(V_IR), 2), [dim(1), dim(2)]), 'LineStyle','none')
scatter3(X, Y, div_val_new(:,1), 1, 'k.')
scatter3(X, Y, div_val_new(:,2), 1, 'k.')
hold off
xlabel('V\_IR')
ylabel('AMPK')
zlabel('pmTORC2 / pmTORC1')
xlim(lim)
view(view_angle)

%pmTORC2 / pmTORC1
subplot(3,2,6)
surf(V_IR,Abundance,1./reshape(div_val(:, 1:length(V_IR), 1), [dim(1), dim(2)]), 'LineStyle','none')
hold on 
surf(V_IR,Abundance,1./reshape(div_val(:, 1:length(V_IR), 2), [dim(1), dim(2)]), 'LineStyle','none')
scatter3(X, Y, 1./div_val_new(:,1), 1, 'k.')
scatter3(X, Y, 1./div_val_new(:,2), 1, 'k.')
hold off
xlabel('V\_IR')
ylabel('AMPK')
zlabel('pmTORC2 / pmTORC1')
xlim(lim)
view(view_angle)

sgtitle(append('Bifurcation Plots of pmTORC1 and pmTORC2 at DEPTOR = ', num2str(data.y0(14))), 'fontsize', 10)