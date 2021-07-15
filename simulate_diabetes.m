%Check effect of diabetes on pmTORC1 level
%by decreasing IR level
clc; clear;

%Set parameters for simulation
param = importdata('modelParameters.txt'); %Import model parameters
param(1) = 0.1 * 3600; % Set V_IR = 0.1 (converted from per s -> per h)
tt = 0:1000; %Set timespan
y0 = [50;0;100;0;0;100;0;250;0;200;0;0;0;350;0;250;0;0;250;0]; %Initial condition

%Solve the system of ODEs at IR = 50
[t1, answer] = ode23s(@(t, x) dR2(t, x, param), tt, y0);
pmTORC1_normal = answer(:,9); %Obtain pmTORC1 output

%Solve the ODE again at decreased IR
y0(1) = 25;

%Solve the system of ODEs at IR = 25
[t2, answer] = ode23s(@(t, x) dR2(t, x, param), tt, y0);
pmTORC1_diabetic = answer(:,9);


%Normalized Experimental data from Ost et al.
%first row = ( pmTORC1 (normal), pmTORC1 (diabetic) )
%second row = corresponding length of error bar
experimental_data = [1 , 3.06/4.08; 0.519/4.08, 0.227/4.08];                                                             
                                                             
%Plot resulting pmTORC1 against experimental data
figure(1)                                                                                                                         
h1 = bar([ 1, 0; pmTORC1_diabetic(end)/pmTORC1_normal(end), 0 ], 'barwidth', 1);
set(gca, 'XTickLabel', {'IR = 50','IR = 25'});
a = get(gca, 'Position');
set(h1, {'DisplayName'}, {'Model prediction', 'Experimental data'}')
h1(1).FaceColor = [0.0745    0.6235    1.0000];

%Get x-coordiantes of the bars to plot the scatter plot
hold on
x = nan(2, 2);
for i = 1:2
    x(i,:) = h1(i).XEndPoints;
end

%Plot scatter plot with error bar
errorbar(x([2,4]), experimental_data (1,:), experimental_data (2,:), 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'Color', 'r',...
            'LineWidth', 0.7, 'displayname', 'Experimental data')

hold off
h=get(gca,'Children'); % grab all the axes handles at once
legend(h([3 1]),{'Model prediction', 'Experimental data'}) 
ylabel('pmTORC1');
text(sum(x(1:2))/2, -0.1, 'non-diabetic', 'HorizontalAlignment','center')
text(sum(x(3:4))/2, -0.1, 'diabetic', 'HorizontalAlignment','center')


