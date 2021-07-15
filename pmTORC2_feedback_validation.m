%Validation of activation strength of AMPK on mTORC2
%by simulating the increase in AKT activation with increasing AMPK
%simulate case of AICAR injection; set V_IR = 1.15e-6
clc; clear; 

%Set model parameters 
param = importdata('modelParameters.txt');
param(1) = 1.15e-6 * 3600; %Set V_IR = 1.15e-6 and convert units from per s -> per h

%Set initial parameters for simulation
tt = 0:1000; %Timespan
y0 = [50;0;100;0;0;100;0;250;0;200;0;0;0;350;0;250;0;0;250;0]; %Initial condition
[t, answer] = ode23s(@(t, x) dR2(t, x, param), tt, y0); %Solve the system of ODEs

%Store simulation results
pAKT_initial = answer(:,7)./answer(:,6);
pAMPK_initial = answer(:,17)./answer(:,16);
pmTORC2_1 = answer(:,11)./answer(:,10);

%Simulate AICAR injection by increasing AMPK activation
param(42:43) = param(42:43) * 2.475; %Increase activation rate of AMPK (factor has been determined by trial-and-error)
[t, answer] = ode23s(@(t, x) dR2(t, x, param), tt, y0); %Solve the system of ODEs

%Store simulation results
pAKT_final = answer(:,7)./answer(:,6);
pAMPK_final = answer(:,17)./answer(:,16);
pmTORC2_2 = answer(:,11)./answer(:,10);


%model predition is normaized by the initial pAMPK and pAKT concentratins 
model_prediction = [1, pAMPK_final(end)/pAMPK_initial(end); 1, pAKT_final(end)/pAKT_initial(end)];
experimental_data = [1 , 2.82; 1, 1.52]; %Experimental data 
experimental_data_errorbars = [0.14, 0.65; 0.09, 0.214]; %Error bars for the experimental data 

%Compare simulation results to experimental data
figure(1) 
hold on

%Create bar plot for model prediction
h1 = bar(model_prediction, 'barwidth', 1);
h1(1).FaceColor = [0.0745    0.6235    1.0000];
h1(2).FaceColor = [1.0000    1.0000    0.0667];
set(gca, 'XTickLabel', {'pAMPK/AMPK','pAKT/AKT'});
set(gca, 'XTick', [1, 2])
set(h1, {'DisplayName'}, {'Control','AICAR'}')
ylim([0, 4])
ylabel('Abundance');


%Plot scatter plot with error bars for experimental results

%Obtain x positions of the bars 
x_position = zeros(2, 2);
for i = 1:2
    x_position(i,:) = h1(i).XEndPoints;
end

%Plot scatter plot 
errorbar(x_position, experimental_data' , experimental_data_errorbars', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'Color', 'r',...
            'LineWidth', 0.7, 'displayname', 'Experimental data')
hold off

%Set appropriate legend
h=get(gca,'Children'); % grab all the axes handles at once
legend(h([4 3 2]),{'Control', 'AICAR', 'Experimental data'}) 