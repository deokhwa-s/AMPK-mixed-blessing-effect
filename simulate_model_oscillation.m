%Simulate model osciilation
clc; clear;

%Import oscillation data from from Holczer et al.
%Each row corresponds to a datapoint in time
%    Data repeated twice from (copied 2 - 24 h to 26 - 48 h)
%    Each row corresponds to a datapoint in time
%    First column  : Time (h) 
%    Second column : Concentration (nM)
%    Third column  : Length of error bar (nM)

pmTORC1_data = importdata('pmTORC1_data.txt'); %Data for pmTORC1
pAMPK_data = importdata('pAMPK_data.txt'); %Data for pAMPK
pULK1_data = importdata('pULK1_data.txt');  %Data for pULK1

%Import parameter set as a column vector
param = importdata('modelParameters.txt');

% Set initial condition
y0 = [45.2678886688302,4.73211133116980,17.9912537197180,10.4134176890683,71.5953285901381,28.9130204869018,71.0869795130983,87.9364799719260,9.14262282859870,11.3168282120687,5.02897769163197,152.920897222202,183.654194080837,13.1620960677460,0.262812643184197,240.279483201148,9.72051679885231,37.0602139519015,199.632839162159,50.3671608378409];

%Solve the ODE using ode23s
tt = 0:500; %Set timespan
[t, answer] = ode23s(@(t, x) dR2(t, x, param), tt, y0); %Solve the ODE

%Find the resulting outputs for pmTORC1, pAMPK, and pULK1
pmTORC1_modelOutput = answer(:, 9); 
pAMPK_modelOutput = answer(:,17);
pULK1_modelOutput = answer(:,20);


%Plot model prediction and experimental data
figure(1)

%1. pmTORC1
subplot(2,2,1)
plot(t, pmTORC1_modelOutput, 'Displayname', 'Model Prediction')
hold on
errorbar(pmTORC1_data(:,1), pmTORC1_data(:,2),pmTORC1_data(:,3), 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'Color', 'r',...
            'LineWidth', 0.7, 'displayname', 'Experimental data')
hold off
xlabel('Time(h)')
ylabel('pmTORC1')
xlim([0, 50])
ylim([0, 100])
legend

%2. pAMPK
subplot(2,2,2) 
plot(t, pAMPK_modelOutput, 'Displayname', 'Model Prediction')
hold on
errorbar(pAMPK_data(:,1), pAMPK_data(:,2), pAMPK_data(:,3), 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'Color', 'r',...
            'LineWidth', 0.7, 'displayname', 'Experimental data')
hold off
xlabel('Time(h)')
ylabel('pAMPK')
xlim([0, 50])
ylim([0, 130])
legend

%3. pULK1
subplot(2,2,3) 
plot(t, pULK1_modelOutput, 'Displayname', 'Model Prediction')
hold on
errorbar(pULK1_data(:,1), pULK1_data(:,2), pULK1_data(:,3), 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'Color', 'r',...
            'LineWidth', 0.7, 'displayname', 'Experimental data')
hold off
xlabel('Time(h)')
ylabel('pULK1')
xlim([0, 50])
ylim([0, 80])
legend

%4. Combined plot
subplot(2,2,4) 
hold on
plot(t, pAMPK_modelOutput, 'displayname', 'pAMPK')
plot(t, pULK1_modelOutput, 'displayname', 'pULK1')
plot(t, pmTORC1_modelOutput, 'displayname', 'pmTORC1')
hold off

xlim([0, 50])
ylim([0, 120])
xlabel('Time(h)')
ylabel('Abundance')
legend