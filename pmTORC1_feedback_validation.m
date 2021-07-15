%Validation of the  inhibition strength of AMPK on mTORC1
%by simulating the increase in pmTORC1 activation with decreasing AMPK
%set V_IR_act = 0.1 to avoid oscillatory regime
clc; clear;

%Set parameters for simulation
param = importdata('modelParameters.txt'); %Import model parameters 
param(1) = 0.1 * 3600; % Set V_IR = 0.1 (converted from per s -> per h)
tt = 0:1000; %Set timespan
y0 = [50;0;100;0;0;100;0;250;0;200;0;0;0;350;0;250;0;0;250;0]; %Initial condition
factors = linspace(0.79,1,10); %Factors to simulate fractional AMPK activation
pAMPK = zeros(length(factors),1); %Vector to store resulting pAMPK/AMPK 
pmTORC1 = zeros(length(factors),1); %Vector to store resulting pmTORC1/mTORC1
    
%Solve the ODE varying the activation rate of AMPK using the multiplication
%factors
for i = 1:length(factors)
    %Reduce K_AMPK and K_AMPK_by_SIRT1  
    param(42:43) = param(42:43) * factors(i); %Change activation strength of AMPK
    [t, answer] = ode23s(@(t, x) dR2(t, x, param), tt, y0); %Solve the ODE
    
    %Store resulting activation levels of AMPK and mTORC1
    pmTORC1(i) = answer(end,9)./answer(end,8); % pmTORC1/mTORC1
    pAMPK(i) = answer(end,17)./answer(end,16); % pAMPK/AMPk
end


%Normalized Experimental data from Pal et el.
%each row = ( [pAMPK]/[AMPK], [pmTORC1]/[mTORC1] )
experimental_data = [
    1	1
    0.549450549	1.207909915
    0.412087912	1.298819006
    0.285714286	1.6918429
    0.098901099	1.810285636     
];

%Plot model prediction and experimental data
figure(1)
hold on

%Plot smooth curve connecting the data points from the model prediction
%using MATLAB function spline
X = 0.09:.01:1;
Y = spline(pAMPK(1:end-1)./pAMPK(1),pmTORC1(1:end-1)./pmTORC1(1),X); %Plot smooth curve normalizing by the first data point
plot(X, Y, 'displayname', 'Model prediction')

%Plot experimental data
scatter(experimental_data(:,1), experimental_data(:,2), 'or', 'filled', 'displayname', 'Experimental data')

hold off
xlabel('pAMPK/AMPK')
ylabel('pmTORC1/mTORC1')
legend
