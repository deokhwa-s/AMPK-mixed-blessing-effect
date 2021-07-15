%Check effect of SIRT1 on pmTORC1, pmTORC2, AMPK, ULK1
%Solving the system of ODEs at two different SIRT1 level
clc; clear;

%Set parameters for simulation
param = importdata('modelParameters.txt');

%Set parameters for simulation
tt = 0:200; %timespan
y0 = [50;0;100;0;0;100;0;250;0;200;0;0;0;350;0;250;0;0;250;0]; %Initial condition
idx = [17, 9, 11]; %Indicies for pmTORC1, pmTORC2, pAMPK, pULK1 in the ouput 

%Define V_IR values to evaulate
V_IR_values = [0.1 * 3600,  0.005 * 3600, 9.673148];
title_name = {'pAMPK', 'pmTORC1', 'pmTORC2'};

%Solve the system of ODEs at the defined V_IR values and plot the results
figure(1)
set(gcf, 'Position',  [100, 100, 800, 700])

for i = 1:3 
    %Set V_IR
    param(1) = V_IR_values(i); 
    
    %Solve the ODE at SIRT = 375
    [t, answer] = ode23s(@(t, x) dR2(t, x, param, 375), tt, y0);
    modelOutput1 = answer(:,idx); %Obtain model output

    %Modify to low SIRT = 37.5 and solve the system of ODEs
    [t, answer] = ode23s(@(t, x) dR2(t, x, param, 37.5), tt, y0);
    modelOutput2 = answer(:,idx);%Obtain model output
    
    %Plot results for pAMPK, pmTORC1 and pmTORC2
    for j = 1:3
        subplot(3,3,3*(j-1)+i);
        hold on
        
        %Plot results
        plot(t, modelOutput1(:,j)./max(modelOutput1(:,j)),'b', 'displayname', 'SIRT1 = 375')
        plot(t, modelOutput2(:,j)./max(modelOutput1(:,j)),'r', 'displayname', 'SIRT1 = 37.5')

        
        hold off
        
        %Add x, y labels
        xlabel('Time (h)')
        ylabel(title_name(j))
        
        %Add title for the plots in the first row
        if j == 1
            text(0.5, 1.1, append('V\_IR = ',sprintf('%0.2e',param(1)/3600)), 'Units', 'normalized',...
                'VerticalAlignment','Bottom', 'HorizontalAlignment','center','FontWeight','bold', 'FontSize', 11);
        end
        
        %Adjust y limits
        if j == 1
            ylim([0 1.7])
        else
            ylim([0 2])
        end
        
        legend
    end
end
