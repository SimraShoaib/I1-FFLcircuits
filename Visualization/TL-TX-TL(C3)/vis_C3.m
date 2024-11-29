clc
clear all

load Exp_C3.mat


With_RBP_AAV = AAV;
With_RBP_ASV = ASV;
With_RBP_LVA = LVA;
With_RBP_0 = noRBP;

tspan=0:0.1:420*60;  %%% seconds
options = odeset('RelTol',1e-11,'AbsTol',1e-11);


% total plasmid concentrations
par.P_z = 1e-9;
par.P_y = 1e-9;
par.P_x = 1e-9; 

par.IPTG = 0.1*10^-3; %M


load fit_C3.mat
p=p;
fit=p(9); 

delay_time = 0;  

values_to_assign = [2*(fit),fit,12*(fit),0];


 for j = 1:4
   
    new_p = p; % new parameter for RBP = 0 
    p(9) = values_to_assign(j);

    if j == 1
        x0 = [0 0 0 0 0 0 0 0 0];
        [t,x] = ode23s(@(t,x)ODE_C3(t,x,p,par,delay_time),tspan,x0);
        x = x.*10^9; 
        Simout_RBP = x(1:6000:end,8).*10^p(16);
        tt = ((0:length(With_RBP_0)-1)*10)+70; 
        Simout_1 = Simout_RBP;
        figure(1)
        plot(tt,Simout_1,'-r',tt,With_RBP_AAV,'*r') 
        ylim([0 2800]);
        xlabel("Time (min)") 
        ylabel ("GFP")
        hold on
    end
    if j == 2
        x0 = [0 0 0 0 0 0 0 0 0];
        [t,x] = ode23s(@(t,x)ODE_C3(t,x,p,par,delay_time),tspan,x0);
        x = x.*10^9; 
        Simout_RBP = x(1:6000:end,8).*10^p(16);
        tt =((0:length(With_RBP_0)-1)*10)+70; 
        Simout_2 = Simout_RBP;
        figure(1)
        plot(tt,Simout_2,'-b',tt,With_RBP_ASV,'*b') 
        ylim([0 2800]);
        xlabel("Time (min)") 
        ylabel ("GFP")
        hold on
    end

    if j == 3
        x0 = [0 0 0 0 0 0 0 0 0];
        [t,x] = ode23s(@(t,x)ODE_C3(t,x,p,par,delay_time),tspan,x0);
        x = x.*10^9; 
        Simout_RBP = x(1:6000:end,8).*10^p(16);
        tt = ((0:length(With_RBP_0)-1)*10)+70; 
        Simout_3 = Simout_RBP;
        figure(1)
        plot(tt,Simout_3,'-k',tt,With_RBP_LVA,'*k') 
        ylim([0 2800]);
        xlabel("Time (min)") 
        ylabel ("GFP")
        hold on
    end

    if j == 4 
        x0 = [0 0 0 0 0 0 0 0 0];
        new_p(3) = 0.17; 
        [t,x] = ode23s(@(t,x)ODE_C3(t,x,new_p,par,delay_time),tspan,x0); %new_p
        x = x.*10^9; 
        Simout_RBP = x(1:6000:end,8).*10^new_p(16);
        tt =((0:length(With_RBP_0)-1)*10)+70; 
        Simout_4 = Simout_RBP;
        figure(1)
        plot(tt,Simout_4,'-m',tt,With_RBP_0,'*m') 
        ylim([0 2800]);
        xlabel("Time (min)") 
        ylabel ("GFP") 
        hold off
        legend ("AAV sim","AAV exp","ASV sim","ASV exp", "LVA sim","LVA exp", "noRBP sim","noRBP exp")
    end
   
 end


objective = sum((Simout_RBP - With_RBP_ASV).^2)/max(With_RBP_ASV);


disp(['SSE Objective: ' num2str(objective)])
disp(['Simu_End: ' num2str(Simout_RBP(end)),  'Exp: ' num2str(With_RBP_ASV(end))])


figure(2)
names = ["STAR","THS","Y","Yact","RBP", "Pzactive", "Z", "GFP","Zrep"];

Pz = par.P_z*10^(9) - x(:,6); 
Z_free = x(:,7) - x(:,9);

for i = 1:9
    subplot(3,4,i)
    plot(x(1:6000:end,i).*10^p(16),'LineWidth',2)
    title(names(i),'HorizontalAlignment','left')
    xlabel("Time (min)")
    set(gca,'FontSize',18)
    set(gca,'FontName','Times New Roman')
end

subplot(3,4,10)
plot(Pz(1:6000:end).*10^p(16),'LineWidth',2)
title('Pzfree','HorizontalAlignment','left')
xlabel("Time (min)")
set(gca,'FontSize',18)
set(gca,'FontName','Times New Roman')


subplot(3,4,11)
plot(Z_free(1:6000:end).*10^p(16),'LineWidth',2)
title('Z_free','HorizontalAlignment','left')
xlabel("Time (min)")
set(gca,'FontSize',18)
set(gca,'FontName','Times New Roman')


figure(3)
p(9) = 12*(fit); 
% New figure for LVA simulations with different delay times
delay = [0, 20, 40, 80, 120, 240].*60; % delay times in seconds

hold on

for d = 1:length(delay)
    delay_time = delay(d); % current delay time
    x0 = [0 0 0 0 0 0 0 0 0]; % initial conditions
    
    % Run the ODE solver with the current delay time
    [t,x] = ode23s(@(t,x)ODE_C3(t,x,p,par,delay_time),tspan,x0);
    x = x .* 10^9; % Convert to nM for visualization
    Simout_RBP = x(1:6000:end, 8) .* 10^p(16); % Get the relevant output
    
    % Visualization
    tt = (0:length(With_RBP_LVA)-1) .* 10; % Time vector for experimental data
    
    % Plot the simulation data
    plot(tt+70, Simout_RBP, '-','DisplayName', sprintf('delay %d min', delay_time./60),'LineWidth',2);
    
    % For the experimental data, plot J23110
    if d == 6
        plot(tt+70, With_RBP_LVA, '*r', 'DisplayName', 'LVA-exp','LineWidth',1.5); 
    end
end


% Finalize the new figure
xlim([0 450])
ylim([0 15000]) 
xlabel("Time (minutes)") 
ylabel("GFP (nM)")
set(gca, 'FontSize', 18)
set(gca, 'FontName', 'Times New Roman')
legend show 
hold off
