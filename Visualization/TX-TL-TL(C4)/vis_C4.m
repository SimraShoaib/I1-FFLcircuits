clc
clear all

load Exp_C4.mat
With_J23110 = J23110;
With_J23116 = J23116;
With_T7PP7 = J23119;
With_noPP7 = noRBP;

tspan=0:0.1:420*60;  %%% seconds 660 => 480 needed
options = odeset('RelTol',1e-11,'AbsTol',1e-11);


% total plasmid concentrations
par.P_z = 1e-9;
par.P_y = 1e-9;
par.P_x = 1e-9; 

par.IPTG = 0.1*10^-3; %M

delay_time = 0*60;

load fit_C4.mat
p = p;

values_to_assign = [3*0.07, 0.07, 11*0.07,0.0]; %j=4, 0


for j = 1:4
    
    options = odeset('RelTol',1e-11,'AbsTol',1e-11);
    x0 = [0 0 0 0 0 0 0 0]; 
    p(2) = values_to_assign(j);
    [t,x] = ode23s(@(t,x)ODE_C4(t,x,p,par,delay_time),tspan,x0);
    x = x.*10^9;
    Simout_RBP = x(1:6000:end,8).*10^p(15); 
    tt = (0:length(With_T7PP7)-1).*10; %%% in minutes for visualization
    if j == 1
        Simout_1 = Simout_RBP;
        figure(1) 
        plot(tt,Simout_1,'-r',tt,With_J23110,'*r') 
        ylim([0 20000]); 
        xlabel("Time") 
        ylabel ("GFP") 
        set(gca,'FontSize',18)
        set(gca,'FontName','Times New Roman')
        hold on
    end
    if j == 2
        Simout_2 = Simout_RBP;
        figure(1) 
        plot(tt,Simout_2,'-b',tt,With_J23116,'*b') 
        ylim([0 20000]); 
        xlabel("Time") 
        ylabel ("GFP")
        set(gca,'FontSize',18)
        set(gca,'FontName','Times New Roman')
        hold on
    end
    if j == 3
        Simout_3 = Simout_RBP;
        figure(1) 
        plot(tt,Simout_3,'-k',tt,With_T7PP7,'*k') 
        xlim([0 450])
        ylim([0 20000])
        xlabel("Time") 
        ylabel ("GFP")
        set(gca,'FontSize',18)
        set(gca,'FontName','Times New Roman')
        hold on
    end
    if j == 4
        Simout_4 = Simout_RBP;
        figure(1)
        plot(tt,Simout_4,'-m',tt,With_noPP7,'*m')
        ylim([0 20000]); 
        xlabel("Time")
        ylabel ("GFP")
        set(gca,'FontSize',18)
        set(gca,'FontName','Times New Roman')
        hold off
        legend ("J23110-sim","J23110-exp","J23116-sim","J23116-exp", "J23119-sim","J23119-exp", "noPP7-sim","noPP7-exp")
    end


end


objective = sum((Simout_2 - With_J23116).^2)/max(With_J23116) + sum((Simout_2 - With_noPP7).^2)/max(With_noPP7);


disp(['SSE Objective: ' num2str(objective)])
disp(['Simu_End: ' num2str(Simout_RBP(end)),  'Exp: ' num2str(With_T7PP7(end))])


figure(2)
names = ["STAR","THS","Pyact","RBP","Z", "Zact", "Zrep", "GFP"];

Py = par.P_y*10^(9) - x(:,3); 

for i = 1:8
    subplot(3,3,i)
    plot(x(1:6000:end,i).*10^p(15),'LineWidth',2)
    title(names(i),'HorizontalAlignment','left')
    xlabel("Time (min)")
    set(gca,'FontSize',18)
    set(gca,'FontName','Times New Roman')
end

subplot(3,3,9)
plot(Py(1:6000:end).*10^p(15),'LineWidth',2)
title('Pyfree','HorizontalAlignment','left')
xlabel("Time (min)")
set(gca,'FontSize',18)
set(gca,'FontName','Times New Roman')
