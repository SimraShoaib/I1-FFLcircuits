clc
clear all

load aTc50_C1.mat

With_TetR_1_10 = J23110;
With_TetR_1_16 = J23116;
With_TetR_1_19 = J23119;
With_TetR_1_0 = noTetR;

tspan=0:0.1:420*60;  %%% seconds % 420 mins (25200 sec) with 0.1 sec apart
options = odeset('RelTol',1e-11,'AbsTol',1e-11);


% total plasmid concentrations
par.P_z = 1e-9;
par.P_y = 1e-9;
par.P_x = 1e-9; 

atc_conv = 0.46822;
par.IPTG = 0.1*10^-3; %M


load fit_C1.mat
p=p;
values_to_assign = [3*0.0043,0.0043,11*0.0043, 0];

% p(2) ratio for promoter strength
% J23110 : J23116 : J23119 = 3 : 1 : >10

for j = 1:4
    
    options = odeset('RelTol',1e-11,'AbsTol',1e-11);
    par.aTc = (50/atc_conv)*10^-9; %M
    p(2) = values_to_assign(j);
    x0 = [0 0 0 par.aTc 0 0 0 0 0 0];
    [t,x] = ode23s(@(t,x)ODE_C1(t,x,p,par),tspan,x0);
    x = x.*10^9; %%% convert to nM for visualization
    Simout_TetR_j = x(1:6000:end,8).*10^p(17); % every 10 mins
    tt = ((0:length(With_TetR_1_0)-1).*10)+70; %%% for visualization

    if j == 1
        Simout_1 = Simout_TetR_j;
        figure(1) 
        plot(tt,Simout_1,'-r',tt,With_TetR_1_10,'*r') 
        ylim([0 2000]); 
        xlabel("Time (min)") 
        ylabel ("GFP") 
        hold on
    end
    if j == 2
        Simout_2 = Simout_TetR_j;
        figure(1) 
        plot(tt,Simout_2,'-b',tt,With_TetR_1_16,'*b') 
        ylim([0 2000]); 
        xlabel("Time (min)") 
        ylabel ("GFP") 
        hold on
    end
    if j == 3
        Simout_3 = Simout_TetR_j;
        figure(1) 
        plot(tt,Simout_3,'-k',tt,With_TetR_1_19,'*k') 
        ylim([0 2000]); 
        xlabel("Time (min)") 
        ylabel ("GFP") 
        hold on
    end
    if j == 4
        Simout_4 = Simout_TetR_j;
        figure(1)
        plot(tt,Simout_4,'-m',tt,With_TetR_1_0,'*m')
        ylim([0 2000]); 
        xlabel("Time (min)")
        ylabel ("GFP")
        hold off 

        legend("J23110-sim","J23110-exp","J23116-sim","J23116-exp","J23119-sim","J23119-exp","noTetR-sim","noTetR-exp")
        
    end
end


objective =  sum((Simout_4 - With_TetR_1_0).^2)/max(With_TetR_1_0) ;




disp(['SSE Objective: ' num2str(objective)])
disp(['Simu1_End: ' num2str(Simout_1(end)),  'Exp1: ' num2str(With_TetR_1_10(end))])
disp(['Simu2_End: ' num2str(Simout_2(end)),  'Exp2: ' num2str(With_TetR_1_16(end))])
disp(['Simu3_End: ' num2str(Simout_3(end)),  'Exp3: ' num2str(With_TetR_1_19(end))])
disp(['Simu4_End: ' num2str(Simout_4(end)),  'Exp4: ' num2str(With_TetR_1_0(end))])


% for individual species

figure(2)
names = ["STAR","THS","TetR","aTc","Pzrep","Z","Zact","GFP","Pyactive", "aTc:TetR"];

Pz = par.P_z*10^(9) - x(:,5); 
Py = par.P_y*10^(9) - x(:,9);

for i = 1:10
    subplot(3,4,i)
    plot(x(1:6000:end,i).*10^p(17),'LineWidth',2)
    title(names(i),'HorizontalAlignment','left')
    xlabel("Time (min)")
    set(gca,'FontSize',18)
    set(gca,'FontName','Times New Roman')
end

subplot(3,4,11)
plot(Pz(1:6000:end).*10^p(17),'LineWidth',2)
title('Pzfree','HorizontalAlignment','left')
xlabel("Time (min)")
set(gca,'FontSize',18)
set(gca,'FontName','Times New Roman')

subplot(3,4,12)
plot(Py(1:6000:end).*10^p(17),'LineWidth',2)
title('Pyfree','HorizontalAlignment','left')
xlabel("Time (min)")
set(gca,'FontSize',18)
set(gca,'FontName','Times New Roman')
