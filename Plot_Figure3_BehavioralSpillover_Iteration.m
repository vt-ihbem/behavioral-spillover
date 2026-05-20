%% Plot trajectories at various spillover scenarios
% Need to change p.R02 to get different values of R_0,B

clearvars; close all;
clc
saveplot = 0;

%% Standard parameters

%%%%%% Standard parameters
p.tauF = 30;
p.k = 100;
p.tauI1  = 7;
p.tauI2  = p.tauI1;
p.w1 = 1;
p.mu = 0;% natural birth/death
p.tauR = 1/100; % waning immunity

%% Two interdependent viruses with spillover

%%%%%% Parameters
p.type = 'imperfect';
p.R01 = 3;
p.R02 = 1.3;
p.w2 = 1;

p.beta1  = p.R01/p.tauI1;
p.beta2  = p.R02/p.tauI2;

svec = [0 .1 0.5 .9 1];

%%%%%% Initial Conditions
IC = InitCond();

%%%%%% Solvers and solutions

colormap('copper')

tmax = 365;

reout1 = NaN(tmax,length(svec));
reout2 = NaN(tmax,length(svec));

for s_it = 1:length(svec)
    p.s = svec(s_it);
    options = odeset('NonNegative',1:8);
    [t,y] = ode15s(@sirb, 0:tmax, IC, options, p);
    i1 = y(:,2); i2 = y(:,6);


    % Dynamic plots for infectious population
    f3 = figure(3);
    plot(t,i1,'--','linewidth',2)

    xlim([0 tmax])
    hold on
    plot(t,i2,'linewidth',2)
    hold off
    set(gca,'fontsize',24,'ytick',0:.01:0.05,'fontname','Times')
    xlabel('Time, days','interpreter','latex')
    ylabel('Infectious population, $I$','interpreter','latex')
    l = legend(sprintf('Disease A, $\\mathcal{R}_{0,A}$ = %0.1f',p.R01) ,sprintf('Disease B, $\\mathcal{R}_{0,B}$ = %0.1f',p.R02),'interpreter','latex');
    set(l,'box','off')
    if p.s ==0
        title('No spillover ($s=0$)','interpreter','latex')
    elseif p.s==1
        title('Perfect spillover ($s=1$)','interpreter','latex')
    else
        title(sprintf('Imperfect spillover ($s=%0.1f$)',p.s),'interpreter','latex')
    end
    ylim([0 0.05])
    
    if saveplot == 1
    figname = strcat(sprintf('Figures/Traj_R0B_%d_s%d',round(p.R02*10),round(p.s*100)),'.png');
    saveas(f3,figname,'png')
    figname1 = strcat(sprintf('Figures/tiff/Traj_R0B_%d_s%d',round(p.R02*10),round(p.s*100)),'.tif');
    saveas(f3,figname1,'tif')
    end
    pause(1)
    
end


%% Threshold calculation

R01 = p.beta1/(p.mu+1/p.tauI1);
R02 = p.beta2/(p.mu+1/p.tauI2);

cutoff = (1-1/R02)/(1-1/R01);
disp(cutoff)

%% ODEs function

function dy = sirb(~,y,p)

s1=y(1); i1=y(2); r1=y(3); f1=y(4);
s2=y(5); i2=y(6); r2=y(7); f2=y(8);

if strcmp(p.type,'none')
    e1 = exp(-p.k*f1);
    e2 = exp(-p.k*f2);
    beta1 = e1*p.beta1;
    beta2 = e2*p.beta2;
elseif strcmp(p.type,'perfect')
    e1 = exp(-p.k*p.w1*f1);
    e2 = exp(-p.k*p.w2*f2);
    beta1 = e1*e2*p.beta1;
    beta2 = e1*e2*p.beta2;
elseif strcmp(p.type,'imperfect')
    e1 = exp(-p.k*p.w1*f1);
    e2 = exp(-p.k*p.w2*f2);
    beta1 = e1*(1-p.s*(1-e2))*p.beta1;
    beta2 = (1-p.s*(1-e1))*e2*p.beta2;
else
    disp('Error')
    return
end

% Differential equations for S1, I1, R1, F1, S2, I2, R2, F2
dy(1,1) = p.mu-beta1*s1*i1-p.mu*s1+r1*p.tauR;
dy(2,1) = beta1*s1*i1-i1/p.tauI1-p.mu*i1;
dy(3,1) = i1/p.tauI1-p.mu*r1-r1*p.tauR;
dy(4,1) = (i1-f1)/p.tauF;
dy(5,1) = p.mu-beta2*s2*i2-p.mu*s2+r2*p.tauR;
dy(6,1) = beta2*s2*i2-i2/p.tauI2-p.mu*i2;
dy(7,1) = i2/p.tauI2-p.mu*r2-r2*p.tauR;
dy(8,1) = (i2-f2)/p.tauF;

end

%% Initial Condtions function

function IC = InitCond()
%%%%%% Initial Conditions

I01 = 0.0001;
S01 = 1-I01;
R01 = 0;
F01 = 0;
I02 = 0.0001;
S02 = 1-I02;
R02 = 0;
F02 = 0;

IC = [S01; I01; R01; F01; S02; I02; R02; F02];

end


%% Effective Reproductive Number

function [re1, re2] = EffRep(y,p)

s1 = y(:,1); f1 = y(:,4);
s2 = y(:,5); f2 = y(:,8);

if strcmp(p.type,'none')
    e1 = exp(-p.k*f1);
    e2 = exp(-p.k*f2);
    beta1 = e1*p.beta1;
    beta2 = e2*p.beta2;
elseif strcmp(p.type,'perfect')
    e1 = exp(-p.k*p.w1*f1);
    e2 = exp(-p.k*p.w2*f2);
    beta1 = e1.*e2*p.beta1;
    beta2 = e1.*e2*p.beta2;
elseif strcmp(p.type,'imperfect')   
    e1 = exp(-p.k*p.w1*f1);
    e2 = exp(-p.k*p.w2*f2);
    beta1 = e1.*(1-p.s*(1-e2))*p.beta1;
    beta2 = (1-p.s*(1-e1)).*e2*p.beta2;
else
    disp('Error')
    return
end


re1 = beta1.*s1*p.tauI1;
re2 = beta2.*s2*p.tauI2;

end
