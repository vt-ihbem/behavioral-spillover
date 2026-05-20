%% Generate data for Simulation Heatmaps 
clearvars; close all;
clc
tic
%% Standard parameters

savedata = 1;
%%%%%% Standard parameters
p.tauF = 30;
p.k = 100;
p.tauI1  = 7;
p.tauI2  = p.tauI1;
p.w1 = 1;
p.mu = 0; % natural birth/death
p.tauR = 1/100; % waning immunity


%% Two interdependent viruses with spillover

%%%%%% Parameters
p.type = 'imperfect';
p.w2 = 1;

p.R01 = 3;
p.beta1  = p.R01/p.tauI1;

R02vec = linspace(1,3,501);
svec = linspace(0,1,501);

%%%%%% Initial Conditions
IC = InitCond();

%%%%%% Solvers and solutions

tmax = 365*1;

reout1 = NaN(length(R02vec),length(svec));
reout2 = NaN(length(R02vec),length(svec));
iout1 = NaN(length(R02vec),length(svec));
iout2 = NaN(length(R02vec),length(svec));
Bdom = NaN(length(R02vec),length(svec));
Bdom1 = NaN(length(R02vec),length(svec));
Bdom1a = NaN(length(R02vec),length(svec));
Bdom2 = NaN(length(R02vec),length(svec));
cutoff = NaN(length(R02vec),1);
options = odeset('NonNegative',1:8);

for r0_it = 1:length(R02vec)
    p.R02 = R02vec(r0_it);
    p.beta2  = p.R02/p.tauI2;

    for s_it = 1:length(svec)
        p.s = svec(s_it);
        
        [t,y] = ode15s(@sirb, 0:tmax, IC, options, p);
        i1 = y(:,2); i2 = y(:,6);
        if sum(i2(i2>i1))==0
            Bdom(r0_it,s_it) = 0;
        elseif sum((i2(i2>i1)-i1(i2>i1)))<1e-9
            Bdom(r0_it,s_it) = NaN;
        else
            Bdom(r0_it,s_it) = sum((i2-i1)>0)/tmax;
        end
        Bdom1(r0_it,s_it) = sum((i2(i2>i1)-i1(i2>i1)));
        Bdom2(r0_it,s_it) = sum((i2-i1)>0)/tmax;
        Bdom1a(r0_it,s_it) = sum((i1(i1>i2)-i2(i1>i2)));
        [re1, re2] = EffRep(y,p);

        reout1(r0_it,s_it) = re1(end);
        reout2(r0_it,s_it) = re2(end);
        iout1(r0_it,s_it) = i1(end);
        iout2(r0_it,s_it) = i2(end);

    end
    % Threshold calculation
    R01 = p.beta1/(p.mu+1/p.tauI1+p.tauR);
    R02 = p.beta2/(p.mu+1/p.tauI2+p.tauR);

    cutoff(r0_it) = (1-1/R02)/(1-1/R01);

end

if p.tauR == 0
    if p.mu == 0
        filename = 'Data_mu0_tauR0.mat';
    else
        filename = sprintf('Data_mu%d_tauR0.mat',round(1/p.mu));
    end
else
    if p.mu == 0
        filename = sprintf('Data_mu0_tauR%d.mat',round(1/p.tauR));
    else
        filename = sprintf('Data_mu%d_tauR%d.mat',round(1/p.mu),round(1/p.tauR));
    end
end

%% Save data
filename = 'Data_mu0_tauR100_long.mat';
if savedata == 1
    save(filename)
end

toc
return




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

dy = NaN(8,1);

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

e1 = exp(-p.k*p.w1*f1);
e2 = exp(-p.k*p.w2*f2);
beta1 = e1.*(1-p.s*(1-e2))*p.beta1;
beta2 = (1-p.s*(1-e1)).*e2*p.beta2;
re1 = beta1.*s1*p.tauI1;
re2 = beta2.*s2*p.tauI2;

end
