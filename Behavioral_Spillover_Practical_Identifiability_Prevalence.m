% This code uses the Monte Carlo approach to perform practical 
% identifiability analysis. 

close all
clear 

tic 

%%%%%% Standard parameters
p.tauP = 30;
p.k = 100;
p.tauI  =7;
%p.tauI2  = p.tauI1;
p.w1 = 1;
p.mu = 0; %11/1000*1/365; % natural birth/death
p.tauR1 = 0.01; % waning period

% if p.tauR~=0
%     p.tauR1 = 1/p.tauR;
% else
%     p.tauR1 = 0;
% end


% EXPERIMENT 1 Independent 
%%%%%% Changing parameters 
p.type = 'step1';
p.R01 = 3;
p.R02 = 2;
p.w2 = 1;

% EXPERIMENT 2 Perfect spillover
%%%%%% Changing parameters 
% p.type = 'step2';
% p.R01 = 3;
% p.R02 = 2;
% p.w2 = 1;

% EXPERIMENT 3 Imperfect spillover
%%%%%% Parameters
% p.type = 'step3';
% p.R01 = 3;
% p.R02 = 2;
% p.w2 = 1;
% p.s = 0.1;


p.beta1  = p.R01/p.tauI;
p.beta2  = p.R02/p.tauI;

%%%%%% Initial Conditions
IC = InitCond();

%%%%%% Solvers and solutions
p.tmax = 1825;

tspan=1:1:p.tmax; 

[t,y_solall] = ode45(@sirb, tspan, IC, [], p);

tspanData=1:5:p.tmax; 
y_solA=y_solall(tspanData,2); 
y_solB=y_solall(tspanData,6);

IAdata=y_solA';
IBdata=y_solB'; 

%Number of iterations for the simulation. 
NumberofIterations=100; 

Fitted_Parameters=[p.beta1 p.beta2 p.k p.tauR1 p.tauI p.tauP]; %betaA, betaB, k, tauR, tauI, tauP

Number_Parameters=length(Fitted_Parameters); %Number of parameters. 
ARE = zeros(6,Number_Parameters); %Storage for Average Relative Estimation error.  

EstiParam = zeros(Number_Parameters,NumberofIterations, 6); %Storage for estimated parameters. 
Levels = [0, 0.01, 0.05, 0.1, 0.2, 0.3]; %Noise levels. 


%% Fitting Check 
% %Bounds for the parameters
% Lowerbounds = [0 0 0 0 0 0];
% Upperbounds=[1 1 120 1 20 50];
% 
% %Optimization 
% options=optimset('Disp','off','TolX',1e-10,'TolFun',1e-10,'MaxIter',10000,'MaxFunEval',10000); 
% [EstimatedParameters,fval,exitflag]=fminsearchbnd(@(Fitted_Parameters)error_independent(tspan, IAdata, IBdata, IC, Fitted_Parameters, p), Fitted_Parameters, Lowerbounds, Upperbounds, options); 
% 
% p.beta1 = EstimatedParameters(1);
% p.beta2 = EstimatedParameters(2);
% p.k = EstimatedParameters(3);
% p.tauR1 = EstimatedParameters(4);
% p.tauI = EstimatedParameters(5);
% p.tauP = EstimatedParameters(6);
% 
% [t,y_solEst] = ode45(@sirb, tspan, IC, [], p);
% 
% figure
% plot(tspan, y_solEst(:,2), 'LineWidth',2)
% hold on
% plot(tspanData, y_solA, '.r', 'MarkerSize',14)
% 
% 
% figure
% plot(tspan, y_solEst(:,6), 'LineWidth',2)
% hold on
% plot(tspanData, y_solB, '.r', 'MarkerSize',14)

%% MC Optimization 
for IterationLevels = 1:6

NoiseLevel = Levels(IterationLevels); %Defines noise level for current iteration. 

parfor i= 1:NumberofIterations

%For IAdata
NoiseA  = NoiseLevel*y_solA; 
IAdata = normrnd(y_solA, NoiseA)'; 

checkA(:,i, IterationLevels)=IAdata; %We can use this and compare it to IADataest. 
valueA=sum(IAdata(:)<0); %Inside the sum is a logical operator (0 or 1).  If the sum of these are non-zero, we enter the while-loop. Otherwise, proceed as normal. 

while (valueA~=0) %This while-loop will produce a positive data set. 
    IAdata = normrnd(y_solA, NoiseA)';
    valueA=sum(IAdata(:)<0);
end

%For IBdata
NoiseB  = NoiseLevel*y_solB; 
IBdata = normrnd(y_solB, NoiseB)'; 

checkB(:,i, IterationLevels)=IBdata; %We can use this and compare it to IBDataest. 
valueB=sum(IBdata(:)<0); %Inside the sum is a logical operator (0 or 1).  If the sum of these are non-zero, we enter the while-loop. Otherwise, proceed as normal. 

while (valueB~=0) %This while-loop will produce a positive data set. 
    IBdata = normrnd(y_solB, NoiseB)';
    valueB=sum(IBdata(:)<0);
end


%Bounds for the parameters
Lowerbounds = [0 0 0 0 0 0];
Upperbounds=[1 1 120 1 20 50];

%Optimization  
options=optimset('Disp','off','TolX',1e-10,'TolFun',1e-10,'MaxIter',10000,'MaxFunEval',10000); 
[EstimatedParameters,fval,exitflag]=fminsearchbnd(@(Fitted_Parameters)error_independent(tspan, IAdata, IBdata, IC, Fitted_Parameters, p), Fitted_Parameters, Lowerbounds, Upperbounds, options); 

EstiParams(:,i) = EstimatedParameters'; %Stores parameters for current noise level to compute ARE.  
EstiParam(:,i, IterationLevels) = EstimatedParameters'; %Stores estimated parameters for each noise level
ExitFlag(:,i, IterationLevels) = exitflag;
Fval(:,i, IterationLevels)=fval; 

IADataest(:,i, IterationLevels) = IAdata';
IBDataest(:,i, IterationLevels) = IBdata';
end

%Computes the ARE Score
ARE_Value = zeros(1,Number_Parameters);  %Storage for ARE calculation. 
    for i = 1:Number_Parameters
        ARE_Value(i) = (100/NumberofIterations) * sum(abs(Fitted_Parameters(i) - EstiParams(i,:)))/abs(Fitted_Parameters(i));
    end

    ARE(IterationLevels,:) = ARE_Value;
end

toc


%% ODEs function 

function dy = sirb(~,y,p)

s1=y(1); i1=y(2); r1=y(3); f1=y(4);
s2=y(5); i2=y(6); r2=y(7); f2=y(8);

if strcmp(p.type,'step1')
    e1 = exp(-p.k*f1);
    e2 = exp(-p.k*f2);
    beta1 = e1*p.beta1;
    beta2 = e2*p.beta2;
elseif strcmp(p.type,'step2')
    e1 = exp(-p.k*p.w1*f1);
    e2 = exp(-p.k*p.w2*f2);
    beta1 = e1*e2*p.beta1;
    beta2 = e1*e2*p.beta2;
elseif strcmp(p.type,'step3')   
    e1 = exp(-p.k*p.w1*f1);
    e2 = exp(-p.k*p.w2*f2);
    beta1 = e1*(1-p.s*(1-e2))*p.beta1;
    beta2 = (1-p.s*(1-e1))*e2*p.beta2;
else
    disp('Error')
    return
end



% Differential equations for S1, I1, R1, F1, S2, I2, R2, F2
dy(1,1) = p.mu-beta1*s1*i1-p.mu*s1+r1*p.tauR1;
dy(2,1) = beta1*s1*i1-i1/p.tauI-p.mu*i1;
dy(3,1) = i1/p.tauI-p.mu*r1-r1*p.tauR1;
dy(4,1) = (i1-f1)/p.tauP;
dy(5,1) = p.mu-beta2*s2*i2-p.mu*s2+r2*p.tauR1;
dy(6,1) = beta2*s2*i2-i2/p.tauI-p.mu*i2;
dy(7,1) = i2/p.tauI-p.mu*r2-r2*p.tauR1;
dy(8,1) = (i2-f2)/p.tauP;

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



%% Sum of Square Errors for Independent Virus

function SSE_Independent= error_independent(tspan, IAdata, IBdata, IC, Fitted_Parameters, p)

p.beta1 = Fitted_Parameters(1);
p.beta2 = Fitted_Parameters(2);
p.k = Fitted_Parameters(3);
p.tauR1 = Fitted_Parameters(4);
p.tauI = Fitted_Parameters(5);
p.tauP = Fitted_Parameters(6);

[~,y] = ode45(@sirb, tspan, IC, [], p);

tspanData=1:5:max(tspan); 

IA=y(tspanData,2)';
IB=y(tspanData,6)';

SSE_Independent = sum((IA- IAdata).^2+(IB-IBdata).^2);
end





