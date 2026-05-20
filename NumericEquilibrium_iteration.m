%% Analytical calculation of Co-existence threshold

clearvars; close all; clc
savedata = 1;
tic

%% Initialize variables and parameters symbolically

syms S1 I1 R1 F1 S2 I2 R2 F2
syms k w1 w2 s b1 b2 tauR tauI1 tauI2 tauF

e1 = exp(-k*w1*F1);
e2 = exp(-k*w2*F2);
beta1 = e1*(1 - s*(1 - e2))*b1;
beta2 = (1 - s*(1 - e1))*e2*b2;
dS1 = - beta1*S1*I1 + tauR*R1;
dI1 = beta1*S1*I1 - I1/tauI1;
dR1 = I1/tauI1 - tauR*R1;
dF1 = (I1 - F1)/tauF;
dS2 = - beta2*S2*I2 + tauR*R2;
dI2 = beta2*S2*I2 - I2/tauI2;
dR2 = I2/tauI2 - tauR*R2;
dF2 = (I2 - F2)/tauF;

%% Solve symbolically a subsystem dependent on I1 and I2

out = solve(dS1==0,dI1==0,dF1==0,dS2==0,dI2==0,dF2==0,S1,R1,F1,S2,R2,F2);

%% Introduce parameter values to subsystem solution


R02_vec = linspace(1,3,101);
sin_vec = linspace(0,1,101);
thres = NaN(length(R02_vec),length(sin_vec));
R0out = NaN(length(R02_vec),length(sin_vec));
sout = NaN(length(R02_vec),length(sin_vec));
for s_index = 1:length(sin_vec)
    sin = sin_vec(s_index);
    for R0_index = 1:length(R02_vec)
        p.R02 = R02_vec(R0_index);

        p.tauF = 30;
        p.k = 100;
        p.tauI1  = 7;
        p.tauI2  = p.tauI1;
        p.w1 = 1;
        p.w2 = p.w1;
        p.tauR = 1/100; % waning
        p.R01 = 3;

        p.b1  = p.R01/p.tauI1;
        p.b2  = p.R02/p.tauI2;

        p.s = sin;

        out1 = subs(out,{tauF,k,tauI1,tauI2,w1,w2,tauR,b1,b2,s},...
            [p.tauF,p.k,p.tauI1,p.tauI2,p.w1,p.w2,p.tauR,p.b1,p.b2,sin]);

        %% Numerically solve for I1 and I2 and plug in
        outI = vpasolve([1==out1.S1+I1+out1.R1,1==out1.S2+I2+out1.R2],[I1,I2]);

        q.S1 = subs(out1.S1,{I1,I2},[outI.I1,outI.I2]);
        q.I1 = outI.I1;
        q.R1 = subs(out1.R1,{I1,I2},[outI.I1,outI.I2]);
        q.F1 = subs(out1.F1,{I1,I2},[outI.I1,outI.I2]);
        q.S2 = subs(out1.S2,{I1,I2},[outI.I1,outI.I2]);
        q.I2 = outI.I2;
        q.R2 = subs(out1.R2,{I1,I2},[outI.I1,outI.I2]);
        q.F2 = subs(out1.F2,{I1,I2},[outI.I1,outI.I2]);

        %%
        thresA = exp(p.k*q.I2)/(exp(p.k*q.I2)-p.s*(exp(p.k*q.I2)-1));
        thresB = exp(p.k*q.I1)/(exp(p.k*q.I1)-p.s*(exp(p.k*q.I1)-1));

        thres(R0_index,s_index) = double(thresB);
        R0out(R0_index,s_index) = p.R02;
        sout(R0_index,s_index) = p.s;


    end

end
toc
%% Save data
if savedata ==1
    save('Data_Numeric_Equilibrium.mat')
end

