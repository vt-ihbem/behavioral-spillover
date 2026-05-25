% This code generates a .mat file needed to run the structural
% identifability analysis. STRIKE-GOLDD is needed to conduct the
% structural analysis. 

% This script is for no spillover and outputs Ia(t) and Ib(t). To run the
% remaining scenarios described in Table 6 of the manuscript, change the
% output equation (line 18), parameters (line 25) and ode system (line 28).

clear 

syms betaA betaB k tauR tauI tauP...
    Sa Ia Ra Iat Sb Ib Rb Ibt 

% State Variables
x=[Sa; Ia; Ra; Iat; Sb; Ib; Rb; Ibt]; 

% Output equation (y(t)=Ia(t) and y(t)=Ib(t))
h=[Ia, Ib];

% Input equations
u=[];
w=[]; 

% Parameters 
p=[betaA; betaB; k; tauR; tauI; tauP]; 

% ODE System
f=[-(exp(-k*Iat))*betaA*Sa*Ia+Ra/tauR; (exp(-k*Iat))*betaA*Sa*Ia-Ia/tauI; Ia/tauI-Ra/tauR; (Ia-Iat)/tauP;...
    -(exp(-k*Ibt))*betaB*Sb*Ib+Rb/tauR; (exp(-k*Ibt))*betaB*Sb*Ib-Ib/tauI; Ib/tauI-Rb/tauR; (Ib-Ibt)/tauP]; 


% Save as .mat file
save('No_Spillover_Prevalence_Res','x','h','u','w','p','f');










