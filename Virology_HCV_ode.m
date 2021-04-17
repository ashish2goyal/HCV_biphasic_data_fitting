%% the ODE equation structure

function dydt= Virology_HCV_ode(t,y) % Virology_HCV_ode will be called in an odesolver

% parameters in the model
global nt; % represents reduction in infectivity
global bt; % represents viral infectivity
global T_0; % represents the constant number/concentration of target cells
global del; % represents the death rate 
global ep; % represents efficacy of the drug in inhibiting the viral production
global p; % representing viral production
global c; % representing viral clearance

%% ODE MODEL

dydt=zeros(2,1);

dydt(1) = (1-nt)*bt*T_0*y(2)-del*y(1); % HCV infected cells (I equation; y(1))

dydt(2) = (1-ep)*p*y(1) - c*y(2); % dually-infected cells (V equation; y(2))
