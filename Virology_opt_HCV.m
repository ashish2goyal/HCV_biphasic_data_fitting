function f=Virology_opt_HCV(v)

global nt;
global bt;
global T_0;
global del;
global ep;
global p;
global c;
global tCm CV; % make observed data global for optimization algorithm in other script to access (Virology_opt_HCV.m)


% untransform our v vector passed to optimization algorithm to allow simulation of ode equations
% we obtain
del=exp(v(1));
c=exp(v(2));
ep=exp(v(3))/(1+exp(v(3)));

tspan1=[0:max(tCm)];

% Since the system is assumed to be in steady state at t=0 (the time of start of antiviral treatment)
% at which we also are often aware of viral loads (V_0), we additionally have the following constraints
I_0=bt*CV(1)*T_0/del;
p=c*del/(T_0*bt); % infectivity HDV (Guedj et al.)

xinit=[I_0;CV(1)]; 

%% run system of ode equations
[t,x]=ode23s(@Virology_HCV_ode,tspan1,xinit);

VSIM=x(:,2); % predicted HCV RNA within tspan1 
% compare VSIM and ALTSIM with V and ALT
iiv=find(tCm>=0 & tCm<=max(tCm));
LOGVCOMP=interp1(t,log10(VSIM),tCm(iiv)); % determine predicted HCV RNA at observed timepoints tCm and if not available, interpolate it using interp1 function
f1=log10(CV(iiv))-LOGVCOMP; % determine residual error (we log tranform first so that errors from data points on scale of 10^7 should roughly have the same weights as data points on scale of 10^3-10^4) 

f=[f1]; 
% replace Infs and Nans
iif=~isfinite(f);
f(iif)=10;