%% Main function estimating parameters, simulating and plotting ODE equations
% to note that is only for datasets excluding censored datasets
function HCV_main

% define all variables that are used by other scripts as global (for example, here we call Virology_HCV_ode.m) 
global nt;
global bt;
global T_0;
global del;
global ep;
global p;
global c;
global tCm CV; % make observed data global for optimization algorithm in other script to access (Virology_opt_HCV.m)

% define parameter values
nt=0.5; % assumed based on literature
bt=7*10^-9 ; %assumed based on literature
T_0=10^7; % represents the constant number/concentration of target cells including at time t=0
% del=0.2; % represents the death rate of infected cells -- will be estimated
% ep=0.998; % represents efficacy of the drug in inhibiting the viral production -- will be estimated
% p=10^3; % representing viral production rate -- will be estimatd
% c=23; % representing viral clearance -- assumed based on literature  -- will be estimated

%% Define observed data here  
tCm=[0 1/24 2/24 4/24 8/24 12/24 1 36/24 3 5 7 10 14 21]'; % measured time points
CV = 1.0e+05 * [1.3356  1.2423 1.6559  1.4978  1.3236  0.1653    0.0079    0.0023   0.0016  0.0009  0.0012    0.0006    0.0006    0.0005]'; % observed viral loads


%% Optimization algorithm starts here to estimate three parameters c, del, ep
% optimization algorithm used -- either 'trust-region-reflective' (default) or 'levenberg-marquardt'.

% (we repeat the search with 100 different initial guesses to avoid local minima)

i=1
while(i<101)
%% pre-define vector of size equal to the number of unknown parameters that will be passed to optimization algorithm
v=zeros(3,1); 

%% our initial guesses
del=0.1 + 0.3* rand(1); % random initial guess for del
c = 6 + 23* rand(1);  % random initial guess for c
ep = 0.9 + 0.09*rand(1);  % random initial guess for ep

% Using following, we first transform our parameters to ensure that c and del are positive whereas ep is between 0 and 1 
v(1)=log(del);
v(2)=log(c);
v(3)=log(ep/(1-ep));
% we will be optimizing v(1), v(2) and v(3) instead of del, c and ep by passing v to optimization routine below

%% Optimization routine
Options=optimset('Display','iter','Algorithm','levenberg-marquardt'); % Choose between 'trust-region-reflective' (default) and 'levenberg-marquardt' -- for more options on tolerance, see Mathworks webpage on lsqnonlin
[v,resnorm]=lsqnonlin('Virology_opt_HCV',v,[],[],Options); % returns optimized v and  the value of the squared 2-norm of the residual at x: sum(fun(x).^2).

%% Optimized Parameters

if(isreal(v)) % only if optimization algorithm return real v, we accept the solution and store it
del_instance(i,1)=exp(v(1));
c_instance(i,1)=exp(v(2));
ep_instance(i,1)=exp(v(3))/(1+exp(v(3)));
resnorm_instance(i,1) = resnorm ; 

% m is the number of unknown parameters and n is the number of data points used in the fits
m=3;
n=length(CV)-1;
AIC_instance(i,1) = (n*log(resnorm/n)) + (2*m*n/(n-m-1)) ; % Corrected AIC

i=i+1;
end

end

%% find the parameter set that best fits the data
[val,ind] =  min(resnorm_instance) ; % this will be the one with the lowest resnorm

c=c_instance(ind);
del=del_instance(ind);
ep=ep_instance(ind);

% Since the system is assumed to be in steady state at t=0 (the time of start of antiviral treatment)
% at which we also are often aware of viral loads (V_0), we additionally have the following constraints
I_0=bt*CV(1)*T_0/del;
p=c*del/(T_0*bt); % infectivity HDV (Guedj et al.)

xinit=[I_0;CV(1)]; 
tspan=[0:max(tCm)];

%% Simulation for optimized parameters
[t,y]=ode23s(@Virology_HCV_ode,tspan,xinit);
figure(1)
semilogy(tCm,CV,'s','MarkerSize',10,'MarkerFaceColor','black')
hold on
semilogy(t,y(:,2),'Color','red','LineStyle','-','LineWidth',1)  % plot viral loads
legend('Observed','Fitted/Predicted');
ylim([10^0 10^7])
xlim([0 max(tCm)+2])
ylabel(' HCV RNA Concentration')
xlabel('time (days)')



%% a simple numerical check to determine parameter identifiablity issues
% m is the number of unknown parameters and n is the number of data points used in the fits
m=3;
n=length(CV)-1;
% find sets of (c, del, ep) that fits without distinction (i.e., within AIC_lowest and AIC_lowest+2)
AIC_lowest=AIC_instance(ind);
c_without_distinction = c_instance(AIC_instance<AIC_lowest+2);
ep_without_distinction = ep_instance(AIC_instance<AIC_lowest+2);
del_without_distinction = del_instance(AIC_instance<AIC_lowest+2);

sets_without_distinction=[c_without_distinction ep_without_distinction del_without_distinction];
disp('Visual check to determine parameter identifiablity issues')
disp(sets_without_distinction)
% if there are multiple ditincts sets of (c, ep, del) or rows in the matrix
% sets_without_distinction, then there are parameter identifiablity issues
% and the user should scatter plot (c, ep), (ep, del) and (c, del) and determine correlation between parameters
% and fix one of the parameter in the combination with the strongest correlation 

figure(2)
subplot(3,1,1)
plot(c_without_distinction,ep_without_distinction,'s','MarkerSize',10,'MarkerFaceColor','red')
ylabel('ep')
xlabel('c')

subplot(3,1,2)
plot(del_without_distinction,ep_without_distinction,'s','MarkerSize',10,'MarkerFaceColor','red')
ylabel('ep')
xlabel('del')

subplot(3,1,3)
plot(del_without_distinction,c_without_distinction,'s','MarkerSize',10,'MarkerFaceColor','red')
ylabel('c')
xlabel('del')
