%% Main Script for MaTH_CAMPS
% Runs everything. 
% Matt Weller, cool dude.
% Nicholas Knezek, also ok.
% March 2018

%% Added Melt and Migration

%%Todo
%       
%       add    Fix melt in lid temperature change
%               temperature dependant viscosity (be aware of 1e-16 rounding errors, will crop up in temp dep visc )
%%
clear all
pm = mantle.parameters(); % mantle parameters
pc = core.parameters(pm); % core parameters

time_end=4600;   %How long to run time series in Myrs
Myr = pm.Myr;

%% Initial Temperatures

% starting temperatures of each layer
%Ta has time series of average temperature in each layer
Ta0 = [pm.Tliq(1), pm.Tliq(2), pm.Tsol(3)-250];
Tc0 = core.utils.adiabat(pm.Tsol(4)+25, pc);
T0 = [Ta0,Tc0];

%% Early melt processing for solid state convection 

% Calculate melt from initial conditions to set initial solid state
% temperature field for ODE solver
[f0,dTm0] = mantle.melt.melt_ini(Ta0,pm);
[dTli,flmi,flvi,Crti,fli] = mantle.melt.migration_ini(f0,Ta0(1),pm);   % calculates inital melting in lid
f0(1) = fli(1);    % induced melt from crystaliztion in lid

% Set initial temperatures to solidus approximation (solid state convection)
for x=1:pm.n
    Ta0(x) = Ta0(x)-dTm0(x);
end

% Solve the system the first time, before melt processing
[t,Ta] = ode45(@(t,Ta)convectionODE(t,Ta,pm,pc), [0/Myr, time_end*Myr], T0');
Ta = Ta';
%% Melting processing
% loops and solves equations over time frame of melt generation

% calculate change in temperature from melt and update Ta
[f,dTm] = mantle.melt(Ta,pm);
Ta = Ta-dTm;

% Splice melt fractions for different steps, corrects because ODE solver
% doesn't include melt, so corrects
for x = 1:n
    f(1,x) = f(1,x)+f0(1,x);
end

% if f is generated from melt in time step (j, where j>t1), use time step j as next Ta0 (solidus temp) in solver.  
% workflow --> ode -> to melt -> is f generated? -> yes, take time step j
% from f, use as new Ta0 in solver, resolve -> no, continue as is. 

ismelt = f>0;           %Finds where melt occurs 1: yes 2: no

% collapse all melt to ismelt(:,1), assumes all melt migrates to top layer
% where melting occurs
for m = 2:(n-1)
    ismelt(:,1) = ismelt(:,1)+ismelt(:,m);            
end

% time series of melt to rerun solver and melt modules - in Myr's
ind=find(ismelt(:,1)>0);
tmelt=t(ind)/Myr;

% Loop to take time step at every f, melt, resolve Ta from ode -- easy
% peasy --> eg., [t(10)/Myr,time_end*Myr],Ta0);, where Ta0 would now be the last
% solved f step, then splice times togther -- WATCH INDEXING
if length(ind) > 1
    for m=2:length(tmelt)
        [t1,Ta1]=ode45(@(t,Ta)convectionODE(t,Ta,pm,pc), [tmelt(m)*Myr,time_end*Myr], Ta(m-1,:));  % uses Temperatures from last time step Ta(m-1) as inital condition for new melt time step tmelt(m)
         [f1,dTm1] = mantle.melt(Ta1,pm);
         Ta1 = Ta1-dTm1;
          % Modify t and Ta such that the new output t is spliced into orginal t
          t=[t(1:m-1);t1(2:end)];                      % new time vector running from last melt instance
          Ta=[Ta(1:m-1,:);Ta1(2:end,:)];               % new Temperature matrix running from last melt instance
          f=[f(1:m-1,:);f1(2:end,:)];                  %new f matrix running from last melt instance    
    end
end

%% Postprocessing
% mantle.post_processing(t,Ta,f,pm,pc)


%% Plot all things in plotting function
% Npl = 1000;
% plotting(t(1:Npl:end),Myr,Qu(1:Npl:end,:),n,Ht(1:Npl:end,:),A, ...
% Ta(1:Npl:end,:),flv(1:Npl:end,:),Tmat(1:Npl:end,:),Rmat(1:Npl:end,:), ...
% Ra(1:Npl:end,:),Crt(1:Npl:end,:),f(1:Npl:end,:),B(1:Npl:end,:),Tb(1:Npl:end,:))


%% Dynamo shutoff time based on termination heat output of ~ 1-2 Tw -- 20-60 mw/m^2
% core.utils.dynamo_timing(Qu)

