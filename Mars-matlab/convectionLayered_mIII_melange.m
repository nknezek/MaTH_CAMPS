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
n = pm.n;
time_end=4;   %How long to run time series in Myrs
Myr = pm.Myr;

%% Initial Temperatures

% starting temperatures of each layer
%Ta has time series of average temperature in each layer
Tm0 = [pm.Tliq(1), pm.Tliq(2), pm.Tsol(3)-250];
Tc0 = core.utils.adiabat(pm.Tsol(4)+25, pc);
T0 = [Tm0,Tc0];

%% Early melt processing for solid state convection 

% Calculate melt from initial conditions to set initial solid state
% temperature field for ODE solver
[f0,dTm0] = mantle.melt.melt_ini(Tm0,pm);
[dTli,flmi,flvi,Crti,fli] = mantle.melt.migration_ini(f0,Tm0(1),pm);   % calculates inital melting in lid
f0(1) = fli(1);    % induced melt from crystaliztion in lid

% Set initial temperatures to solidus approximation (solid state convection)
for x=1:pm.n
    Tm0(x) = Tm0(x)-dTm0(x);
end

% Solve the system the first time, before melt processing
[t,Tall] = ode23s(@(t,Tall)convectionODE(t,Tall,pm,pc), [0/Myr, time_end*Myr], T0');
Tm = Tall(:,1:pm.n);
Tc = Tall(:,pm.n+1:end);
%% Melting processing
% loops and solves equations over time frame of melt generation

% calculate change in temperature from melt and update Ta
[f,dTm] = mantle.melt.melt(Tm,pm);
Tm = Tm-dTm;
% Splice melt fractions for different steps, corrects because ODE solver
% doesn't include melt, so corrects
for x = 1:n
    f(1,x) = f(1,x)+f0(1,x);
end

% if f is generated from melt in time step (j, where j>t1), use time step j as next Tm0 (solidus temp) in solver.  
% workflow --> ode -> to melt -> is f generated? -> yes, take time step j
% from f, use as new Tm0 in solver, resolve -> no, continue as is. 

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
% peasy --> eg., [t(10)/Myr,time_end*Myr],Tm0);, where Tm0 would now be the last
% solved f step, then splice times togther -- WATCH INDEXING
Tall = [Tm,Tc];
if length(ind) > 1
    for m=2:length(tmelt)
        [t1,Tall1]=ode45(@(t,T)convectionODE(t,T,pm,pc), [tmelt(m)*Myr,time_end*Myr], Tall(m-1,:)');  % uses Temperatures from last time step Ta(m-1) as inital condition for new melt time step tmelt(m)
        Tm1 = Tall1(:,1:pm.n);
         [f1,dTm1] = mantle.melt.melt(Tm1,pm);
         Tm1 = Tm1-dTm1;
         Tall1 = [Tm1,Tall1(:,pm.n+1:end)];
          % Modify t and Ta such that the new output t is spliced into orginal t
          t=[t(1:m-1);t1(2:end)];                      % new time vector running from last melt instance
          Tall=[Tall(1:m-1,:);Tall1(2:end,:)];               % new Temperature matrix running from last melt instance
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

