%% Main Script for MaTH_CAMPS
% Runs everything with Euler timestepping
% Nicholas Knezek, mainly.
% Matt Weller, cool dude.
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

%%
Myr = pm.Myr;

time_start = 0; % [Myr] start time
time_init = 1; % [Myr] length of initial burn-in time
dt_myr_init = 0.001; % [Myr] size of initial burn-in Euler timestep
time_end = 4500;   % [Myr] How long to run time series
dt_myr = 0.01; % [Myr] size of Euler timestep for most of Earth history

Nt_init = round(time_init/dt_myr_init);
times_init = linspace(0,time_init*Myr,Nt_init);

Nt = round((time_end-time_init)/dt_myr);
dt = dt_myr*Myr;
pc.dt = dt;
times = [times_init,linspace(times_init(end)+dt, time_end*Myr, Nt)];

Ntkeep_approx = 1000; % approx number of timesteps to keep
dtkeep = max(1,ceil(Nt/Ntkeep_approx));
Ntkeep = ceil(Nt/dtkeep);
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
T0(1:3) = Tm0;
%% Solve the system using Euler stepping
T = T0;
Tvec = zeros(Ntkeep,length(T0));
tvec = zeros(Ntkeep,1);
fvec = zeros(Ntkeep,pm.n);

i = 1;
Tvec(i,:) = T;
tvec(i) = times(i);
fvec(i,:) = f0;
i = i+1;
for it=2:Nt
    dTdt = convectionODE(times(it),T',pm,pc);
    T = T + dTdt'*(times(it)-times(it-1));
    % Melting processing
    % calculate change in temperature from melt and update Ta
    [f,dTm] = mantle.melt.melt(T(1:pm.n),pm);
    T(1:pm.n) = T(1:pm.n)-dTm;

    if mod(it,dtkeep) == 0
        Tvec(i,:) = T;
        tvec(i) = times(it);
        fvec(i,:) = f;
        i = i+1;
    end
end
Tm = Tvec(:,1:pm.n);
Tc = Tvec(:,pm.n+1:end);

%% Plot things to test
hold on
core.plot.temp_v_r(Tc(1,:),pc)
core.plot.temp_v_r(Tc(47,:),pc)
core.plot.temp_v_r(Tc(114,:),pc)
core.plot.temp_v_r(Tc(225,:),pc)
core.plot.temp_v_r(Tc(end,:),pc)
legend('init','198Myr','500Myr','1Byr','final')
title('Core Temperature')
grid on
figure()
hold on
title('Temp evolution')
plot(tvec/Myr,Tm(:,1))
plot(tvec/Myr,Tm(:,2))
plot(tvec/Myr,Tm(:,3))
plot(tvec/Myr,Tc(:,end))
plot(tvec/Myr,Tc(:,1))
xlabel('time')
ylabel('Temperature (K)')
legend('lithosphere','mantle','lower mantle layer','cmb','core center')
grid on
hold off


%% Postprocessing
pp = mantle.post_processing(t,Ta,f,pm,pc);


%% Plot all things in plotting function
% Npl = 1000;
% plotting(t(1:Npl:end),Myr,Qu(1:Npl:end,:),n,Ht(1:Npl:end,:),A, ...
% Ta(1:Npl:end,:),flv(1:Npl:end,:),Tmat(1:Npl:end,:),Rmat(1:Npl:end,:), ...
% Ra(1:Npl:end,:),Crt(1:Npl:end,:),f(1:Npl:end,:),B(1:Npl:end,:),Tb(1:Npl:end,:))


%% Dynamo shutoff time based on termination heat output of ~ 1-2 Tw -- 20-60 mw/m^2
% core.utils.dynamo_timing(Qu)

