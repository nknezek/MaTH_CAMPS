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

% wtpS = 5; % 0 - 25 wt% S supported
% param_case = 2; % 0: no layer, 1: hot case, 2: cold case

for wtpS = [5,15,25]
for param_case = 1:2
    
folder_casenames = ["nolayer/","hot/","cold/"];
basefolder = './results/';
folder = [char(folder_casenames(param_case+1)),sprintf('%.0fwtpS/', wtpS)];
mkdir([basefolder,folder])

pm = mantle.parameters(param_case, wtpS); % mantle parameters

pc = core.parameters(pm); % core parameters
n = pm.n;
%%
Myr = pm.Myr;

time_start = 0; % [Myr] start time
time_init = 1; % [Myr] length of initial burn-in time
dt_myr_init = 0.001; % [Myr] size of initial burn-in Euler timestep
% time_end = 500;   % [Myr] How long to run time series % FIX THIS WRONG TODO
time_end = 4600;   % [Myr] How long to run time series
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


% Initial Temperatures

% starting temperatures of each layer
%Ta has time series of average temperature in each layer
Tm0 = [pm.Tliq(1), pm.Tliq(2), pm.Tsol(3)-250];
Tc0 = core.utils.adiabat(pm.Tsol(4)+25, pc);
T0 = [Tm0,Tc0];


% Early melt processing for solid state convection 

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


% Solve the system using Euler stepping
T = T0;
Tvec = zeros(Ntkeep,length(T0));
tvec = zeros(Ntkeep,1);
fvec = zeros(Ntkeep,pm.n);
melt_mass_cumulative = zeros(1,pm.n);
melt_mass_vec = zeros(Ntkeep,pm.n);
crust_thickness_vec = zeros(Ntkeep,pm.n);
crust_thickness_cumulative = 0;
dT_lid = 0;
melt_mass_lid = 0;
melt_volume_lid = 0;
thickness_lid = 0;

i = 1;
Tvec(i,:) = T;
tvec(i) = times(i);
fvec(i,:) = f0;
i = i+1;
for it=2:Nt
    dTdt = convectionODE(times(it),T',pm,pc);
    Tprev = T;
    T = T + dTdt'*(times(it)-times(it-1));
    % Melting processing
    % calculate change in temperature from melt and update Ta
    [f,dTm] = mantle.melt.melt(Tprev(1:pm.n),T(1:pm.n),pm);
    T(1:pm.n) = T(1:pm.n)-dTm;
    melt_mass_cumulative = melt_mass_cumulative + f*pm.M;
    
    % melt migrates to lid and changes lid temp
    [dT_lid, melt_mass_lid, melt_volume_lid, thickness_lid] = mantle.melt.migration(f,T,pm);
%     T(1) = T(1)+dT_lid;
    crust_thickness_cumulative = crust_thickness_cumulative + thickness_lid;
    
    if mod(it,dtkeep) == 0
        Tvec(i,:) = T;
        tvec(i) = times(it);
        fvec(i,:) = f;
        melt_mass_vec(i,:) = melt_mass_cumulative;
        crust_thickness_vec(i,:) = crust_thickness_cumulative;        
        i = i+1;
    end
end
Tm = Tvec(:,1:pm.n);
Tc = Tvec(:,pm.n+1:end);


% Postprocessing
pp = mantle.post_processing(tvec,Tvec,fvec,pm,pc);
pp.flv = melt_mass_vec;
pp.Crt = crust_thickness_vec;
pp.pm = pm;
pp.pc = pc;
save([basefolder,folder,'pp.mat'],'pp')

%% Plot all things in plotting function

% wtpS = 15;
% param_case = 2;
for wtpS = [5,15,25]
for param_case = 1:2
folder_casenames = ["nolayer/","hot/","cold/"];
basefolder = './results/';
folder = [char(folder_casenames(param_case+1)),sprintf('%.0fwtpS/', wtpS)];
mkdir([basefolder,folder])

load([basefolder,folder,'pp.mat'])
pm = mantle.parameters(param_case, wtpS); % mantle parameters
pc = core.parameters(pm); % core parameters
pp.pm = pm;
pp.pc = pc;
save([basefolder,folder,'pp.mat'],'pp')

fig = mantle.plot.run_summary(pp,pm,pc);
saveas(fig,[basefolder,folder,'run_summary.png'])

fig = mantle.plot.Urey_Buoyancy(pp,pm,pc);
saveas(fig,[basefolder,folder,'Urey.png'])

fig = mantle.plot.melt(pp,pm,pc);
saveas(fig,[basefolder,folder,'melt.png'])

fig = core.plot.temp_layers(pp,pm,pc);
saveas(fig,[basefolder,folder,'temp_layers.png'])

fig = core.plot.temp_profiles(pp,pm,pc);
saveas(fig,[basefolder,folder,'temp_profiles.png'])


close all 
end
end
%%
end
end