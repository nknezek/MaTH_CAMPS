%% Test Thermal Evolution Code 
% Adiabatic convection up and thermal conduction down
% Nicholas Knezek


%% Run parameters
core.parameters
T_cmb = 2400; % K
T_ad = core.utils.adiabat(T_cmb); % K
N_myr = 305;
end_time = N_myr*pc.Myr;
Nt = 10*N_myr;
Q_cmb = 1e11;


%% Euler Method
debug = 0;
debugerr=1;
t = linspace(0, end_time, Nt);
dt = t(2)-t(1);
T_cmb = zeros(1,Nt);
T_cen = zeros(1,Nt);
T = T_ad;
T_start = T;
for i=1:Nt
    T = core.therm.convect_and_conduct(T, dt, Q_cmb, debug,debugerr);
%     dTdt = dconvect_and_conduct_dt(T, r, dr, D, A, dt, k, rho_cp_dV, Q_cmb);
%     T = T+dTdt*dt;
    T_cmb(i) = T(end);
    T_cen(i) = T(1);    
end
T_end = T;

%% Plot Results
figure()
core.plot.temp_v_t(t,T_cmb)
hold on
core.plot.temp_v_t(t,T_cen)
legend('cmb','center')

figure()
core.plot.temp_v_r(T_start)
hold on
core.plot.temp_v_r(T_end)
legend('init','end')

%% Verify Energy Balance
E_i = core.utils.energy(T_ad);
E_f = core.utils.energy(T);
dEq = Q_cmb*Nt*dt
dE_model = E_f-E_i

%% ODE Solver
N_myr = 450;
end_time = N_myr*pc.Myr;
Nt = 450;
dt = 0.1*pc.Myr;
T = T_ad;
odef = @(t,T) core.therm.dconvect_and_conduct_dt(T', dt, Q_cmb)';
% odef = @(t,y) dTdt_core(y, Q_cmb);
tspan = linspace(0,end_time,Nt);
[t, T_t] = ode113(odef, tspan, T_ad');

%% Plot ODE
T_cmb = T_t(:,end);
T_cen = T_t(:,1);
T_start = T_ad;
T_end = T_t(end,:);

figure()
core.plot.temp_v_t(t,T_cmb)
hold on
core.plot.temp_v_t(t,T_cen)
legend('cmb','center')

figure()
core.plot.temp_v_r(r,T_start)
hold on
core.plot.temp_v_r(r,T_end)
legend('init','end')

