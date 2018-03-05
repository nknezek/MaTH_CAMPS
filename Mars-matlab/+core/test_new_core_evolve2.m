%% Test Thermal Evolution Code 
%
% Adiabatic convection up and thermal conduction down
% Nicholas Knezek

dt = 365.25*24*3600*1e6; % s (1 Myr)
Nt = 4500;

T_cmb = 2300; % K
T_ad = core.utils.adiabat(T_cmb); % K

Q_cmb_t = ones(1,Nt);
Q_cmb_t(1:500) = linspace(8e12, -0.25e12,500);
Q_cmb_t(501:1500) = linspace(-0.25e12,0.75e12,1000);
Q_cmb_t(1501:end) = linspace(0.75e12,0.25e12,Nt-1500);

%%
T = T_ad;
T_cmb = zeros(1,Nt);
T_cen = zeros(1,Nt);
for it=1:Nt
%     T = convect_and_conduct(T, r, dr, D, A, dt, k, rho_cp_dV, Q_cmb_t(it));
    dTdt = dconvect_and_conduct_dt(T, r, dr, D, A, dt, k, rho_cp_dV, Q_cmb_t(it));
    T = T+dTdt*dt;
    T_cmb(it) = T(end);
    T_cen(it) = T(1);
end
%%
plot(r,T_ad)
hold on
plot(r,T)
figure()
plot(T_cen)
hold on
plot(T_cmb)

%%
E_i = energy(T_ad,rho_cp_dV);
E_o = energy(T,rho_cp_dV);
dEt = E_o-E_i
dEq = -sum(Q_cmb_t(1:Nt)*dt)
%%
figure()
plot(Q_cmb_t(1:Nt))
hold on
plot(Q_ad_c_out(1:Nt))

figure()
plot(T_cmb)
hold on
plot(T_cen)

figure()
plot(E_core)

figure()
plot(dR_TBL)
%% Test ode solver
Myr = 365.25*24*3600*1e5; % s (1 Myr)
dt = 1*Myr;

odef = @(t,T) dconvect_and_conduct_dt(T', r, dr, D, A, dt, k, rho_cp_dV, 1e12)';
% odef = @(t,y) dTdt_core(y, Q_cmb);
tspan = linspace(1,4500,4500)*Myr;
[t, T_t] = ode15s(odef, tspan, T_ad');
%%
plot(t,T_t(:,end))
hold on
plot(t,T_t(:,1))
legend('cmb','center')

figure()
plot(r,T_t(1,:))
hold on
plot(r,T_t(int16(length(t)/2),:))
plot(r,T_t(end,:))
legend('init','mid','final')

