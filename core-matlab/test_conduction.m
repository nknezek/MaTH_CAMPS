%% Test adiavatic temperature gradient, then conduction down above adiabat
% Nicholas Knezek
% March 2, 2018

% function test_conduction()

%%
N = 100; 
R = 2000e3; % [m]
dr = R/N; % [m]
dt = 365.25*24*3600*1e6; % s (1 Myr)

cp = 840; %J/kg-K = m^2/s^2-K
rho = linspace(12e4, 9e4, N); % kg/m^3
k = 135; % J/m-K-s = kg-m/K-s^3

re = linspace(0,R,N+1); % m
r = linspace(dr/2, R-dr/2, N); %m
dV = 4*pi*dr*r.^2; % m^3
A = 4*pi*re.^2; % m^2

T_cmb = 5000; % K
T_ad = adiabat(T_cmb, r, R); % K
T_p = linspace(0,0,N); % K
T_p(71:end) = linspace(0,100,30);

%% Plot Adiabat
plot(r,T_ad)
hold on
plot(r, T_ad+T_p)

%% Conduct heat downwards
T = T_p;
rho_cp_dV = rho.*cp.*dV;
for i=1:100
    T = conduct(T, dr, dt, k, A, rho_cp_dV);
end
plot(r,T_p)
hold on
plot(r,T)
T_p = T;
%% Pull heat up by convection
Q_cmb = 1e1; % W out of core
T_old = T_ad+T_p;
E_p = energy(T_p,rho_cp_dV);
if Q_cmb*dt >= E_p
    T_p = zeros(N);
    Q_cmb_extra = Q_cmb-E_p;
    minfun_ad  = @(T_cmb) minfun(T_cmb, Q_cmb_extra, T_ad, rho_cp_dV, dt);
    Tcmb_new = fminbnd(minfun_ad,0,T_ad(end));
    T_ad = adiabat(Tcmb_new);
else
    minfun_p  = @(T_cmb) minfun(T_cmb, Q_cmb, T_p, rho_cp_dV, dt);
    Tpcmb_new = fminbnd(minfun_p, 0, Tmax);
    T_p(T_p>Tpcmb_new) = Tpcmb_new;
end
T_new = T_ad+T_p;

T_ad_new = adiabat(T_new(end),re(2:end),R);
%%
hold on
plot(r,T_ad)
plot(r,T_ad_new)
plot(r,T_old,'--')
plot(r,T_new,'--')
ylim([5000,5400])
%%
Tmax = T(end);
Tctest = linspace(80,Tmax,N);
errs = zeros(1,N);
for i = 1:N
    errs(i) = minfun_curr(Tctest(i));
end

plot(Tctest,errs)
line([Tcmb_new, Tcmb_new], [min(errs),max(errs)])
grid on
%%
T_new = T;
T_new(T>Tcmb_new) = Tcmb_new;
plot(r, T)
hold on
plot(r,T_new)
grid on
% end