%% Test Thermal Evolution Code 
% Adiabatic convection up and thermal conduction down
% Nicholas Knezek

N = 50;
R = 2000e3; % [m]
dr = R/N; % [m]
dt = 365.25*24*3600*1e6; % s (1 Myr)
D = 6203e3; % m, adiabatic lengthscale

Nt = 4000;

cp = 840; %J/kg-K = m^2/s^2-K
rho = linspace(12e4, 9e4, N); % kg/m^3
k = 135; % J/m-K-s = kg-m/K-s^3

re = linspace(0,R,N+1); % m
r = linspace(dr/2, R-dr/2, N); %m
dV = 4*pi*dr*r.^2; % m^3
A = 4*pi*re.^2; % m^2

rho_cp_dV = rho.*cp.*dV; % J/K

T_cmb = 5000; % K
T_ad = adiabat(T_cmb, r, R, D); % K
T_p = linspace(0,0,N); % K
T_0 = T_ad;

Q_cmb_t = ones(1,Nt);
Q_cmb_t(1:500) = 2e12;
Q_cmb_t(500:1000) = -1e12;
Q_cmb_t(1000:end) = 5e12;

%%
for i=1:Nt
[T_ad, T_p] = evolve_core(T_ad, T_p, Q_cmb_t(i));
if mod(i,100)==0
    hold on
    plot(r,T_ad+T_p)
end
end
