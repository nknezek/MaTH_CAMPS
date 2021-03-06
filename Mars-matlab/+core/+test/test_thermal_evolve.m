%% Test Thermal Evolution Code 
% Adiabatic convection up and thermal conduction down
% Nicholas Knezek

N = 50;

R = 2000e3; % [m]
dr = R/N; % [m]
dt = 365.25*24*3600*1e6; % s (1 Myr)
D = 6203e3; % m, adiabatic lengthscale

Nt = 4500;

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

T_vec0 = horzcat(T_p,T_cmb);

Q_cmb_t = ones(1,Nt);
Q_cmb_t(1:500) = 1e12;
Q_cmb_t(500:1000) = -1e12;
Q_cmb_t(1000:end) = 10e12;
E_ad = zeros(1,Nt);
E_p = zeros(1,Nt);
T_vec_all = zeros(Nt, N+1);
%%
T_vec = T_vec0;
for i=1:Nt
    T_vec = evolve_core(T_vec, Q_cmb_t(i));
    T_vec_all(i,:) = T_vec;
    T_p = T_vec(1:end-1);
    T_ad = adiabat(T_vec(end),r,R,D);
    E_ad(i) = energy(T_ad,rho_cp_dV);
    E_p(i) = energy(T_p, rho_cp_dV);
    if mod(i,100)==0
        hold on
        plot(r,T_ad+T_p,'.-')
    end
end
%%
t = linspace(0,Nt,Nt);
% plot(t,E_ad)
plot(t,T_vec_all(:,end))
% plot(t,E_p)
% plot(t,E_ad+E_p)
% plot(t,E_ad+E_p+Q_cmb_t*dt)


%%
T_vec = T_vec0;
Q_cmb = 3.5e12;
odef = @(t,y) dTdt_core(y, Q_cmb);
[t, T_v] = ode45(odef, [0, Nt*dt], T_vec0);
T_cmb = T_v(:,end)+T_v(:,end-1); 
plot(t,T_cmb)