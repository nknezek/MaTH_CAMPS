%% Test Thermal Evolution Code 
% Adiabatic convection up and thermal conduction down
% Nicholas Knezek

function [T_ad, T_p] = evolve_core(T_ad, T_p, Q_cmb)
N = length(T_ad);
R = 2000e3; % [m]
dr = R/N; % [m]

dt = 365.25*24*3600*1e6; % s (1 Myr)

D = 6203e3; % m, adiabatic lengthscale
cp = 840; %J/kg-K = m^2/s^2-K
rho = linspace(12e4, 9e4, N); % kg/m^3
k = 135; % J/m-K-s = kg-m/K-s^3

re = linspace(0,R,N+1); % m , edge radial locations of layers (0,dr,...,R)
r = linspace(dr/2, R-dr/2, N); % m , center radial locations of layers
dV = 4*pi*dr*r.^2; % m^3 , volume of each layer
A = 4*pi*re.^2; % m^2 , surface area of each radial layer edge

rho_cp_dV = rho.*cp.*dV; % J/K , Joules per Kelvin of each layer

E_p = energy(T_p,rho_cp_dV); % J, Energy stored in perturbation off adiabat
Q_ad_c = Q_ad_cmb(T_ad(end), k, R, D); % W

if Q_cmb >= Q_ad_c
    if Q_cmb*dt >= E_p
        T_p = linspace(0,0,N);
        Q_cmb_extra = Q_cmb-E_p/dt;
        minfun_ad_c  = @(T_c) minfun_ad(T_c, Q_cmb_extra, T_ad, rho_cp_dV, dt, r, R, D);
        T_cmb_new = fminbnd(minfun_ad_c,0,T_ad(end));
        T_ad = adiabat(T_cmb_new, r, R, D);
    else 
        minfun_p  = @(T_c) minfun(T_c, Q_cmb, T_p, rho_cp_dV, dt);
        Tpcmb_new = fminbnd(minfun_p, 0, T_p(end));
        T_p(T_p>Tpcmb_new) = Tpcmb_new;
    end
else
    Q_p_cmb = Q_cmb-Q_ad_c;
    T_p = conduct(T_p, dr, dt, k, A, rho_cp_dV, Q_p_cmb);
end

end
