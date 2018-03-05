%% Test Thermal Evolution Code 
% Adiabatic convection up and thermal conduction down
% Nicholas Knezek

function dTdt_vec = dTdt_core(T_vec, Q_cmb)
T_vec = T_vec';
N = length(T_vec)-1;
T_adcmb = T_vec(end);
R = 2000e3; % [m]
dr = R/N; % [m]

dt = 365.25*24*3600*1e6; % s (10 kyr)

D = 6203e3; % m, adiabatic lengthscale
cp = 840; %J/kg-K = m^2/s^2-K
rho = linspace(12e4, 9e4, N); % kg/m^3
k = 135; % J/m-K-s = kg-m/K-s^3

re = linspace(0,R,N+1); % m , edge radial locations of layers (0,dr,...,R)
r = linspace(dr/2, R-dr/2, N); % m , center radial locations of layers
dV = 4*pi*dr*r.^2; % m^3 , volume of each layer
A = 4*pi*re.^2; % m^2 , surface area of each radial layer edge

rho_cp_dV = rho.*cp.*dV; % J/K , Joules per Kelvin of each layer

T_p = T_vec(1:end-1);
T_ad = adiabat(T_adcmb,r,R,D);

E_p = energy(T_p,rho_cp_dV); % J, Energy stored in perturbation off adiabat
Q_ad_c = Q_ad_cmb(T_adcmb, k, R, D); % W

Q_cmb_extra = Q_cmb-Q_ad_c;
minfun_ad_c  = @(T_c) minfun_ad(T_c, Q_cmb_extra, T_ad, rho_cp_dV, dt, r, R, D);
T_adcmb_new = fminbnd(minfun_ad_c,0,T_ad(end));

if Q_cmb_extra >= 0.
    
end
if Q_cmb >= Q_ad_c
    if Q_cmb*dt >= E_p
%         disp('1')
        T_p_new = linspace(0,0,N);
        Q_cmb_extra = Q_cmb-E_p/dt;
        minfun_ad_c  = @(T_c) minfun_ad(T_c, Q_cmb_extra, T_ad, rho_cp_dV, dt, r, R, D);
        T_cmb_new = fminbnd(minfun_ad_c,0,T_ad(end));
    else 
%         disp('2')
        minfun_p  = @(T_c) minfun(T_c, Q_cmb, T_p, rho_cp_dV, dt);
        Tpcmb_new = fminbnd(minfun_p, 0, T_p(end));
        T_p_new = T_p;
        T_p_new(T_p>Tpcmb_new) = Tpcmb_new;
        minfun_ad_c  = @(T_c) minfun_ad(T_c, Q_cmb_extra, T_ad, rho_cp_dV, dt, r, R, D);
        T_cmb_new = fminbnd(minfun_ad_c,0,T_ad(end));
    end
else
%     disp('3')
    Q_p_cmb = Q_cmb-Q_ad_c;
    T_p_new = conduct(T_p, dr, dt, k, A, rho_cp_dV, Q_p_cmb);
    T_cmb_new = T_adcmb;
end
dTdt_vec = (horzcat(T_p_new,T_cmb_new)-T_vec)'/dt;
end
