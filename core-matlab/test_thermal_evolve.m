%% Test Thermal Evolution Code 
% Adiabatic convection up and thermal conduction down
% Nicholas Knezek

N = 100;
R = 2000e3; % [m]
dr = R/N; % [m]
dt = 365.25*24*3600*1e6; % s (1 Myr)
D = 6203e3; % m, adiabatic lengthscale

Nt = 1000;

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
Q_cmb_t(1:200) = 2e12;
Q_cmb_t(200:500) = -1e11;
Q_cmb_t(500:end) = 1e13;
%%
Q_ad_cmb = Q_ad_cmb(T_ad(end), k, R, D); % W

for i=1:Nt
    Q_cmb = Q_cmb_t(i);
    E_p = energy(T_p,rho_cp_dV);
    if Q_cmb >= Q_ad_cmb
        if Q_cmb*dt >= E_p
            disp("1")
            T_p = linspace(0,0,N);
            Q_cmb_extra = Q_cmb-E_p/dt;
            minfun_ad_c  = @(T_c) minfun_ad(T_c, Q_cmb_extra, T_ad, rho_cp_dV, dt, r, R, D);
            T_cmb_new = fminbnd(minfun_ad_c,0,T_ad(end));
            T_ad = adiabat(T_cmb_new, r, R, D);
        else 
            disp("2")
            minfun_p  = @(T_c) minfun(T_c, Q_cmb, T_p, rho_cp_dV, dt);
            Tpcmb_new = fminbnd(minfun_p, 0, T_p(end));
            T_p(T_p>Tpcmb_new) = Tpcmb_new;
        end
    else
        disp("3")
        Q_p_cmb = Q_cmb-Q_ad_cmb;
        T_p = conduct(T_p, dr, dt, k, A, rho_cp_dV, Q_p_cmb);
    end
    if mod(i,100)==0
        hold on
        plot(r,T_ad+T_p)
    end
end

%%
plot(r,T_0)
hold on
% plot(r,T_1)
plot(r,T_ad+T_p)