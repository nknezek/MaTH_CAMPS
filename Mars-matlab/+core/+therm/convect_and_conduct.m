function T_out = convect_and_conduct(T, dt, Q_cmb, pc, debug, debugerr)
% Computes one timestep convecting and conducting temperature in core
% 
% Parameters
% T [K] - (1,N) vector of temperatures in each layer
% dt [s] - timestep
% Q_cmb [W = J/s] - heat out of top of core

switch nargin
    case 4
        debug = 0;
        debugerr = 0;
end
r = pc.r;
N = pc.N;
E_i = core.utils.energy(T,pc); % Initial Energy of Core

% Conduct heat down gradients
qe = -pc.k*(T(2:end)-T(1:end-1))./pc.dr; % W/m^2, conductive heat flow across layer boundaries
Qe = zeros(1, length(T)); % W, conductive heat flow across layer boundaries 
Qe(1:end-1) = qe.*pc.A(2:end-1); % W
Qe(end) = Q_cmb; % W , heat flow out top of core 
Q = zeros(1,length(T)); % W, net conductive heat flux out of each layer
Q(2:end) = Qe(2:end)-Qe(1:end-1); % W
Q(1) = Qe(1); % W, heat flow from central point only goes up
dT_cond = -Q.*dt ./pc.rho_cp_dV; % K
T_cond = T+dT_cond; % K

if max(dT_cond)>1
    if debugerr==1 || debug==1
        core.plot.temp_v_r(dT_cond, pc,'dTemperature')
        figure()
        hold on
        core.plot.temp_v_r(T, pc, 'init')
        core.plot.temp_v_r(T_cond,pc, 'cond')
        legend('init','cond')
    end
%     error('conduction dT = %.2f > 1 K, use a smaller timestep. dt=%.2f Myr currently. error',max(dT_cond), dt/pc.Myr)
end

% Find regions in which convection occurs, then convect them
Ch_ad_total = 0; % total heat capacity of current adiabatic region
T_out = zeros(1,N); % final temperatures after timestep
in_an_adiabat = 0; % is the current layer in a convective area
T_regions = T_cond; % temperatures to use to define convective/conductive regions
Q_typical = 1e11;
epsilonTp1 = 1+min(Q_typical*dt/mean(pc.rho_cp_dV)*1e-4,1e-4);
for i=N:-1:2
    T_im1a = core.utils.adiabat_from(T_regions(i),r(i),r(i-1), pc);
    if T_im1a <= T_regions(i-1)*epsilonTp1 % in convective regime
        if in_an_adiabat == 0 % in top layer of convective regime
            i_ad_top = i;
            Q_ad_top = Qe(i);
        end
        in_an_adiabat = 1;
        Ch_ad_total = Ch_ad_total + pc.rho_cp_dV(i)*exp((r(i_ad_top)^2 - r(i)^2)/pc.D^2);
    elseif T_im1a > T_regions(i-1)*epsilonTp1 % in conductive regime
        if in_an_adiabat == 1% in top layer of conductive regime
            % deal with adiabat that we just left
            Ch_ad_total = Ch_ad_total + pc.rho_cp_dV(i)*exp((r(i_ad_top)^2 - r(i)^2)/pc.D^2);
            T_ad_top_new = T(i_ad_top)-Q_ad_top*dt/Ch_ad_total;
            for j=i:i_ad_top
                T_out(j) = core.utils.adiabat_from(T_ad_top_new, r(i_ad_top), r(j), pc);
            end
            in_an_adiabat = 0;
            Ch_ad_total = 0;
        else % in at least second layer of conductive region
            T_out(i) = T_cond(i);
        end
    end
end
if in_an_adiabat == 1
    i = 1;
    Ch_ad_total = Ch_ad_total + pc.rho_cp_dV(i)*exp((r(i_ad_top)^2 - r(i)^2)/pc.D^2);
    T_ad_top_new = T(i_ad_top)-Q_ad_top*dt/Ch_ad_total;
    for j=i:i_ad_top
        T_out(j) = core.utils.adiabat_from(T_ad_top_new, r(i_ad_top), r(j), pc);
    end
end
if debug==1
    h = figure();
    hold on
    core.plot.temp_v_r(T, 'init')
    core.plot.temp_v_r(T_cond, 'cond')
    core.plot.temp_v_r(T_out, 'final')
    legend('init','cond','final')
    waitfor(h)
    hold off
end

% Solidifcation energy balance
% s0 = check_solid(T,r);
% s1 = check_solid(T_out,r);
% alpha_c = 1;
% r_cmb = r(end)+dr/2;
% rho_cen = 7e3;
% psi_r = grav_pot(r, r_cmb, rho_cen);
% dMsolid = sum((s1-s0).*rho.*dV);
% dEg = dMsolid*sum(rho.*dV.*psi_r)*alpha_c;
% L_h = 750e3; % J/kg
% dEl = dMsolid*L_h;

% Correction to make sure you satisfy energy conservation
E_o = core.utils.energy(T_out, pc); % final energy of core, J
dEq = Qe(end)*dt; % Heat pulled out of core, J
dEt = E_i-E_o; % change in internal energy, J
% dE_corr = dEt-dEq+dEg+dEl;
dE_corr = dEt-dEq;
dT_corr = dE_corr/sum(pc.rho_cp_dV.*exp(-r.^2./pc.D^2));
if max(abs(dT_corr)) > 1
%     error('dT_corr = %.2f K is larger than 1 K. Maybe use a smaller timestep?',dT_corr)
end
T_out = T_out + dT_corr.*exp(-r.^2./pc.D^2);
% *exp(-r.^2./D^2)
end