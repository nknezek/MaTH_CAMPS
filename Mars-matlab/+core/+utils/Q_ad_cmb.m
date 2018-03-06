function Q = Q_ad_cmb(T_cmb, pc)
% Finds conduction down adiabat at CMB
%
% Parameters:
% T [K] - Temp. at CMB


Q = 4*pi*pc.r_cmb^3*pc.k*2*T_cmb/pc.D^2;
end