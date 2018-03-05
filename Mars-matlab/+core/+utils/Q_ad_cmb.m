function Q = Q_ad_cmb(T_cmb)
% Finds conduction down adiabat at CMB
%
% Parameters:
% T [K] - Temp. at CMB

global pc

Q = 4*pi*pc.r_cmb^3*pc.k*2*T_cmb/pc.D^2;
end