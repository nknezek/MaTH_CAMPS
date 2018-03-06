function T_out = conduct(T, dt, Q_cmb, pc)
% conduct energy between radial layers in core

qe = -pc.k*(T(2:end)-T(1:end-1))./pc.dr; % W/m^2
Qe = zeros(1, length(T)); % W
Qe(1:end-1) = qe.*pc.A(2:end-1); % W
Qe(end) = Q_cmb; % W
Q = zeros(1,length(T)); % W
Q(2:end) = Qe(2:end)-Qe(1:end-1); % W
Q(1) = Qe(1);
dT = -Q.*dt ./pc.rho_cp_dV; % K
T_out = T+dT; % K
end