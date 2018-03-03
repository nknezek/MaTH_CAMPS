function T_out = conduct(T, dr, dt, k, A, rho_cp_dV, Q_cmb)

qe = -k*(T(2:end)-T(1:end-1))./dr; % W/m^2
Qe = zeros(1, length(T)); % W
Qe(1:end-1) = qe.*A(2:end-1); % W
Qe(end) = Q_cmb; % W
Q = zeros(1,length(T)); % W
Q(2:end) = Qe(2:end)-Qe(1:end-1); % W
dT = -Q.*dt ./rho_cp_dV; % K
T_out = T+dT; % K
end