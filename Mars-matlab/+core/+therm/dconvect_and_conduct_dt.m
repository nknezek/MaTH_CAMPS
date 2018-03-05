function dTdt = dconvect_and_conduct_dt(T, r, dr, D, A, dt, k, rho_cp_dV, Q_cmb)
T_n = convect_and_conduct(T, r, dr, D, A, dt, k, rho_cp_dV, Q_cmb);
dTdt = (T_n-T)/dt;
end
