function err = minfun(T_cmb, Q_cmb, T, rho_cp_dV, dt)
    T_new = T;
    T_new(T>T_cmb) = T_cmb;
    E_old = energy(T,rho_cp_dV);
    E_new = energy(T_new,rho_cp_dV);
    err = abs(E_old-(E_new+Q_cmb*dt));
end
