function err = minfun_ad(T_c, Q_cmb, T_ad, rho_cp_dV, dt, r, R, D)
    T_ad_new = adiabat(T_c, r, R, D);
    E_old = energy(T_ad,rho_cp_dV);
    E_new = energy(T_ad_new,rho_cp_dV);
    err = abs(E_old-(E_new+Q_cmb*dt));
end
