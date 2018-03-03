function E = energy(T,rho_cp_dV)
E = sum(T.*rho_cp_dV);
end