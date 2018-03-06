function E = energy(T, pc)
% compute total internal energy in core
% T [K] - (1,N) vector of temperatures in each radial layer

E = sum(T.*pc.rho_cp_dV);
end