function E = energy(T)
% compute total internal energy in core
% T [K] - (1,N) vector of temperatures in each radial layer

global pc
E = sum(T.*pc.rho_cp_dV);
end