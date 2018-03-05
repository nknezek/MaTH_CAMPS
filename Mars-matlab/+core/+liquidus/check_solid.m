function solid = check_solid(T,r)
rho_cen = 7e3; % kg/m^3
P_cmb = 18e9; % GPa
r_cmb = r(end)+(r(2)-r(1))/2;
P = pressure(r, r_cmb, P_cmb, rho_cen);
solid = zeros(1,length(T));
T_cmb_freeze = 2000;
T_freeze = T_cmb_freeze*(P/P_cmb).^0.5;
solid(T<T_freeze) = 1;
end