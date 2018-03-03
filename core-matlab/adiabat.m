function T = adiabat(T_cmb,r, R, D)
T = T_cmb*exp((R.^2-r.^2)./(D.^2));
end