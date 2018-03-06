function grav_pot = grav_pot(r, pc)

psi = 2/3*pi*pc.G*pc.rho_cen*r.^2.*(1-3*r.^2/(10*pc.L.^2));
psi_cmb = 2/3*pi*pc.G*pc.rho_cen*pc.r_cmb.^2.*(1-3*pc.r_cmb.^2/(10*pc.L.^2));
grav_pot = psi-psi_cmb;
end