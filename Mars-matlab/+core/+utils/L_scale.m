function L = L_scale(rho_cen, pc)
% Density Length Scale in Core
% used in rho(r) = rho_cen*exp(-r^2/L^2)
switch nargin
    case 0
        rho_cen = pc.rho_cen;
end       
L = sqrt(3*pc.K_0*(log(rho_cen/pc.rho_0+1)/(2*pi*pc.G*pc.rho_0*rho_cen)));
end