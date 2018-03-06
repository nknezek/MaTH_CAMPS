function D = D_scale(rho_cen, pc)
% Adiabatic Lengthscale in Core
% used in T_a = T_cen*exp(-r^2/D^2)

switch nargin
    case 0
        rho_cen = pc.rho_cen;
end
D = sqrt(3*pc.cp/(2*pi*pc.alpha*rho_cen*pc.G));
end