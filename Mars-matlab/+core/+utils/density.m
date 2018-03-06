function rho = density(r, rho_cen, D, pc)
% function to compute density 
%
% parameters
% r [m] (1,N) - radial centers of layers
% rho_cen [kg/m^3] - density at center of core
% D [m] - adiabatic lengthscale

switch nargin
    case 0
        r = pc.r;
        rho_cen = pc.rho_cen;
        D = pc.D;
    case 1
        rho_cen = pc.rho_cen;
        D = pc.D;
    case 2
        D = pc.D;        
end
  
rho = rho_cen*exp(-r.^2/D^2);
end