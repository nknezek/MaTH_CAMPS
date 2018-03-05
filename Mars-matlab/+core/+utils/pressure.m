function P = pressure(r)
% compute pressure at radial locaitons r in core
%
% parameter
% r [m] - radial location 

global pc
switch nargin
    case 0
        r = pc.r;
end
P = pc.P_cmb+ 4*pi*pc.G*pc.rho_cen^2/3.*((3*pc.r_cmb.^2/10 - pc.L^2/5).*exp(-pc.r_cmb.^2/pc.L^2) - (3*r.^2/10-pc.L^2/5).*exp(-r.^2/pc.L^2));
end