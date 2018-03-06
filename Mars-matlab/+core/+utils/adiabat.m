function T = adiabat(T_cmb,pc, r)
% calculate adiabat from T_cmb
%
% params
% T_cmb [K] - CMB temp
% pc [struct] - core parameters
% r [m] - radial locations
switch nargin
    case 2
        r = pc.r;
end
T = T_cmb.*exp((pc.r_cmb.^2-r.^2)./(pc.D.^2));
end