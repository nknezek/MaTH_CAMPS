function T = adiabat(T_cmb,r)
% calculate adiabat from T_cmb

global pc
switch nargin
    case 1
        r = pc.r;
end
T = T_cmb.*exp((pc.r_cmb.^2-r.^2)./(pc.D.^2));
end