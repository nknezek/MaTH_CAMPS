function T = adiabat_from(T_1,r_1,r_2, pc)
% calculate adiabat from one radius to another
%
% parameters
% T_1 [K] - temperature to calculate from
% r_1 [m] - radius to calculate from
% r_2 [m] - radius to calculate new temperature
% pc

T = T_1.*exp((r_1.^2-r_2.^2)./(pc.D.^2));
end