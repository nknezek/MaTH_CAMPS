function density_v_r()
% Plot density vs radius
%
% parameters
% none

global pc

plot(pc.r/1e3,pc.rho)
xlabel('radius (km)')
ylabel('density (kg/m^3)')
title('Density')
end