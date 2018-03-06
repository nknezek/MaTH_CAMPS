function density_v_r(pc)
% Plot density vs radius
%
% parameters
% none

plot(pc.r/1e3,pc.rho)
xlabel('radius (km)')
ylabel('density (kg/m^3)')
title('Density')
end