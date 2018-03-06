function pressure_v_r(pc)
% Plot pressure vs radius
%
% parameters
% none


plot(pc.r/1e3,pc.P/1e9)
xlabel('radius (km)')
ylabel('pressure (GPa)')
title('Pressure')
end