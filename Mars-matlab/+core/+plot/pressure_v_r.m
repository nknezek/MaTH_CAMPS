function pressure_v_r()
% Plot pressure vs radius
%
% parameters
% none

global pc

plot(pc.r/1e3,pc.P/1e9)
xlabel('radius (km)')
ylabel('pressure (GPa)')
title('Pressure')
end