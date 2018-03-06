function temp_v_r(T, pc, title_text)
% Plot temperature vs radius
%
% parameters
% T [K] - (1,N) temperature at each radial layer 

switch nargin
    case 1
        title_text = 'Temperature';
end

plot(pc.r/1e3,T)
xlabel('radius (km)')
ylabel('temperature (K)')
title(title_text)
end