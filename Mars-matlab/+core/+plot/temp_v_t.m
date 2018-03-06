function temp_v_t(t,T,pc,Nplt)
% Plot temperature vs radius
%
% parameters
% T [K] - (1,N) temperature at each radial layer 
% t [s] - 

switch nargin
    case 3
        Nplt = 200;
end
di = round(length(t)/Nplt);
plot(t(1:di:end)/pc.Myr,T(1:di:end))
xlabel('time (Myr)')
ylabel('temperature (K)')
title('temperature')
end