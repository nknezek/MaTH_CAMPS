function dTdt = dconvect_and_conduct_dt(T, dt, Q_cmb, pc, debug, debugerr)
% wrapper function to compute derivatives based on the Euler core
% timestepping
% 
% parameters
% T [K] (1,N) - temperatures in each layer
% dt [s] - timestep size to compute dTdt
% Q_cmb [W] - heat flow out of core
% debug - flag to plot at each timestep
% debugerr - flag to plot when error occurs

switch nargin
    case 4
        debug = 0;
        debugerr = 0;
end

T_n = core.therm.convect_and_conduct(T, dt, Q_cmb, pc, debug, debugerr);
dTdt = (T_n-T)/dt;
end
