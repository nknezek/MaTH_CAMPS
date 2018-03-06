function model_params(T_cmb, pc)
% Plot density, pressure, adiabat
%
% parameters
% T_cmb [K] - Temperature at CMB

switch nargin
    case 0
        T_cmb = 2400; %[K]
end
T_ad = core.utils.adiabat(T_cmb);
subplot(221)
core.plot.density_v_r()
grid on
subplot(222)
core.plot.pressure_v_r()
grid on
subplot(223)
core.plot.temp_v_r(T_ad)
grid on
end