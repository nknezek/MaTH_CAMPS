%% Plot adiabat vs P

core.parameters
T_cmb = 2300;
P = core.utils.pressure()/1e9;
hold on
for Tc=1400:100:2400
    Ta = core.utils.adiabat(Tc);
    plot(P,Ta)
end
ylim([1300,2700])
xlim([20,40])
grid on
title('adiabat in Mars, 2000km core')
xlabel('P (GPa)')
ylabel('T (K)')

%%
dP = P(end)-P(1);
dT = Ta(end)-Ta(1);
dT/dP