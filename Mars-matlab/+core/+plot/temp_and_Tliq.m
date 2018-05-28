function fig = temp_and_Tliq(pp,pm,pc)
fig = figure();
clf;
Tc = pp.Tc;
hold on
core.plot.temp_v_r(Tc(end,:),pc)
fliq = core.liquidus.liquidus_polyfit(pc.wtpS);
Tliq = fliq(pc.P/1e9);
core.plot.temp_v_r(Tliq,pc)
legend('core temperature at present','liquidus')
title('Core Temperature vs Liquidus')
grid on
