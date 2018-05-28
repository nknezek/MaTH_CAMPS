function fig = temp_profiles(pp,pm,pc)
fig = figure();
clf;
fliq = core.liquidus.liquidus_polyfit(pc.wtpS);
Tliq = fliq(pc.P/1e9);
Tc = pp.Tc;
time1 = 200;
time2 = 500;
time3 = 1000;
i1 = find(pp.t_plt>time1,1);
i2 = find(pp.t_plt>time2,1);
i3 = find(pp.t_plt>time3,1);
hold on

core.plot.temp_v_r(Tc(1,:),pc)
core.plot.temp_v_r(Tc(i1,:),pc)
core.plot.temp_v_r(Tc(i2,:),pc)
core.plot.temp_v_r(Tc(i3,:),pc)
core.plot.temp_v_r(Tc(end,:),pc)
plot(pc.r/1e3,Tliq,'k--')
legend('init',string(pp.t_plt(i1))+'Myr',string(pp.t_plt(i2))+'Myr',string(pp.t_plt(i3))+'Myr','final','liquidus')
title('Core Temperature')
grid on
