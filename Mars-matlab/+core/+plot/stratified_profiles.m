function fig = stratified_profiles(pp,pm,pc)
fig = figure();
clf;
Tc = pp.Tc;
time1 = 200;
time2 = 500;
time3 = 1000;
i1 = find(pp.t_plt>time1,1);
i2 = find(pp.t_plt>time2,1);
i3 = find(pp.t_plt>time3,1);
hold on
core.plot.temp_v_r(core.therm.check_stratified(Tc(1,:),pc),pc)
core.plot.temp_v_r(core.therm.check_stratified(Tc(i1,:),pc),pc)
core.plot.temp_v_r(core.therm.check_stratified(Tc(i2,:),pc),pc)
core.plot.temp_v_r(core.therm.check_stratified(Tc(i3,:),pc),pc)
core.plot.temp_v_r(core.therm.check_stratified(Tc(end,:),pc),pc)
legend('init',string(pp.t_plt(i1))+'Myr',string(pp.t_plt(i2))+'Myr',string(pp.t_plt(i3))+'Myr','final')
title('Portion of Core Stratified')
grid on
