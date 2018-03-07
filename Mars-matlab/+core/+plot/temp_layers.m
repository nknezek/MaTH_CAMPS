function fig = temp_layers(pp,pm,pc)
Tm = pp.Tm;
Tc = pp.Tc;
fig = figure();
clf
hold on
title('Temp evolution')
plot(pp.t_plt,Tm(:,1))
plot(pp.t_plt,Tm(:,2))
plot(pp.t_plt,Tm(:,3))
plot(pp.t_plt,Tc(:,end))
plot(pp.t_plt,Tc(:,1))
xlabel('time')
ylabel('Temperature (K)')
legend('lithosphere','mantle','lower mantle layer','cmb','core center')
grid on
hold off
end