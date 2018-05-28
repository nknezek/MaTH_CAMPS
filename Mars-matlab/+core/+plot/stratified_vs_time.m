function fig = stratified_vs_time(pp,pm,pc)
fig = figure();
clf;
Tc = pp.Tc;
N = length(Tc);
r_strat = zeros(N,1);
for i=1:N
    r_strat(i) = core.therm.radius_stratified(Tc(i,:),pc);
end
plot(pp.t_plt,(pc.r(end)-r_strat)/1e3)
title('Thickness of Stratified Layer in core over time ')
ylabel('km')
xlabel('time (Myrs)')
grid on
