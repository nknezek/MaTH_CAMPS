function fig = melt(pp,pm,pc)
% plots the melt over time

Myr = pm.Myr;
t = pp.tvec;
n = pm.n;
Tm = pp.Tm;
fig = figure(); 
clf;
subplot 221;
plot(pp.t_plt,pp.Crt','b','linewidth', 2.5);
hold on;
if t(end)/Myr >=4.56e3 
    x=[4.56e3,4.56e3]; y=[0,max(pp.Ht./pm.A(1)./(pp.Qu(:,1)))]; plot(x,y,'k'); % current time
end
if t(end)/Myr >=0.56e3 
    x=[0.56e3,0.56e3]; y=[0,max(pp.Ht./pm.A(1)./(pp.Qu(:,1)))]; plot(x,y,'r'); % Dynamo shutoff
end

xlabel('Time (Myr)'); ylabel('Crustal thickness from melt (km)');


subplot 222;
plot(pp.t_plt, pp.f(:,1)', 'b', 'linewidth', 2.5);
hold on
plot(pp.t_plt, pp.f(:,2)', 'g', 'linewidth', 2.5);
plot(pp.t_plt, pp.f(:,3)', 'r', 'linewidth', 2.5);
if n >3
    plot(pp.t_plt, pp.f(:,4)', 'c', 'linewidth', 2.5);
end
xlabel('Time (Myr)');
ylabel('Melt Fraction - layer');
title('Melt Fraction Over Time');
%ylabel('Melt Present - layer');
legend(num2str([1:n]'));
%ylim([0,1]); xlim([min(pp.t_plt),max(pp.t_plt)]);
ylim([0,0.25]); 
if t(end)/Myr >=4.56e3 
    x=[4.56e3,4.56e3]; y=[0,max(max(Tm'))]; plot(x,y,'k'); % current time
end
if t(end)/Myr >=0.56e3; 
    x=[0.56e3,0.56e3]; y=[0,max(max(Tm'))]; plot(x,y,'r'); % Dynamo shutoff
end
end
