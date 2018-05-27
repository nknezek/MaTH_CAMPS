function fig = melt(pp,pm,pc)
% plots the melt over time

Myr = pm.Myr;
t = pp.tvec;
n = pm.n;
Tm = pp.Tm;
fig = figure(); 
clf;
subplot 121;
plot(pp.t_plt,pp.Crt','b','linewidth', 2.5);
hold on;
if t(end)/Myr >=4.56e3 
    x=[4.56e3,4.56e3]; 
    y=[0,max(pp.Ht./pm.A(1)./(pp.Qu(:,1)))]; 
    plot(x,y,'k'); % current time
end
if t(end)/Myr >=0.56e3 
    x=[0.56e3,0.56e3]; 
    y=[0,max(pp.Ht./pm.A(1)./(pp.Qu(:,1)))]; 
    plot(x,y,'r'); % Dynamo shutoff
end

xlabel('Time (Myr)'); 
ylabel('Crustal thickness from melt (km)');


subplot 122;
meltf_per_90myr = zeros(1000/50,pp.pm.n);
for i=1:pp.pm.n
    meltf_per_90myr(:,i) = sum(reshape([diff(pp.flv(:,1));0],50,[]))./pp.pm.M(i);
tplt = pp.t_plt(1:50:end);
plot(tplt,meltf_per_90myr(:,1), 'b', 'linewidth', 2.5);
hold on
plot(tplt,meltf_per_90myr(:,2), 'g', 'linewidth', 2.5);
if n>2
    plot(tplt,meltf_per_90myr(:,3), 'r', 'linewidth', 2.5);
elseif n >3
    plot(pp.t_plt, pp.f(:,4)', 'c', 'linewidth', 2.5);
end
xlabel('Time (Myr)');
ylabel('Mass Fraction of Layer Melted per 90 Myr ');
title('Mass Fraction of Layer Melted per 90 Myr');
%ylabel('Melt Present - layer');
legend(num2str([1:n]'));

%ylim([0,1]); xlim([min(pp.t_plt),max(pp.t_plt)]);
% ylim([0,0.25]); 

% if t(end)/Myr >=4.56e3 
%     x=[4.56e3,4.56e3]; y=[0,max(max(Tm'))]; plot(x,y,'k'); % current time
% end
% if t(end)/Myr >=0.56e3
%     x=[0.56e3,0.56e3]; y=[0,max(max(Tm'))]; plot(x,y,'r'); % Dynamo shutoff
% end
end
