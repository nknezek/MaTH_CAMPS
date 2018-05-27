function fig = run_summary(pp,pm,pc)
Myr = pm.Myr;
A = pm.A;
n = pm.n;

fig = figure(); 
clf

hold on
subplot 221;
plot(pp.t_plt,pp.Qu'*1e3, 'linewidth', 2.5)
xlabel('Time (Myr)');
ylabel('Heat flux (mW/m^2)');
ylim([-20, 80]);
title('Top of each layer');
legend(num2str([1:n]'));
hold on;
plot(pp.t_plt,pp.Ht./pm.A(1)*1e3,'k','linewidth',3.5);
if pp.t_plt(end) >=4.56e3
    x=[4.56e3,4.56e3]; y=[-50,150]; plot(x,y,'k'); % current time
end
if pp.t_plt(end) >=0.56e3
    x=[0.56e3,0.56e3]; y=[-50,150]; plot(x,y,'r'); % Dynamo shutoff
end

subplot 222;
plot(pp.t_plt,pp.Tm', 'linewidth', 2.5);
xlabel('time (Myr)');
ylabel('temperature (K)');
title('Average for each layer (with melt)');
legend(num2str([1:n]'));
hold on;
if pp.t_plt(end) >=4.56e3; 
    x=[4.56e3,4.56e3]; y=[0,max(max(pp.Tm'))]; plot(x,y,'k'); % current time
end
if pp.t_plt(end) >=0.56e3; 
    x=[0.56e3,0.56e3]; y=[0,max(max(pp.Tm'))]; plot(x,y,'r'); % Dynamo shutoff
end

% subplot 223;
% plot(pp.t_plt, f(:,1)', 'b', 'linewidth', 2.5);
% hold on
% plot(pp.t_plt, f(:,2)', 'g', 'linewidth', 2.5);
% plot(pp.t_plt, f(:,3)', 'r', 'linewidth', 2.5);
% if n >3;
%     plot(pp.t_plt, f(:,4)', 'c', 'linewidth', 2.5);
% end
% xlabel('Time (Myr)');
% ylabel('Melt Fraction - layer');
% title('Melt Fraction Over Time');
% %ylabel('Melt Present - layer');
% legend(num2str([1:n]'));
% %ylim([0,1]); xlim([min(pp.t_plt),max(pp.t_plt)]);
% ylim([0,(max(max(f)))]); 
% if pp.t_plt(end) >=4.56e3; 
%     x=[4.56e3,4.56e3]; y=[0,max(max(Ta'))]; plot(x,y,'k'); % current time
% end
% if pp.t_plt(end) >=0.56e3; 
%     x=[0.56e3,0.56e3]; y=[0,max(max(Ta'))]; plot(x,y,'r'); % Dynamo shutoff
% end

subplot 223;
flv_dt = diff(pp.flv,1,1)./diff(pp.t_plt);

semilogy(pp.t_plt(1:end-1), flv_dt', 'b', 'linewidth', 2.5);
hold on
xlabel('Time (Myr)');
ylabel('Melt production (kg/Myr)');
%ylim([0,1]); xlim([min(pp.t_plt),max(pp.t_plt)]);
ylim([0,(max(max(flv_dt)))]); 
if pp.t_plt(end) >=4.56e3
    x=[4.56e3,4.56e3]; y=[0,max(max(flv_dt'))]; plot(x,y,'k'); % current time
end
if pp.t_plt(end) >=0.56e3
    x=[0.56e3,0.56e3]; y=[0,max(max(flv_dt'))]; plot(x,y,'r'); % Dynamo shutoff
end



subplot 224;
plot(pp.t_plt,pp.Tmat', 'linewidth', 2.5)
xlabel('time (Myr)');
ylabel('temperature (K)');
title('Temperatures (without melt) (actually with melt with euler)');
hold on;
if pp.t_plt(end) >=4.56e3
    x=[4.56e3,4.56e3];y=[0,max(max(pp.Tm'))]; plot(x,y,'k'); % current time
end
if pp.t_plt(end) >=0.56e3
    x=[0.56e3,0.56e3];y=[0,max(max(pp.Tm'))]; plot(x,y,'r'); % Dynamo shutoff
end
end