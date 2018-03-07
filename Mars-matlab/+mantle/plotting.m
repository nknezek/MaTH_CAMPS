function []=plotting(t,Myr,Qu,n,Ht,A,Ta,flv,Tmat,Rmat,Ra,Crt,f,B,Tb);

%% Plotting

%% Figure
hold on
figure(2); clf
subplot 221;
plot(t/Myr,Ht./A(1)./(Qu(:,1)),'b','linewidth', 2.5);
hold on;
if t(end)/Myr >=4.56e3 
    x=[4.56e3,4.56e3]; y=[0,max(Ht./A(1)./(Qu(:,1)))]; plot(x,y,'k'); % current time
end
if t(end)/Myr >=0.56e3 
    x=[0.56e3,0.56e3]; y=[0,max(Ht./A(1)./(Qu(:,1)))]; plot(x,y,'r'); % Dynamo shutoff
end
xlabel('Time (Myr)'); ylabel('Urey Ratio');

subplot 222;
plot (B, Tb(:,3),'*'); xlabel('Buoyancy number'); ylabel('Temperature[k]');

subplot 223;
plot(Tmat',Rmat'/1000);
xlabel('Temperature (K)');
ylabel('Radius (km)');
ylim([1000,max(max(Ta'))]);

subplot 224;
semilogy(t/Myr,Ra);
%legend(num2str([1:n]'));
xlabel('Time (Myr)');
ylabel('Ra');


%% Figure
figure(3); clf
subplot 221;
plot(t/Myr,Crt','b','linewidth', 2.5);
hold on;
if t(end)/Myr >=4.56e3 
    x=[4.56e3,4.56e3]; y=[0,max(Ht./A(1)./(Qu(:,1)))]; plot(x,y,'k'); % current time
end
if t(end)/Myr >=0.56e3 
    x=[0.56e3,0.56e3]; y=[0,max(Ht./A(1)./(Qu(:,1)))]; plot(x,y,'r'); % Dynamo shutoff
end

xlabel('Time (Myr)'); ylabel('Crustal thickness from melt (km)');

subplot 222;
plot(t/Myr, f(:,1)', 'b', 'linewidth', 2.5);
hold on
plot(t/Myr, f(:,2)', 'g', 'linewidth', 2.5);
plot(t/Myr, f(:,3)', 'r', 'linewidth', 2.5);
if n >3;
    plot(t/Myr, f(:,4)', 'c', 'linewidth', 2.5);
end
xlabel('Time (Myr)');
ylabel('Melt Fraction - layer');
title('Melt Fraction Over Time');
%ylabel('Melt Present - layer');
legend(num2str([1:n]'));
%ylim([0,1]); xlim([min(t/Myr),max(t/Myr)]);
ylim([0,0.25]); 
if t(end)/Myr >=4.56e3; 
    x=[4.56e3,4.56e3]; y=[0,max(max(Ta'))]; plot(x,y,'k'); % current time
end
if t(end)/Myr >=0.56e3; 
    x=[0.56e3,0.56e3]; y=[0,max(max(Ta'))]; plot(x,y,'r'); % Dynamo shutoff
end


% figure
% plot(t/Myr, f(:,1)', 'b', 'linewidth', 2.5);
% hold on
% plot(t/Myr, f(:,2)', 'g', 'linewidth', 2.5);
% plot(t/Myr, f(:,3)', 'r', 'linewidth', 2.5);
% if n >3;
%     plot(t/Myr, f(:,4)', 'c', 'linewidth', 2.5);
% end
% xlabel('Time (Myr)');
% ylabel('Melt Fraction - layer');
% title('Melt Fraction Over Time');
% %ylabel('Melt Present - layer');
% legend(num2str([1:n]'));
% %ylim([0,1]); xlim([min(t/Myr),max(t/Myr)]);
% ylim([0,0.25]); 
% if t(end)/Myr >=4.56e3; 
%     x=[4.56e3,4.56e3]; y=[0,max(max(Ta'))]; plot(x,y,'k'); % current time
% end
% if t(end)/Myr >=0.56e3; 
%     x=[0.56e3,0.56e3]; y=[0,max(max(Ta'))]; plot(x,y,'r'); % Dynamo shutoff
% end
