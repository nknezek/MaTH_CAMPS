function fig = Urey_Buoyancy(pp,pm,pc)
Myr = pm.Myr;

hold on
fig = figure();
clf;
subplot 221;
plot(pp.t_plt,pp.Ht./pm.A(1)./(pp.Qu(:,1)),'b','linewidth', 2.5);
hold on;
if pp.tvec(end)/Myr >=4.56e3 
    x=[4.56e3,4.56e3]; y=[0,max(pp.Ht./pm.A(1)./(pp.Qu(:,1)))]; plot(x,y,'k'); % current time
end
if pp.tvec(end)/Myr >=0.56e3 
    x=[0.56e3,0.56e3]; y=[0,max(pp.Ht./pm.A(1)./(pp.Qu(:,1)))]; plot(x,y,'r'); % Dynamo shutoff
end
xlabel('Time (Myr)'); ylabel('Urey Ratio');

subplot 222;
plot(pp.B, pp.Tb(:,3),'*'); xlabel('Buoyancy number'); ylabel('Temperature[k]');

subplot 223;
plot(pp.Tmat',pp.Rmat'/1000);
xlabel('Temperature (K)');
ylabel('Radius (km)');
ylim([1000,max(max(pp.Tm'))]);

subplot 224;
semilogy(pp.t_plt,pp.Ra);
%legend(num2str([1:n]'));
xlabel('Time (Myr)');
ylabel('Ra');
end
