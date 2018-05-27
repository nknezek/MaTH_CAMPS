function fig = Urey_Buoyancy(pp,pm,pc)
Myr = pm.Myr;

fig = figure();
clf
hold on

subplot 221;
plot(pp.t_plt,pp.Ht./pm.A(1)./(pp.Qu(:,1)),'b','linewidth', 2.5);
hold on;
if pp.tvec(end)/Myr >=4.56e3 
    x=[4.56e3,4.56e3]; 
    y=[0,max(pp.Ht./pm.A(1)./(pp.Qu(:,1)))]; 
    plot(x,y,'k'); % current time
end
if pp.tvec(end)/Myr >=0.56e3 
    x=[0.56e3,0.56e3]; 
    y=[0,max(pp.Ht./pm.A(1)./(pp.Qu(:,1)))]; 
    plot(x,y,'r'); % Dynamo shutoff
end
xlabel('Time (Myr)'); 
ylabel('Urey Ratio');

subplot 222;
plot(pp.B, pp.Tb(:,3),'*'); 
xlabel('Buoyancy number'); 
ylabel('Temperature[k]');


subplot 223;
Tplt = pp.Tvec(1:100:end,:);
Tplt(:,4:end) = pp.Tvec(1:100:end,end:-1:4);
Tplt = Tplt';
Rplt = [pm.R(1:3)',pc.r(end:-1:1)]/1000;
plot(Tplt,Rplt);
xlabel('Temperature (K)');
ylabel('Radius (km)');
ylim([0,3500]);
xlim([500,3000]);

subplot 224;
semilogy(pp.t_plt,pp.Ra);
%legend(num2str([1:n]'));
xlabel('Time (Myr)');
ylabel('Ra');
end
