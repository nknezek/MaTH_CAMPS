t_solid = zeros(3,6);
i=1;
% for param_case = [1,]
% for wtpS = [5,]
wtpSvals = [0,5,10,15,20,25];
for param_case = [0,1,2]
j=1;
for wtpS = wtpSvals

folder_casenames = ["nolayer/","hot/","cold/"];
basefolder = './results/';
folder = [char(folder_casenames(param_case+1)),sprintf('%.0fwtpS/', wtpS)];
% mkdir([basefolder,folder])

load([basefolder,folder,'pp.mat'])

pm = pp.pm;
pc = pp.pc;
pc.wtpS = wtpS;
n = pm.n;

tplt = pp.t_plt;

fliq = core.liquidus.liquidus_polyfit(pc.wtpS);
Tliq = fliq(pc.P/1e9);
Tcs = pp.Tvec(:,pm.n+1:end);
%
solid = (sign(Tliq.*ones(1000,150)-Tcs)+1)/2;
ind = min(find(sum(solid,2)));
if isempty(ind)
    t_solid(i,j) = nan;
else
    t_solid(i,j) = tplt(ind);
end
j = j+1;
end
i = i+1;
end

%%
fig = plot(wtpSvals,t_solid','.-');
legend('nolayer','hot','cold')
ylabel('time (Myrs)')
xlabel('wt% S in core')
title('inner-core nucleation timing')
% saveas(fig, 'ic_nuc_time.png')

%%
wtpSvals = [0,5,10,15,20,25];
Tliqs = zeros(6,150);
for i = 1:6
fliq = core.liquidus.liquidus_polyfit(wtpSvals(i));
Tliqs(i,:) = fliq(pc.P/1e9);
end