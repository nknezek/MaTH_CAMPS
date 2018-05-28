% for wtpS = [5,15,25]
% for param_case = 1:2
% for param_case = [0,]

wtpS = 5;
param_case = 0;

folder_casenames = ["nolayer/","hot/","cold/"];
basefolder = './results/';
folder = [char(folder_casenames(param_case+1)),sprintf('%.0fwtpS/', wtpS)];
mkdir([basefolder,folder])

load([basefolder,folder,'pp.mat'])

pm = pp.pm;
pc = pp.pc;
pc.wtpS = wtpS;
n = pm.n;


% fig = mantle.plot.run_summary(pp,pm,pc);
% saveas(fig,[basefolder,folder,'run_summary.png'])
% 
% fig = mantle.plot.Urey_Buoyancy(pp,pm,pc);
% saveas(fig,[basefolder,folder,'Urey.png'])
% 
% fig = mantle.plot.melt(pp,pm,pc);
% saveas(fig,[basefolder,folder,'melt.png'])
% 
% fig = core.plot.temp_layers(pp,pm,pc);
% saveas(fig,[basefolder,folder,'temp_layers.png'])
% 
fig = core.plot.temp_profiles(pp,pm,pc);
saveas(fig,[basefolder,folder,'temp_profiles.png'])

% fig = core.plot.temp_and_Tliq(pp,pm,pc);
% saveas(fig,[basefolder,folder,'temp_vs_liq.png'])

% fig = core.plot.stratified_vs_time(pp,pm,pc);
% saveas(fig,[basefolder,folder,'stratified_vs_time.png'])

% close all 
% end
% end