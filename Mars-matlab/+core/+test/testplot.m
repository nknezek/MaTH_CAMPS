x = linspace(0,1,100);
fig = plot(x,x);

basefolder = './results/';
wtpS = 10;
paramcase = 'hot/';
folder = [paramcase,sprintf('%.0fwtpS/', wtpS)];
mkdir([basefolder,folder])
saveas(fig,[basefolder,folder,'testfig.png'])
close all