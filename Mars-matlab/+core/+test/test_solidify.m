%%
P = pp.pc.P;
T = pp.Tvec(end,4:end);
plot(P,T)

pf = polyfit(P,T);

