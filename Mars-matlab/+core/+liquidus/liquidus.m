% function pfl, pfh = liquidus_polyfit(wtpS)

fliq = core.liquidus.liquidus_polyfit(pp.pc.wtpS);
%%
Pfit = pp.pc.P/1e9;
Tliq = fliq(Pfit);
plot(Pfit,Tliq)

%%
Tc = pp.Tvec(end,4:end)-600;

plot(Pfit,Tliq)
hold on
plot(Pfit, Tc)

%%
solid = core.liquidus.check_solid_layers(Tc, pp.pc);
plot(Pfit,solid)

%% 
Ph0 = 30;
Pl0 = 21;
Psolid_h = fzero(@(P) polyval(pfliq_h,P)-interp1(Pfit,Tc,P),Ph0);
Psolid_l = fzero(@(P) polyval(pfliq_l,P)-interp1(Pfit,Tc,P,'linear',1e6),Pl0);

%%

