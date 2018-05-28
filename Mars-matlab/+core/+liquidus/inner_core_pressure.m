function Pic = inner_core_pressure(Tc,pc)
Ph0 = 30;
fliq = core.liquidus.liquidus_polyfit(pc.wtpS);
try
    Pic = fzero(@(P) fliq(P)-interp1(pc.P/1e9,Tc,P),Ph0);
catch
    Pic = nan;
end