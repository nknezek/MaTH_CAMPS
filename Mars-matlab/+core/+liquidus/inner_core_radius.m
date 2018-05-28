function Ric = inner_core_radius(Tc,pc)
% returns inner-core radius in (m)
Pic = core.liquidus.inner_core_pressure(Tc,pc);
if isnan(Pic)
    Ric = 0;
elseif Pic <= core.utils.pressure(pc.r(end),pc)/1e9
    Ric = pc.r(end);
else
    Ric = fzero(@(r) core.utils.pressure(r,pc)/1e9-Pic,1000e3);
end