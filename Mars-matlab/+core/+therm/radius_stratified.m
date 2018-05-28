function r_strat = radius_stratified(T, pc)
stratified = core.therm.check_stratified(T, pc);
sind = find(stratified);
if isempty(sind)
    r_strat = pc.r(end);
else
    r_strat = pc.r(min(sind));
end