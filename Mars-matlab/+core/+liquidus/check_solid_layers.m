function solid = check_solid_layers(T,pc)
% check if a particular layer in the core is below or above the liquidus
% value

fliq = core.liquidus.liquidus_polyfit(pc.wtpS);
Tliq = fliq(pc.P/1e9);
solid = (sign(Tliq-T)+1)/2;
end
