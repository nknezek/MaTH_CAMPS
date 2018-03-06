function heat = H(t,pm)
% heat production each layer
heat = pm.H0.*exp(-pm.lambda*t);
end

