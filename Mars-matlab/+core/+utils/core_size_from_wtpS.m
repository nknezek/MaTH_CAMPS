function R_core = core_size_from_wtpS(wtpS)
% gives a core radius in (m), given wt% sulfur in the core
% values from Fei

wtpS0 = [0,5,10,15,20,25]; 
D_core0 = [2073.2, 2029.1, 1981.6, 1929.9, 1874.9, 1814.9];

pf = polyfit(wtpS0, D_core0, 2);
D_core = polyval(pf,wtpS)*1e3;
R_core = 3400e3-D_core;
end
