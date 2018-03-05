function liquidus()
% Low Sulfur Case
S_l = [0.,3.,6.,9.,12.];
P_l = [0.0,10.,14.,23.,40.];

S_14GPa_l = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,18.27,];
T_14GPa_l = [2219,2071,1954,1863,1793,1743,1708,1685,1669,1658,1648,1635,1615,1585,1542,1482,1400,1294,1160,1119,];
S_21GPa_l = [0.,3.,6.,9.,12.,15.4];
T_21GPa_l = [2271.,2185.,2081.,1949.,1758.,1348.];
S_23GPa_l = [0.,3.,6.,9.,12.,16.];
T_23GPa_l = [2308.,2211.,2106.,1985.,1819.,1435.];
S_40GPa_l = [0.,3.,6.,9.,12.];
T_40GPa_l = [2498.,2372.,2215.,1997.,1565.];

% High Sulfur Case
S_14GPa_mid = [19,20,20.8];
T_14GPa_mid = [1144,1160,1173];
S_14GPa_h = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36];
T_14GPa_h = [1236,1394,1471,1528,1576,1618,1656,1690,1722,1752,1781,1808,1834,1858,1882,1905];
S_21GPa_h = [16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,];
T_21GPa_h = [1370,1405,1428,1443,1454,1461,1467,1472,1477,1480,1532,1586,1637,1685,1731,1775,1818,1859,1898,1937,1974,];
S_23GPa_h = [16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,];
T_23GPa_h = [1435,1450,1461,1470,1485,1492,1555,1614,1669,1720,1768,1812,1852,1888,1921,1950,1976,1997,2015,2029,2040,];
S_30GPa_h = [14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,];
T_30GPa_h = [1500,1575,1646,1707,1759,1804,1842,1876,1906,1934,1958,1982,2003,2023,2042,2060,2077,2093,2109,2124,2138,2151,2164,];
S_40GPa_h = [13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,];
T_40GPa_h = [1647,1730,1787,1833,1873,1908,1941,1971,1999,2026,2051,2075,2098,2120,2142,2162,2182,2202,2220,2239,2256,2274,2291,2308,];

S_21GPa_h1 = S_21GPa_h(010);
S_21GPa_h2 = S_21GPa_h(10-1);
T_21GPa_h1 = T_21GPa_h(010);
T_21GPa_h2 = T_21GPa_h(10-1);
S_23GPa_h1 = S_23GPa_h(05);
S_23GPa_h2 = S_23GPa_h(5-1);
T_23GPa_h1 = T_23GPa_h(05);
T_23GPa_h2 = T_23GPa_h(5-1);

% Fit splines at each pressure
p14lsp = spline(S_14GPa_l,T_14GPa_l);
p21lsp = spline(S_21GPa_l,T_21GPa_l);
p23lsp = spline(S_23GPa_l,T_23GPa_l);
p30lsp = spline(S_30GPa_l,T_30GPa_l);
p40lsp = spline(S_40GPa_l,T_40GPa_l);
p14msp = spline(S_14GPa_mid,T_14GPa_mid);
p14hsp = spline(S_14GPa_h,T_14GPa_h);

p21hsp1 = spline(S_21GPa_h1,T_21GPa_h1);
p21hsp2 = spline(S_21GPa_h2,T_21GPa_h2);
p23hsp1 = spline(S_23GPa_h1,T_23GPa_h1);
p23hsp2 = spline(S_23GPa_h2,T_23GPa_h2);

p30hsp = spline(S_30GPa_h,T_30GPa_h);
p40hsp = spline(S_40GPa_h,T_40GPa_h);

%%
Pfit = [14.,21.,23.,30.,40.];
Tfit = [p14lsp(S),p21lsp(S),p23lsp(S),p30lsp(S),p40lsp(S)];
T_SPl = spline(Pfit, Tfit);

% High S side of eutectic
function T = T_p14h(S)
    if S<21.
        T =  p14msp(S);
    else
        T =  p14hsp(S);
    end
end
function T = T_p21h(S)
    if S<25.
        T =  p21hsp1(S);
    else
        T =  p21hsp2(S);
    end
end
function T = T_p23h(S)
    if S<21.
        T =  p23hsp1(S);
    else
        T =  p23hsp2(S);
    end
end
function T = T_p30h(S)
    T = p30hsp(S);
end
function T = T_p40h(S)
    T = p40hsp(S);
end

function T = T_SPh(S,P)
    Pfit = [14.,21.,23.,30.,40.];
    Tfit = [T_p14h(S),T_p21h(S),T_p23h(S),T_p30h(S),T_p40h(S)];
    spl = spline(Pfit, Tfit);
    T = spl(P);
end 

function T_SP(S,P)
    P = P/1e9;
    S = S*100;
    assert(S>=0. && S<=35.)
    assert(P>=0. && P<=45.)
    if S<12.
        T = T_SPl(S,P);
    elseif S>18.27
        T = T_SPh(S,P);
    else
        T =  max(T_SPl(S,P),T_SPh(S,P));
    end
end






end

