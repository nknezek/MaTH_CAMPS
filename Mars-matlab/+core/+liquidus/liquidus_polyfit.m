function liq_fun = liquidus_polyfit(wtpS)
% returns liquidus function taking in pressure (GPa) and returning
% liquidus temperature (K)
if (wtpS >= 13) && (wtpS <= 17)
    data_l = csvread(['+core/+liquidus/liquidus_',num2str(wtpS),'S_low.csv'],1);
    Tliq_l = data_l(:,1); % K
    Pliq_l = data_l(:,2); % GPa
    pfliq_l = polyfit(Pliq_l,Tliq_l,2);

    data_h = csvread(['+core/+liquidus/liquidus_',num2str(wtpS),'S_high.csv'],1);
    Tliq_h = data_h(:,1); % K
    Pliq_h = data_h(:,2); % GPa
    pfliq_h = polyfit(Pliq_h,Tliq_h,2);
    liq_fun = @(P) max(polyval(pfliq_l,P),polyval(pfliq_h,P));
else
    data = csvread(['+core/+liquidus/liquidus_',num2str(wtpS),'S.csv'],1);
    Tliq = data(:,1); % K
    Pliq = data(:,2); % GPa
    pfliq = polyfit(Pliq,Tliq,2);
    liq_fun = @(P) polyval(pfliq,P);
end