function [dT,flm,flv,Crt]=migration(f,Ti,pm)

% Calculates mass, energy, volume change in lithosphere layer

% Note: As of now, assume all melt (f) is carried into lithosphere and
% crystalizes there, releasing heat.

%Note: assume all excess T goes into melt.

% ADD adiabat for melt to lower temperature through ascendance
% ADD Melt migration and freezing energy loss to calculate dT of layer
% ADD bouyancy of melt (neutral level)

n = pm.n;
Cp = pm.Cp;
rho = pm.rho;
M = pm.M;
L = pm.L;
R = pm.R;
Tsol = pm.Tsol;
Tliq = pm.Tliq;

%% Migration post proc
fm=zeros(size(f));         % melt mass matrix in kg
flv=zeros(size(f));         % melt volume in km^3
Q0=zeros(size(fm));         % energy matrix from crystallize material within lid, in J

for m=1:n
    fm(:,m)=f(:,m)*M(m);
end

%Add in Cpm of 1500 from di genova et al 2014
%Qo is incomplete, change in T from superheated melt, change in T from
%phase change, change in T from subliquiudius to subsolidus?, make sure Q0 from latent heat of crystalizaion is ~ same Q0 from melt - seems to not be at moment. 

for m=1:n
%   Q0(:,1)=Q0(:,1)+(Tsol(1)-Tliq(m))*Cp(m)*fm(:,m);      % take mass of melt from layer, + melt temperature in layer, move to lid and crystallize, releasing energy to layer [this formulation may not be right, need melt temperature]
    Q0(:,1)=Q0(:,1)+L(m)*fm(:,m); %+(Tliq(m)-Tsol(1))*Cp(m)*fm(:,m); %latent heat released from crystallization of melt, energy transferred from melt temperature decrease to lid solidus  -- No adaiabet inclusded yet
end

dT=Q0./(M(1)*Cp(1));     % temperature change in lid

if (dT(1) + Ti(1)) > Tsol(1)
    Tt=(dT(1) + Ti(1))-Tsol(1);     % dT above solidus
    f0 = (Tt-Tsol(1))/(Tliq(1)-Tsol(1));
       
    %This part needs some work ... think on
    if f0 > 1
        f0=Q0(1)/(M(1)*L(1)*f0);                     % assume entire layer melts, rest of energy is lost radiatively
    else
        f0=f + Q0(1)/(M(1)*L(1)*f0);
    end
end

if (dT(1) + Ti(1)) < Tsol(1)
    f0=0;
end


Ti(1)=Tsol(1);

fm(1,1)=f0*M(1);     % this sometimes causes problems as f0 may apparently not always exist --- FIX!

flm=fm(:,1);
for m=2:n
    flm=flm+fm(:,m);             %all mass of melt into layer 1 - "lid"
end

flv=flm/(rho(1))/(1e3)^3;        % convert to volume km^3

Crt=0.2*flv/(pi*(R(1)/1e3)^2);   %crustal thickness from extrusion km^3

end