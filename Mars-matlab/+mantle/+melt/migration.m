function [dT_lid, melt_mass_lid, melt_volume_lid, thickness_lid]=migration(f,Ti,pm)

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

%% Migration during model
Q0 = zeros(1,n);         % energy matrix from crystallize material within lid, in J
melt_mass = zeros(1,n);         % melt mass matrix in kg

for m=1:n
    melt_mass(m) = f(m)*M(m);
end

%Add in Cpm of 1500 from di genova et al 2014
%Qo is incomplete, change in T from superheated melt, change in T from
%phase change, change in T from subliquiudius to subsolidus?, make sure Q0 from latent heat of crystalizaion is ~ same Q0 from melt - seems to not be at moment. 

% account for energy release of crystallization of melt in
% crust/lithosphere layer
for m=1:n 
    Q0(1)=Q0(1)+L(m)*melt_mass(m); %+(Tliq(m)-Tsol(1))*Cp(m)*melt_mass(:,m); %latent heat released from crystallization of melt, energy transferred from melt temperature decrease to lid solidus  -- No adiabat included yet
end

dT_lid = Q0(1)./(M(1)*Cp(1));     % temperature change in lid from melt crystallization

% how much crust do you melt from the transported melt/heat from below
if (dT_lid + Ti(1)) > Tsol(1)
    Tt=(dT_lid + Ti(1))-Tsol(1);     % dT above solidus
    f0 = (Tt-Tsol(1))/(Tliq(1)-Tsol(1));
    
    %This part needs some work ... think on
    if f0 > 1
        f0=Q0(1)/(M(1)*L(1)*f0);                     % assume entire layer melts, rest of energy is lost radiatively
    else
        f0=f(1) + Q0(1)/(M(1)*L(1)*f0);
    end
elseif (dT_lid(1) + Ti(1)) <= Tsol(1)
    f0=0;
end
Ti(1)=Tsol(1);

melt_mass_lid = f0(1)*M(1);     % this sometimes causes problems as f0 may apparently not always exist --- FIX!
for m=2:n
    melt_mass_lid = melt_mass_lid + melt_mass(:,m);             %all mass of melt into layer 1 - "lid"
end

melt_volume_lid = melt_mass_lid/(rho(1))/(1e3)^3;        % convert to volume km^3
thickness_lid= 0.2*melt_volume_lid/(pi*(R(1)/1e3)^2);   %crustal thickness from extrusion km^3

end