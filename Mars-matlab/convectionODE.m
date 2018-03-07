function dTalldt = convectionODE(t,Tall,pm,pc)
%% TODO: Add in core code in this function
Tm = Tall(1:pm.n);
Tc = Tall(pm.n+1:end);
Ts = pm.Ts;
A = pm.A;
M = pm.M;
C = pm.Cp;
n = pm.n;

th = pm.Theta.^pm.gamma;
R = (pm.R(1)-pm.R)/1e3;   % flip R(m --> km) so that surface is 0 for solidus/liquidus calculation

% Determine boundary layer temperatures
Tu = mantle.utils.Tuf(Tm,pm);
Tl = mantle.utils.Tlf(Tm,pm);
Tb = Tl.*NaN;
Tb(1) = Ts; % surface temperature

% convect
for i=2:n
    Tb(i)=(th(i).*Tu(i)+th(i-1).*Tl(i-1))./(th(i)+th(i-1));
end

Tb(n+1)=Tc(end); % CMB temperature from core code

% change in temperature across boundary layers
dTu=Tu-Tb(1:end-1);
dTl=Tb(2:end)-Tl;

% heat flow across boundary layers [W/m^2]
qt=sign(dTu).*pm.Theta.*(abs(dTu)).^(4/3);
qb=sign(dTl).*pm.Theta.*(abs(dTl)).^(4/3);

dTmdt=((A(2:end).*qb... heat flux from next layer below
    -A(1:end-1).*qt... heat flux toward next layer above
    )./M... spread over the mass of each layer
    +mantle.utils.H(t,pm)... adding heat production W/kg; assumed constant layer
    )./C;... convert to temperature change

% call core processing
dt = pc.dt;
dTcdt = core.therm.dconvect_and_conduct_dt(Tc', dt, qb(end)*A(end), pc);

dTalldt = [dTmdt',dTcdt]';
end
    
