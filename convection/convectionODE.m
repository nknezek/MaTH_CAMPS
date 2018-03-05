function dTadt=convectionODE(t,Ta,Tuf,Tlf,Ts,Theta,gamma,A,M,H,C,n,R);
%%
th=Theta.^gamma;
R=(R(1)-R)/1e3;   % flip R(m --> km) so that surface is 0 for solidus/liquidus calculation
% Ta=@(T,mu)mu.*T; %average temperature in each layer (vector)
% Tl=@(Tu)Tu; %temperature at the top of the bottom boundary layer (vector)

% Ta=mu.*Tu;
Tu=Tuf(Ta);
Tl=Tlf(Ta);
% Determine boundary temperatures
Tb=Tl.*NaN;
Tb(1)=Ts; % surface temperature
%f(1) = (Tb(1)-Tsol(1))/(Tliq(1)-Tsol(1)); %surface melt
for i=2:n
    Tb(i)=(th(i).*Tu(i)+th(i-1).*Tl(i-1))./(th(i)+th(i-1));
    %f(i) = (Tb(i-1)-Tsol(i-1))/(Tliq(i-1)-Tsol(i-1));    % melt fraction
end


%F=horzcat(f,delTm)

Tb(n+1)=Tl(end); %made-up temperature at the center of the Earth;	

dTu=Tu-Tb(1:end-1);
dTl=Tb(2:end)-Tl;

qt=sign(dTu).*Theta.*(abs(dTu)).^(4/3);
qb=sign(dTl).*Theta.*(abs(dTl)).^(4/3);

dTadt=((A(2:end).*qb... heat flux from next layer below
    -A(1:end-1).*qt... heat flux toward next layer above
    )./M... spread over the volume of each layer
    +H(t)... adding heat production W/kg; assumed constant layer
    )./C;... convert to temperature change
    
% dTadt=(H+(...
%     (A(2:end).*qb... Heat flux from next layer below
%     -A(1:end-1).*qt) ... Heat flux toward next layer above
%     ./M)./C);