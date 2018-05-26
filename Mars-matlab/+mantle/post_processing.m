function pp = post_processing(t,Tvec,f,pm,pc)
% Perform various post-processing routines after a run of the thermal
% evolution code
pp = struct;
pp.tvec = t;
pp.Tvec = Tvec;
pp.f = f;
pp.t_plt = t/pm.Myr;

rho = pm.rho;
eta = pm.eta;
R = pm.R;
a = pm.a;
n = pm.n;
Ts = pm.Ts;
g = pm.g;
K = pm.K;
Myr = pm.Myr;
time_end = t(end)/Myr;

Tm = Tvec(:,1:pm.n);
pp.Tm = Tm;
pp.Tc = Tvec(:,pm.n+1:end);

% migrate melt to the lid to recrystallize (assume all moves to lid) 
% only used for melt production calculations, not temp change
[dTl,flm,flv,Crt]=mantle.melt.migration(f,Tm,pm);   % need to fix, thinks first time step is 8000k
%Tm=dTl+Tm;
pp.Crt = Crt;
pp.flv = flv;
pp.dTl = dTl;
pp.flm = flm;
% ======================
% temperature post-proc
% ======================
%% Calculate full temperatures within layer at all timesteps (not reported
% from ODE)
nt=size(Tm,1);
Tu=Tm./repmat(pm.mu',[nt,1]); %base of the upper BL
Tl=2*Tm-Tu; %top of the lower BL
Tb=NaN(nt,n+1); % Boundary temperatures
Tb(:,1)=Ts; % surface temperature
ThetaMat=repmat(pm.Theta',[nt,1]);
for i=2:n
    Tb(:,i)=(ThetaMat(:,i).*Tu(:,i)+ThetaMat(:,i-1).*Tl(:,i-1))./(ThetaMat(:,i)+ThetaMat(:,i-1));
end

Tb(:,n+1)=Tl(:,end); %made-up temperature at the center of the planet; only used in post-processing
pp.Tu = Tu;
pp.Tl = Tl;
pp.Tb = Tb;

%Temperature drop (K)
DTu=Tu-Tb(:,1:end-1);
DTl=Tb(:,2:end)-Tl;
pp.DTu = DTu;
pp.DTl = DTl;
% Heat flux (W/m^2) across layer boundaries
Qu=sign(DTu).*ThetaMat.*(abs(DTu)).^(4/3);
Ql=sign(DTl).*ThetaMat.*(abs(DTl)).^(4/3);
pp.Qu = Qu;
pp.Ql = Ql;
% boundary Layer thickness
kMat=repmat(pm.k',[nt,1]);
du=kMat.*DTu./Qu;
dl=kMat.*DTl./Ql;
pp.du = du;
pp.dl = dl;

% Calculate temperature at each layer boundary and thermal boundary layer Temp and radius
Ru=repmat([R(1:end-1)]',[nt,1])-du; %radius at the bottom of the upper BL
Rl=repmat([R(2:end)]',[nt,1])+dl; %radius at the top of the lower BL
pp.Ru = Ru;
pp.Rl = Rl;

Tmat=NaN(nt,3*n+1);
Rmat=NaN(nt,3*n+1);
for i=1:n
    Tmat(:,i*3-2)=pp.Tb(:,i);
    Tmat(:,i*3-1)=pp.Tu(:,i);
    Tmat(:,i*3-0)=pp.Tl(:,i);
    Rmat(:,i*3-2)=repmat(R(i),[nt,1]);
    Rmat(:,i*3-1)=Ru(:,i);
    Rmat(:,i*3-0)=Rl(:,i);
end
Rmat(:,n*3+1)=repmat(R(n+1),[nt,1]);
Tmat(:,i*3+1)=Tb(:,n+1);
pp.Rmat = Rmat;
pp.Tmat = Tmat;
%% buoyancy number between bottom layer and UM  --> 
% e.g. B = (rhol-rho0)/(alpha*rho0*delT)
pp.B=(rho(end-1)-rho(end-2))./(a(end-2)*rho(end-2)*(Tb(:,end)-Tb(:,1)));

%% System Ra (uses averages) ****NOTE: Currently setup for n=3,4 system
Llr=(R(n-1)-R(n))/((R(n-2)-R(end-1)));      % percent of lower layer to entire convecting mantle (Radial)
Ulr=(R(n-2)-R(n-1))/((R(n-2)-R(end-1)));    % percent of upper layer to entire convecting mantle (Radial)
rho_a=Llr*rho(end-1)+Ulr*rho(end-2);
eta_a=Llr*eta(end-1)+Ulr*eta(end-2);

Ra=zeros(length(Tb),2);
for i=1:length(Tb)
    Ra(i,1)=mean(a)*rho_a*mean(g)*(Tb(i,end-1)-Tb(i,1))*(R(1)-R(end-1))^3/(mean(K)*eta_a);            % system Ra
    Ra(i,2)=a(n-1)*rho(n-1)*g(n-1)*(Tb(i,end-1)-Tb(i,end-2))*(R(n-1)-R(n))^3/(K(n-1)*eta(n-1));           % ll Mantle Ra
    if max(Tb(i,:))-min(Tb(i,:)) < 700
        Ra(i,1)=1; Ra(i,2)=1;
    end
    if Ra(i,2) < 0
        Ra(i,2)=NaN;  %temperature inversion
    end
end
pp.Ra = Ra;

%% radiactive heating in each layer over time
HMat=repmat(pm.H0',[nt,1]).*exp(-repmat(pm.lambda',[nt,1]).*repmat(t,[1,n]));
VMat=repmat(pm.V',[nt,1]); % volume 
rhoMat=repmat(pm.rho',[nt,1]); % density

pp.Ht=sum(HMat.*VMat.*rhoMat,2); % total heat production

%% Eruption rate
try
display(['Cumulative Eruptive volume of' ' ' num2str(0.2*sum(flv)/1e10) 'e10 km^3'])     % assumes 80% intrusive, 20% extrusive
display(['Cumulative Eruptive volume in' ' ' num2str(0.2*sum(flv(2:end))/(pi*(R(1)/1e3)^2*0.25*100*.2)) '  '  'Tharsis Volumes (ignoring initial volume of crust formation)'])     % assumes 80% intrusive, 20% extrusive
display(['Cumulative Crustal Thickness ' ' ' num2str(sum(Crt)) ' km'])     % assumes 80% intrusive, 20% extrusive
display(['Cumulative Crustal Thickness ' ' ' num2str(sum(Crt(2:end))) ' km'   '  (ignoring initial volume of crust formation)'])     % assumes 80% intrusive, 20% extrusive

display(['Cumulative Eruption Rate of ' ' ' num2str(0.2*sum(flv)/((t(end)/Myr)*1e6)) ' km^3 yr^-1'])     % assumes 80% intrusive, 20% extrusive

inda=find (t/Myr > (time_end-3e3)); 
indh=find(t/Myr < (time_end-3e3) & t/Myr > (time_end-3.8e3)); 
indn=find(t/Myr < (time_end-3.8e3));
Amaz=flv(inda); Hesp=flv(indh); Noac=flv(indn); 

display(['Noachian Eruption Rate of ' ' ' num2str(0.2*sum(flv(indn))/((t(indn(end))/Myr)*1e6)) ' km^3 yr^-1'])     % assumes 80% intrusive, 20% extrusive 
display(['Hesperian Eruption Rate of ' ' ' num2str(0.2*sum(flv(indh))/((t(indh(end))/Myr)*1e6)) ' km^3 yr^-1'])     % assumes 80% intrusive, 20% extrusive
display(['Amazonian Eruption Rate of ' ' ' num2str(0.2*sum(flv(inda))/((t(inda(end))/Myr)*1e6)) ' km^3 yr^-1'])     % assumes 80% intrusive, 20% extrusive
catch
    
end

end