%% Main Script for MaTH_CAMPS
% Runs everything. 
% Matt Weller, cool dude.

%% Added Melt and Migration

%%Todo
%       
%       add    Fix melt in lid temperature change
%               temperature dependant viscosity (be aware of 1e-16 rounding errors, will crop up in temp dep visc )
%%
clear all

R=[3400;3300;2200;1830;100]*1e3; %Radius boundaries between layers; m
%rho=[3.4; 3.4; 5.066; 7]*1e3; % density; kg/m^3
rho=[3.2;3.4;4.0;8]*1e3; % density of layers, from mix run

Ts=250; % surface temperature, K

time_end=4600;   %How long to run time series in Mys

%mark=[[Ts Ta(end,[2:4])]' R(1:4)];

%partition coeff for radiogenics
%al = 0; be = 0; ga = 0;
%al = 1; be = 1; ga = 1;
%al = 1.95; be = 0.05; ga = 1.;
al = 1.25; be = 0.25; ga = 1.;
%al = 2.45; be = 0.05; ga = 0.5;

H0=[al*8.16e-8; be*8.16e-8; ga*8.16e-8; 0]./rho; %heat production W/m^3 (input)
lambda=[1.37e-17;1.37e-17; 1.37e-17; 0];         % decay constants

%Temperature for non-magma melt ocean
%Tbase=[2100 +- 200 K]
% Hot is 1100 - 2100 K
% Cold is 1100 - 1700 K


%% Geometry
A=4*pi*R.^2; %Surface Area of layer boundaries; m^2
Vc=4*pi/3.*R.^3; %cumulative volume of layer; m^3;
V=diff(-Vc); %volume in each layer
n=numel(V); %number of layers

%%
% Material properties

%Longhi 1992  Mars bible
%Tsol= -2.1474E-10*((R(1)-R)/1e3).^4 + 1.0833E-06*((R(1)-R)/1e3).^3 - 0.0018821*((R(1)-R)/1e3).^2 + 1.6835*((R(1)-R)/1e3) + 1304.8;
%Tliq= -1.22E-10*((R(1)-R)/1e3).^4 + 5.763E-07*((R(1)-R)/1e3).^3 - 0.00088614*((R(1)-R)/1e3).^2 + 8.0731E-01*((R(1)-R)/1e3) + 1765;

%Combination of SFB (LPSC abstract '01) 
Tsol= 2.02E-07*((R(1)-R)/1e3).^3 - 0.001*((R(1)-R)/1e3).^2 + 1.68E+00*((R(1)-R)/1e3) + 1.34E+03;
Tliq= 9.36E-08*((R(1)-R)/1e3).^3 - 5.28e-4*((R(1)-R)/1e3).^2 + 1.06E+00*((R(1)-R)/1e3) + 1.79E+03;

Ta0=[Tliq(1), Tliq(2), Tsol(3)-250, Tsol(4)+25];

temp=ones(size(V));mu=ones(size(V)); %mu is conversion factor in each layer (vector); STO table 13.2
Rac=658*temp; %critical Rayleigh number source Breuer and Moore 2008

k=3*temp; %thermal conductivity; W/m/K
k(end-1)=4.8;  % Sets k for lower layer. Note: Make sure core is basal unit!
k(end)=40; % thermal conductivity core; W/m/K 
K=(1e-6)*temp; %thermal diffusivity; m^2/s
K(end-1)=0.5e-6; % IMPORTANT: ONLY WORKS FOR 4-layer case
eta=0.1e21*temp; %viscosity structure Pa.s
eta(end-1)=1e16; %viscosity lower layer Pa.s
%eta(end)=1e16; %viscosity of core Pa.s
eta(1)=1e21;     %viscosity  of 'lid' Pa.s
L=6e5*temp;   % Latent heat of fusion J/Kg
L(end)=6e5;   % Latent heat of fusion J/Kg for core -- 600 - 1200 kJ/kg Poirier (1994) Anderson & Duba (1997)

nu=eta./rho; %kinematic viscosity; m^2/s
Cp=k./(rho.*K).*temp; %heat capacity;
Cp(end-1)=1000;       %heat capacity lower layer Pa.s
Cp(end)=800;   % heat capicity of core -- 800 � 80 J/(kg K) Poirier (1994)

a=(2e-5)*temp; %thermal expansion; 1/K
g=3.71*temp; %acceleration of gravity; m/s^2;
beta=(1/3)*temp; %parameter; Nu=Ra^(1/beta);
gamma=1./(1+beta); 
Theta=(k.*((a.*g./(K.*nu.*Rac)).^beta)); %temperature drop to heat flux conversion factor
M=rho.*V; %mass in each layer
Myr=1e6*3e7; %convert s to Myr

% dependent variable:
%% Early melt processing for solid state convection 
%Ta has time series of average temperature in each layer

Tuf=@(Ta)Ta./mu;
Tlf=@(Ta)2*Ta-Tuf(Ta);
H=@(t)H0.*exp(-lambda*t);

% Calculate melt from initial conditions to set initial solid state
% temperature field for ODE solver
[f0,dTm0]=melt_ini(Tsol,Tliq,Ta0,Cp,L,M,n);
[dTli,flmi,flvi,Crti,fli]=migration_ini(f0,n,Cp,Tsol,Tliq,Ta0(1),rho,M,L,R);   % calculates inital melting in lid
f0(1)=fli(1);    % induced melt from crystaliztion in lid

% Set initial temperatures to solidus approximation (solid state convection)
for x=1:n
    Ta0(x)=[Ta0(x)-dTm0(x)];
end

% Solve the system the first time, before melt processing
% [t,Ta]=ode45(@(t,Ta)convectionODE(t,Ta,Tuf,Tlf,Ts,Theta,gamma,A,M,H,Cp,n,R), [0/Myr,time_end*Myr],Ta0);
%%
dt = 1; %time step in Myr
t = 1:dt:time_end;
Ta = zeros(length(t),n);
t = t'*Myr; % Convert to seconds
dt = dt*Myr; % Convert to seconds
Ta(1,:) = Ta0;
%%
for i=1:(length(t)-1)
    dTa_dt = convectionODE(t(i),Ta(i,:)',Tuf,Tlf,Ts,Theta,gamma,A,M,H,Cp,n,R);
    Ta(i+1,:) = Ta(i,:) + dTa_dt'*dt;
end
% [t,Ta]=dTadt = convectionODE(t,Ta,Tuf,Tlf,Ts,Theta,gamma,A,M,H,Cp,n,R), [0/Myr,time_end*Myr],Ta0);


%% Melting processing
% loops and solves equations over time frame of melt generation

% nt=size(Ta,1);
% ThetaMat=repmat(Theta',[nt,1]);


[f,dTm]=melt(Tsol,Tliq,Ta,Cp,L,M,n);
Ta=Ta-dTm;


% Splice melt fractions for different steps
for x=1:n
    f(1,x)=f(1,x)+f0(1,x);
end


%if f is generated from melt in time step (j, where j>t1), use time step j as next Ta0 (solidus temp) in solver.  
% workflow --> ode -> to melt -> is f generated? -> yes, take time step j
% from f, use as new Ta0 in solver, resolve -> no, continue as is. 

ismelt=f>0;           %Finds where melt occurs 1: yes 2: no

for m=2:(n-1)
    ismelt(:,1)=ismelt(:,1)+ismelt(:,m);             %collapses all melt to ismelt(:,1) 
end

ind=find(ismelt(:,1)>0);
tmelt=t(ind)/Myr;              %time series of melt to rerun solver and melt modules - in Myr's


% Loop to take time step at every f, melt, resolve Ta from ode -- easy
% peasy --> eg., [t(10)/Myr,time_end*Myr],Ta0);, where Ta0 would now be the last
% solved f step, then splice times togther -- WATCH INDEXING
if length(ind) > 1
    for m=2:length(tmelt)
        [t1,Ta1]=ode45(@(t,Ta)convectionODE(t,Ta,Tuf,Tlf,Ts,Theta,gamma,A,M,H,Cp,n,R), [tmelt(m)*Myr,time_end*Myr],Ta(m-1,:));  % uses Temperatures from last time step Ta(m-1) as inital condition for new melt time step tmelt(m)
         [f1,dTm1]=melt(Tsol,Tliq,Ta1,Cp,L,M,n);
         Ta1=Ta1-dTm1;
%           for x=1:n
%               f1(1,x)=f1(1,x)+f0(1,x);
%           end
          % Modify t and Ta such that the new output t is spliced into
          % orginal t
          t=[t(1:m-1);t1(2:end)];                      % new time vector running from last melt instance
          Ta=[Ta(1:m-1,:);Ta1(2:end,:)];               % new Temperature matrix running from last melt instance
          f=[f(1:m-1,:);f1(2:end,:)];                  %new f matrix running from last melt instance    
    end
end

[dTl,flm,flv,Crt]=migration(f,n,Cp,Tsol,Tliq,Ta,rho,M,L,R);   % need to fix, thinks first time step is 8000k
%Ta=dTl+Ta;

%% Postprocessing
nt=size(Ta,1);
% more temperatures
Tu=Ta./repmat(mu',[nt,1]); %base of the upper BL
Tl=2*Ta-Tu; %top of the lower BL
% Boundary temperatures
Tb=NaN(nt,n+1);
Tb(:,1)=Ts; % surface temperature
ThetaMat=repmat(Theta',[nt,1]);

for i=2:n;
    Tb(:,i)=(ThetaMat(:,i).*Tu(:,i)+ThetaMat(:,i-1).*Tl(:,i-1))./(ThetaMat(:,i)+ThetaMat(:,i-1));
 
end

Tb(:,n+1)=Tl(:,end); %made-up temperature at the center of the planet;



%% temperature post-proc

%Temperature drop (K)
DTu=Tu-Tb(:,1:end-1);
DTl=Tb(:,2:end)-Tl;
% Heat flux (W/m^2)
Qu=sign(DTu).*ThetaMat.*(abs(DTu)).^(4/3);
Ql=sign(DTl).*ThetaMat.*(abs(DTl)).^(4/3);
% boundary Layer thickness
kMat=repmat(k',[nt,1]);
du=kMat.*DTu./Qu;
dl=kMat.*DTl./Ql;

%buoyancy number between bottom layer and UM  --> B =
%(rhol-rho0)/(alpha*rho0*delT)
B=(rho(end-1)-rho(end-2))./(a(end-2)*rho(end-2)*(Tb(:,end)-Tb(:,1)));

%System Ra (uses averages) ****NOTE: Currently setup for n=3,4 system
Llr=(R(n-1)-R(n))/((R(n-2)-R(end-1)));      %percent of lower layer to entire convecting mantle (Radial)
Ulr=(R(n-2)-R(n-1))/((R(n-2)-R(end-1)));    %percent of upper layer to entire convecting mantle (Radial)
rho_a=Llr*rho(end-1)+Ulr*rho(end-2);
eta_a=Llr*eta(end-1)+Ulr*eta(end-2);

Ra=zeros(length(Tb),2);
for i=1:length(Tb);
    Ra(i,1)=mean(a)*rho_a*mean(g)*(Tb(i,end-1)-Tb(i,1))*(R(1)-R(end-1))^3/(mean(K)*eta_a);            % system Ra
    Ra(i,2)=a(n-1)*rho(n-1)*g(n-1)*(Tb(i,end-1)-Tb(i,end-2))*(R(n-1)-R(n))^3/(K(n-1)*eta(n-1));           % ll Mantle Ra
    if max(Tb(i,:))-min(Tb(i,:)) < 700;
        Ra(i,1)=1; Ra(i,2)=1;
    end
    if Ra(i,2) < 0
        Ra(i,2)=NaN;  %temperature inversion
    end
end

Ru=repmat([R(1:end-1)]',[nt,1])-du; %radius at the bottom of the upper BL
Rl=repmat([R(2:end)]',[nt,1])+dl; %radius at the top of the lower BL

HMat=repmat(H0',[nt,1]).*exp(-repmat(lambda',[nt,1]).*repmat(t,[1,n]));
VMat=repmat(V',[nt,1]);
rhoMat=repmat(rho',[nt,1]);
Ht=sum(HMat.*VMat.*rhoMat,2);

%%
Tmat=NaN(nt,3*n+1);
Rmat=NaN(nt,3*n+1);
for i=1:n;
    Tmat(:,i*3-2)=Tb(:,i);
    Tmat(:,i*3-1)=Tu(:,i);
    Tmat(:,i*3-0)=Tl(:,i);
    Rmat(:,i*3-2)=repmat(R(i),[nt,1]);
    Rmat(:,i*3-1)=Ru(:,i);
    Rmat(:,i*3-0)=Rl(:,i);
end
Rmat(:,n*3+1)=repmat(R(n+1),[nt,1]);
Tmat(:,i*3+1)=Tb(:,n+1);


%% Plotting
Npl = 1000;
plotting(t(1:Npl:end),Myr,Qu(1:Npl:end,:),n,Ht(1:Npl:end,:),A, ...
Ta(1:Npl:end,:),flv(1:Npl:end,:),Tmat(1:Npl:end,:),Rmat(1:Npl:end,:), ...
Ra(1:Npl:end,:),Crt(1:Npl:end,:),f(1:Npl:end,:),B(1:Npl:end,:),Tb(1:Npl:end,:))


%% Dynamo shutoff time based on termination heat output of ~ 1-2 Tw -- 20-60 mw/m^2


Dynamoh=find(Qu(:,end)*1e3 > 60); %Dynamo shutoff [High end]
Dynamol=find(Qu(:,end)*1e3 > 15); %Dynamo shutoff [low end]

if isempty(Dynamoh) || isempty(Dynamol)   
    display ('No upper or lower bound possible for dynamo')
else
    if Dynamoh < 2
        display(['Thermally Driven Dynamo Unlikely'])
    end
    if length(t) >= (Dynamol(end)+1)
        DTimeh=4.45-t(Dynamoh(end)+1)/Myr/1000;   % dynamo shut down time [High end]
        DTimel=4.45-t(Dynamol(end)+1)/Myr/1000;   % dynamo shut down time [low end]
        display(['Dynamo Critical at ' num2str(DTimeh) ' to ' num2str(DTimel) ' Bya'])
    else
        display(['Dynamo trucking along at' ' ' num2str(t(end)/Myr) ' Myr'])
    end
end

%% Eruption rate

display(['Cumulative Eruptive volume of' ' ' num2str(0.2*sum(flv)/1e10) 'e10 km^3'])     % assumes 80% intrusive, 20% extrusive
display(['Cumulative Eruptive volume in' ' ' num2str(0.2*sum(flv(2:end))/(pi*(R(1)/1e3)^2*0.25*100*.2)) '  '  'Tharsis Volumes (ignoring initial volume of crust formation)'])     % assumes 80% intrusive, 20% extrusive
display(['Cumulative Crustal Thickness ' ' ' num2str(sum(Crt)) ' km'])     % assumes 80% intrusive, 20% extrusive
display(['Cumulative Crustal Thickness ' ' ' num2str(sum(Crt(2:end))) ' km'   '  (ignoring initial volume of crust formation)'])     % assumes 80% intrusive, 20% extrusive

display(['Cumulative Eruption Rate of ' ' ' num2str(0.2*sum(flv)/((t(end)/Myr)*1e6)) ' km^3 yr^-1'])     % assumes 80% intrusive, 20% extrusive

inda=find (t/Myr > (time_end-3e3)); indh=find(t/Myr < (time_end-3e3) & t/Myr > (time_end-3.8e3)); indn=find(t/Myr < (time_end-3.8e3));
Amaz=flv(inda); Hesp=flv(indh); Noac=flv(indn); 

display(['Noachian Eruption Rate of ' ' ' num2str(0.2*sum(flv(indn))/((t(indn(end))/Myr)*1e6)) ' km^3 yr^-1'])     % assumes 80% intrusive, 20% extrusive 
display(['Hesperian Eruption Rate of ' ' ' num2str(0.2*sum(flv(indh))/((t(indh(end))/Myr)*1e6)) ' km^3 yr^-1'])     % assumes 80% intrusive, 20% extrusive
display(['Amazonian Eruption Rate of ' ' ' num2str(0.2*sum(flv(inda))/((t(inda(end))/Myr)*1e6)) ' km^3 yr^-1'])     % assumes 80% intrusive, 20% extrusive

