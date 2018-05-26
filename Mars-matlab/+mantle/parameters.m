function pm = parameters()
%% Parameters for Martian Mantle
% Nicholas Knezek
%
% returns
% pm = parameters struct of mantle parameters

pm = struct;

% Number of layers for thermal convection (crust, mantle, lower mantle)
pm.n = 3; % [-]
ones_vec = ones(pm.n,1);

% Geometry
pm.R = [3400; 3300; 2200; 1830]*1e3; %Radius boundaries between layers; m
pm.A = 4*pi*pm.R.^2; %Surface Area of layer boundaries; m^2
Vc = 4*pi/3.*pm.R.^3; %cumulative volume of layer; m^3;
pm.V = diff(-Vc); %volume in each layer

pm.rho=[3.2; 3.4; 4.0]*1e3; % density of layers, from mix run
pm.M = pm.rho.*pm.V; %mass in each layer
pm.g = 3.71*ones_vec; % [m/s^2] - acceleration of gravity

% Thermal Parameters
% mu is conversion factor in each layer between volume average and upper/lower thermal boundary layeres (vector); 
% STO table 13.2 -- dunno what STO is?? (Nick)
pm.Ts = 250; % surface temperature, K
pm.mu = ones_vec; 
pm.Rac = 658*ones_vec; % [-] - critical Rayleigh number (Breuer and Moore 2008)
pm.k = 3*ones_vec; % [W/m-K] - thermal conductivity
pm.k(3) = 4.8;  % [W/m-K] - thermal conductivity for lower layer. 
pm.K = (1e-6)*ones_vec; % [m^2/s] - thermal diffusivity; 
pm.K(3) = 0.5e-6; % IMPORTANT: ONLY WORKS FOR 4-layer case

% Viscosity
pm.eta = 0.1e21*ones_vec;  % [Pa-s] - viscosity structure 
pm.eta(3) = 1e16; % [Pa-s] - viscosity lower layer 
pm.eta(1) = 1e21;     % [Pa-s] - viscosity  of 'lid' 
pm.nu = pm.eta./pm.rho; % [m^2/s] - kinematic viscosity for each layer; 

% heat 
pm.L = 6e5*ones_vec;   % [J/kg] - Latent heat of fusion 
pm.Cp = pm.k./(pm.rho.*pm.K); % [J/kg-K] heat capacity
pm.Cp(3) = 1000; % [J/kg-K] heat capacity lower layer 
pm.a = (2e-5)*ones_vec; % [1/K] - thermal expansivity (alpha)

% == convection parameters
pm.beta = (1/3)*ones_vec; %parameter; Nu=Ra^(1/beta);
pm.gamma = 1./(1+pm.beta); 
pm.Theta = (pm.k.*((pm.a.*pm.g./(pm.K.*pm.nu.*pm.Rac)).^pm.beta)); %temperature drop to heat flux conversion factor
% ==
pm.Myr = 365.25*24*60*60*1e6; %convert s to Myr

% ======= Radiogenics ====== 
% partition coeff for radiogenics (al = crust, be = mantle, 
% ga = lower mantle layer
%al = 0; be = 0; ga = 0;
%al = 1; be = 1; ga = 1;
%al = 1.95; be = 0.05; ga = 1.;
al = 1.25; be = 0.25; ga = 1.;
%al = 2.45; be = 0.05; ga = 0.5;

pm.H0 = [al*8.16e-8; be*8.16e-8; ga*8.16e-8]./pm.rho; %heat production W/m^3 (input)
pm.lambda = [1.37e-17;1.37e-17; 1.37e-17];         % decay constants

% ===== Solidus and Liquidus =======
%Longhi 1992  Mars bible
%Tsol= -2.1474E-10*((R(1)-R)/1e3).^4 + 1.0833E-06*((R(1)-R)/1e3).^3 - 0.0018821*((R(1)-R)/1e3).^2 + 1.6835*((R(1)-R)/1e3) + 1304.8;
%Tliq= -1.22E-10*((R(1)-R)/1e3).^4 + 5.763E-07*((R(1)-R)/1e3).^3 - 0.00088614*((R(1)-R)/1e3).^2 + 8.0731E-01*((R(1)-R)/1e3) + 1765;

%Combination of (Schmerr, Fei, Bertka 2001 LPSC abstract) 
%R = pm.R;
%pm.Tsol = 2.02E-07*((R(1)-R)/1e3).^3 - 0.001*((R(1)-R)/1e3).^2 + 1.68E+00*((R(1)-R)/1e3) + 1.34E+03;
%pm.Tliq = 9.36E-08*((R(1)-R)/1e3).^3 - 5.28e-4*((R(1)-R)/1e3).^2 + 1.06E+00*((R(1)-R)/1e3) + 1.79E+03;


% New formulation Duncan et al., 2018 (in press... soonish)
R = pm.R;
pm.Tsol = 2.17E-07*((R(1)-R)/1e3).^3 - 9.153E-04*((R(1)-R)/1e3).^2 + 1.501E+00*((R(1)-R)/1e3) + 1.365E+03;
pm.Tliq = -1.075E-05*((R(1)-R)/1e3).^2 + 3.730E-01*((R(1)-R)/1e3) + 1.938E+03;

%Temperature for non-magma melt ocean
%Tbase=[2100 +- 200 K]
% Hot is 1100 - 2100 K
% Cold is 1100 - 1700 K

end
