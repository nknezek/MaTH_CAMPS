%% Parameters for Martian Core
global pc
pc = struct;

%% Number of layers
pc.N = 150; % [-]

%% Radius of Core
pc.r_cmb = 2000e3; % [m] (guess)

%% Pressure at CMB
pc.P_cmb = 20e9; % [Pa] (guess)

%% Density at center of Martian Core
pc.rho_cen = 7e3; % [kg/m^3] (guess)

%% Thermal expansivity [/K]
% alpha = 1e-5; %[/K] thermal expansivity (Buffett 1996)
% alpha = 0.89e-5; %[/K] thermal expansivity (Roberts et al. 2003 low)
% alpha = 1.77e-5; %[/K] thermal expansivity (Roberts et al. 2003 high)
% alpha = 1.25e-5; %[/K] thermal expansivity (Labrosse 2003) low
% alpha = 1.7e-5; %[/K] thermal expansivity (Labrosse 2003) high
% alpha = 1.02e-5; %[/K] thermal expansivity (Pozzo et al. 2012) low
% alpha = 1.95e-5; %[/K] thermal expansivity (Pozzo et al. 2012) high
alpha = 1.25e-5; %[/K] thermal expansivity (Nimmo 2015)
pc.alpha = alpha;

%% heat capacity [J/kg-K]
% cp = 800; % (Buffett 1996)
% cp = 819; % (Roberts et al. 2003 low)
% cp = 850; % (Roberts et al. 2003 high)
% cp = 850; % (Labrosse 2003)
% cp = 715; % (Pozzo et al. 2012)
cp = 840; % (Nimmo 2015)
pc.cp = cp;

%% Thermal Conductivity [W/m-K]
% k = 35; % (Buffett 1996)
% k = 46; % (Roberts et al. 2003)
% k = 50; % (Labrosse 2003)
% k = 100; % (Pozzo et al. 2012)
k = 130; % (Nimmo 2015)
pc.k = k;

%% Compostional Expansivity [-]
% alpha_c = 0.93; % (Buffett 1996)
% alpha_c = 1.0; % (Roberts et al. 2003)
% alpha_c = 1.1; % (Pozzo et al. 2012)
alpha_c = 1.1; % (Nimmo 2015)
pc.alpha_c = alpha_c;

%% Latent Heat of Solification of Earth's Inner Core [J/kg]
% L_H = 600e3; % (Buffett 1996)
% L_H = 1560e3; % (Roberts et al. 2003)
% L_H = 660e3; % (Labrosse 2003)
% L_H = 750e3; % (Pozzo et al. 2012)
L_H = 750e3; % (Nimmo 2015)
pc.L_H = L_H;

%% Compressibility of iron at 0 Pa
K_0 = 500e9; % [Pa] (Nimmo 2015)
pc.K_0 = K_0;

%% Density of iron at 0 Pa
rho_0 = 7900; % [kg/m^3] (Nimmo 2015)
pc.rho_0 = rho_0;

%% Constants of Nature
pc.G = 6.6743e-11; % [m^3/kg-s^2] - Gravitational Constant
pc.Myr = 365.25*24*3600*1e6; % [s] - 1 Myr

%% Adiabatic Length Scale
pc.D = core.utils.D_scale(pc.rho_cen); % [m]

%% Pressure Length Scale
pc.L = core.utils.L_scale(pc.rho_cen); % [m]

%% set up radial location vectors
pc.dr = pc.r_cmb/pc.N; % [m]
pc.r = linspace(pc.dr/2,pc.r_cmb-pc.dr/2,pc.N); % [m] (1,N)
pc.re = linspace(0,pc.r_cmb,pc.N+1); % [m] (1,N+1)
pc.A = 4*pi*pc.re.^2; % [m^2] (1,N+1)
pc.dV = 4*pi*pc.r.^2.*pc.dr; % [m^3] (1,N)

%% Compute density through core
pc.rho = core.utils.density(pc.r, pc.rho_cen, pc.D); % [m^3] (1,N)

%% Compute Pressure through core (Pa)
pc.P = core.utils.pressure(pc.r);

%% Calculate helper parameters
pc.rho_cp_dV = pc.rho.*pc.cp.*pc.dV; % J/K (1,N)

