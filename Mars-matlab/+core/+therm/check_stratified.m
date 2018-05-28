function stratified = check_stratified(T, pc)
% Computes one timestep convecting and conducting temperature in core
% 
% Parameters
% T [K] - (1,N) vector of temperatures in each layer

r = pc.r;
N = pc.N;
dt = 365.25*24*60*60*1e6*0.01;

% Find regions in which convection occurs
stratified = zeros(N,1); % is the current layer stratified
T_regions = T; % temperatures to use to define convective/conductive regions
Q_typical = 1e11;
epsilonTp1 = 1+min(Q_typical*dt/mean(pc.rho_cp_dV)*1e-4,1e-4);
for i=N:-1:2
    T_im1a = core.utils.adiabat_from(T_regions(i),r(i),r(i-1), pc);
    if T_im1a <= T_regions(i-1)*epsilonTp1 % in convective regime
        stratified(i) = 0;
    elseif T_im1a > T_regions(i-1)*epsilonTp1 % in conductive regime
        stratified(i) = 1;
    end
end
end