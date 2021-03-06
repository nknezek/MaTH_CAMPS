%% Test Thermal Evolution Code 
% Adiabatic convection up and thermal conduction down
% Nicholas Knezek

N = 100;

R = 2000e3; % [m]
dr = R/N; % [m]
dt = 365.25*24*3600*1e6; % s (1 Myr)
D = 6203e3; % m, adiabatic lengthscale

Nt = 4500;

cp = 840; %J/kg-K = m^2/s^2-K
rho = linspace(12e4, 9e4, N); % kg/m^3
k = 135; % J/m-K-s = kg-m/K-s^3

re = linspace(0,R,N+1); % m
r = linspace(dr/2, R-dr/2, N); %m
dV = 4*pi*dr*r.^2; % m^3
A = 4*pi*re.^2; % m^2

rho_cp_dV = rho.*cp.*dV; % J/K

T_cmb = 2300; % K
T_ad = adiabat(T_cmb, r, R, D); % K
T_p = linspace(0,0,N); % K
T_0 = T_ad;

T_vec0 = horzcat(T_p,T_cmb);

Q_cmb_t = ones(1,Nt);
% Q_cmb_t(1:end) = 2e12;
Q_cmb_t(1:500) = linspace(8e12, -0.25e12,500);
Q_cmb_t(501:1500) = linspace(-0.25e12,0.75e12,1000);
Q_cmb_t(1501:end) = linspace(0.75e12,0.25e12,Nt-1500);

E_ad = zeros(1,Nt);
E_p = zeros(1,Nt);
T_vec_all = zeros(Nt, N+1);

%%
% Q_cmb = 5e12;
% Nt = 1500;
T = T_ad;
T_cmb = zeros(1,Nt);
T_cen = zeros(1,Nt);
E_core = zeros(1,Nt);
Q_ad_c_out = zeros(1,Nt);
dR_TBL = zeros(1,Nt);
% plot(r,T,'.-')
% hold on
for it=1:Nt
    E_i = energy(T,rho_cp_dV); % Initial Energy of Core
    
    % Conduct heat down gradients
    qe = -k*(T(2:end)-T(1:end-1))./dr; % W/m^2
    Qe = zeros(1, length(T)); % W, heat flow across layer boundaries 
    Qe(1:end-1) = qe.*A(2:end-1); % W
    Qe(end) = Q_cmb_t(it); % W , heat flow out top of core 
    Q = zeros(1,length(T)); % W
    Q(2:end) = Qe(2:end)-Qe(1:end-1); % W
    Q(1) = Qe(1); % W, heat flow from central point only goes up
    dT_cond = -Q.*dt ./rho_cp_dV; % K
    T_cond = T+dT_cond; % K

    % Convective Step
    conv = zeros(1,N);
    cond = zeros(1,N);
    Ch_ad_total = 0;
    T_out = zeros(1,N);
    in_an_adiabat = 0;
    T_regions = T_cond;
    for i=N:-1:2
        T_im1a = adiabat_from(T_regions(i),r(i),r(i-1),D);
        if T_im1a <= T_regions(i-1)*1.0001 % in convective regime
            conv(i) = 1;
            if in_an_adiabat == 0 % in top layer of convective regime
                i_ad_top = i;
                Q_ad_top = Qe(i);
            end
            in_an_adiabat = 1;
            Ch_ad_total = Ch_ad_total + rho_cp_dV(i)*exp((r(i_ad_top)^2 - r(i)^2)/D^2);
        elseif T_im1a > T_regions(i-1)*1.0001 % in conductive regime
            cond(i) = 1;
            if in_an_adiabat == 1% in top layer of conductive regime
                % deal with adiabat that we just left
                Ch_ad_total = Ch_ad_total + rho_cp_dV(i)*exp((r(i_ad_top)^2 - r(i)^2)/D^2);
                T_ad_top_new = T(i_ad_top)-Q_ad_top*dt/Ch_ad_total;
                for j=i:i_ad_top
                    T_out(j) = adiabat_from(T_ad_top_new, r(i_ad_top), r(j), D);
                end
%                 if i_ad_top ~= N % if top of adiabat is not the CMB, conduct heat into top of convective layer to grow TBL
%                     T_out(i_ad_top) = T_out(i_ad_top) + dT_cond(i_ad_top);
%                 end

                in_an_adiabat = 0;
                Ch_ad_total = 0;
            else % in at least second layer of conductive region
                T_out(i) = T_cond(i);
            end
        end
    end
    if in_an_adiabat == 1
        i = 1;
        conv(i) = 1;
        Ch_ad_total = Ch_ad_total + rho_cp_dV(i)*exp((r(i_ad_top)^2 - r(i)^2)/D^2);
        T_ad_top_new = T(i_ad_top)-Q_ad_top*dt/Ch_ad_total;
        for j=i:i_ad_top
            T_out(j) = adiabat_from(T_ad_top_new, r(i_ad_top), r(j), D);
        end
%         if i_ad_top ~= N % if top of adiabat is not the CMB, conduct heat into top of convective layer to grow TBL
%             T_out(i_ad_top) = T_out(i_ad_top) - dT_cond(i_ad_top);
%         end
    end
    T = T_out;
    % Correction to make sure you satisfy energy conservation
    E_o = energy(T,rho_cp_dV); % final energy of core, J
    dEq = Qe(end)*dt; % Heat pulled out of core, J
    dEt = E_i-E_o; % change in internal energy, J
    dE_corr = dEt-dEq;
    dT_corr = dE_corr/sum(rho_cp_dV.*exp(-r.^2./D^2));
    T = T + dT_corr*exp(-r.^2./D^2);

    % Storing things for output
    T_cmb(it) = T(end);
    T_cen(it) = T(1);
    Q_ad_c_out(it) = Q_ad_cmb(T(end),k,R,D);
    E_core(it) = energy(T,rho_cp_dV);
    if sum(cond)>0
        dR_TBL(it) = (R-min(r(cond==1)))/1e3;
    else
        dR_TBL(it) = 0;
    end
    
    if mod(it,Nt/10)==1
        hold on
        plot(r,T)
        grid on
%         disp(it)
    end
end
%%
plot(r,T_ad)
hold on
plot(r,T)

figure()
plot(r,cond)
hold on
plot(r,conv)
legend('cond','conv')

disp((R-min(r(cond==1)))/1e3)

%%
E_i = energy(T_ad,rho_cp_dV);
E_o = energy(T,rho_cp_dV);
dEt = E_o-E_i
dEq = -sum(Q_cmb_t(1:Nt)*dt)
%%
figure()
plot(Q_cmb_t(1:Nt))
hold on
plot(Q_ad_c_out(1:Nt))

figure()
plot(T_cmb)
hold on
plot(T_cen)

figure()
plot(E_core)

figure()
plot(dR_TBL)
%% find convective and conductive regions 
% T = T_cond;
% conv = zeros(1,N);
% cond = zeros(1,N);
% for i=N:-1:2
%     T_im1a = adiabat_from(T(i),r(i),r(i-1),D);
%     if T_im1a <= T(i-1)*1.001
%         conv(i) = 1;
%     elseif T_im1a > T(i-1)*1.001
%         cond(i) = 1;
%     end
% end
% plot(r,conv)
% hold on
% plot(r,cond)
% figure()
% plot(r,T)
%% Convection Step in convective regions
conv = zeros(1,N);
cond = zeros(1,N);
Ch_ad_total = 0;
T_out = zeros(1,N);
in_an_adiabat = 0;
for i=N:-1:2
    T_im1a = adiabat_from(T(i),r(i),r(i-1),D);
    if T_im1a <= T(i-1)*1.0001 % in convective regime
        conv(i) = 1;
        if in_an_adiabat == 0 % in top layer of convective regime
            i_ad_top = i;
            Q_ad_top = Qe(i);
        end
        in_an_adiabat = 1;
        Ch_ad_total = Ch_ad_total + rho_cp_dV(i_ad_top)*exp((r(i_ad_top)^2 - r(i)^2)/D^2);
    elseif T_im1a > T(i-1)*1.0001 % in conductive regime
        cond(i) = 1;
        if in_an_adiabat == 1% in top layer of conductive regime or at center
            % deal with adiabat that we just left
            Ch_ad_total = Ch_ad_total + rho_cp_dV(i_ad_top)*exp((r(i_ad_top)^2 - r(i)^2)/D^2);
            T_ad_top_new = T(i_ad_top)-Q_ad_top*dt/Ch_ad_total;
            disp(i_ad_top)
            disp(T(i_ad_top))
            disp(T_ad_top_new)
            for j=i:i_ad_top
                T_out(j) = adiabat_from(T_ad_top_new, r(i_ad_top), r(j), D);
            end
        in_an_adiabat = 0;
        Ch_ad_total = 0;
        else % in at least second layer of conductive region
            T_out(i) = T_cond(i);
        end
    end
end
if in_an_adiabat == 1
    disp(i_ad_top)
    disp(T(i_ad_top))
    disp(T_ad_top_new)
    i = 1;
    Ch_ad_total = Ch_ad_total + rho_cp_dV(i_ad_top)*exp((r(i_ad_top)^2 - r(i)^2)/D^2);
    T_ad_top_new = T(i_ad_top)-Q_ad_top*dt/Ch_ad_total;
    for j=i:i_ad_top
        T_out(j) = adiabat_from(T_ad_top_new, r(i_ad_top), r(j), D);
    end
end
plot(r,conv)
hold on
plot(r,cond)
figure()
plot(r,T_cond ,'.-')
hold on
plot(r,T_out,'.-')

legend('before','after')





%%
% for i=N:-1:1
%     T_cond(i)
% end

% T_vec = T_vec0;
% for i=1:Nt
%     T_vec = evolve_core(T_vec, Q_cmb_t(i));
%     T_vec_all(i,:) = T_vec;
%     T_p = T_vec(1:end-1);
%     T_ad = adiabat(T_vec(end),r,R,D);
%     E_ad(i) = energy(T_ad,rho_cp_dV);
%     E_p(i) = energy(T_p, rho_cp_dV);
%     if mod(i,100)==0
%         hold on
%         plot(r,T_ad+T_p,'.-')
%     end
% end
%%
t = linspace(0,Nt,Nt);
% plot(t,E_ad)
plot(t,T_vec_all(:,end))
% plot(t,E_p)
% plot(t,E_ad+E_p)
% plot(t,E_ad+E_p+Q_cmb_t*dt)


%%
T_vec = T_vec0;
Q_cmb = 3.5e12;
odef = @(t,y) dTdt_core(y, Q_cmb);
[t, T_v] = ode45(odef, [0, Nt*dt], T_vec0);
T_cmb = T_v(:,end)+T_v(:,end-1); 
plot(t,T_cmb)