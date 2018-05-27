function [f,dT]=melt(Tprev,Ti,pm)

%Note: As of now, this is just the phase change energy loss, no energy
%is being expended in heating these parcels up. Will be added functionality
%as needed (at it should be added)

%Note: need to calculate T carried away from layer by melt
% need to account for Tu, Tl, and Ti.

%Note: Tacitly assumes heating up from solidus/sub-solidus, not cooling down from super-solidus/liquidus

% ADD change in T if latent heat does not remove all energy : Q=mcpdT

% Some apparent instability in melt calculations (f may oscillate), check
% these later.
Cp = pm.Cp;
L = pm.L;
M = pm.M;
n = pm.n;
Tsol = pm.Tsol;
Tliq = pm.Tliq;

%% Melting post proc
f0=zeros(size(Ti));       % melt fraction based soley on T (gross over estimation)
Q0=zeros(size(Ti));       % energy avaiable based on melt fraction based f0
Q=zeros(size(Ti));        % energy difference accounting for latent heat of fusion
f=zeros(size(Ti));        % melt fraction based on conservation of energy/phase change (fusion)
Q1=zeros(size(Ti));       % energy balance assuming Q1 =0 and f, a check
Q2=zeros(size(Ti));       % energy balance assuming Q2 =0 and f, a check
fa=zeros(size(Ti));       % melt fraction based soley on T (gross over estimation
Tt=Ti;
count = 0;
countl=0;
counts=0;

layer=n;   % neglect melt in core --> solidii and liquidii are wrong here
for i=1:layer
    if (Ti(i) - Tliq(i) >= 0) && (Ti(i)>Tprev(i)) % find Ti from solver above the liquidus
        Tt(i) = (Ti(i)-Tprev(i)) + Tliq(i);                  % temperature value change over the solidus from preceding step - thermal history
        f0(i) = (Tt(i)-Tsol(i))/(Tliq(i)-Tsol(i));          % melt fraction (not real)
        %Q0(i)=(Ti(i) - Tliq(i))*Cp(i)*M(i);                 % CURRENTLY UNDERESTIMATES MELT energy with no phase changes (not real), complete melt above liquidus
        Q0(i)=(Tt(i) - Tliq(i))*Cp(i)*M(i)*f0(i);         % energy with no phase changes (not real), complete melt above liquidus
        Q(i)=Q0(i)-M(i)*L(i);                              % System assumes Q0 is energy max, Q accounts for phase change (IS THIS CORRECT?)
        if Q(i) < 0                                         % if energy is negative, not enough T increase to melt entire layer
            f(i)=Q0(i)/(M(i)*L(i));                       % Find f that will give Q=0
            Q1(i)=Q0(i)-M(i)*L(i)*f(i);
            if (Q1(i)~=0) && (Q1(i)/Q(i) < 1e-15)        % rounding errors in 16th decimal place common in matlab (what this looks like)
                Q1(i)=0;                                     % sets likely rounding error to 0; Q < ~ 1e27, results in less than 1K change in T
            end
        end
        if Q(i) > 0 && countl == 0
            %                 f(i)=f(i)+ Q0(i)/(M(i)*L(i));   %cpMdt
            %                  display ([num2str(Q(i)) ' ' num2str(Tt(i)-Tsol(i))])
            display ('Energy is greater than 0, can no longer neglect temperature change of parcel (e.g. include cp*M*dT + L*M)')
            display ('Fault is in liquidus approximation: assume either f increases, or melt takes away excess heat')
            countl=countl+1;
        end
        Tt(i)=Tliq(i);                                       % NOTE: this is currently incorrect, assumes energy carried away brings T to liquidius/or solidus... think on
        if f(i) > 1
            display ('Energy balance requires the melt fraction to exceed 1, violates conditions - melt module is invalid')
            display ('Fault is in liquidus approximation')
            pause
        end
    end
    if (Ti(i) <= Tliq(i)) && (Ti(i) >= Tsol(i))             % find Ti from solver above the solidus
        f0(i) = (Ti(i)-Tsol(i))/(Tliq(i)-Tsol(i));          % melt fraction(not real)
        Q0(i) = (Ti(i)-Tsol(i))*Cp(i)*M(i); % energy with no phase changes (not real), fractional melt above solidus
        Q(i) = Q0(i) - M(i)*L(i)*f0(i);      % System assumes Q0 is energy max, Q accounts for phase change (IS THIS CORRECT?)
        if Q(i) < 0                                          % if energy is negative, not enough T increase to melt entire layer
            f(i) = Q0(i)/(M(i)*L(i));                       % Find f that will give Q=0
%             disp(f)
            Q1(i) = Q0(i) - M(i)*L(i)*f(i);
            if (Q1(i)~=0) && (Q1(i)/Q(i) < 1e-15)        % rounding errors in 16th decimal place common in matlab (what this looks like)
                Q1(i)=0;                                     % sets likely rounding error to 0; Q < ~ 1e27, results in less than 1K change in T
            end
        end
        if Q(i) > 0 && counts == 0
            disp(Q)
            % follow example above, calcualte a dt, plug into energy,
            % find a new dt that results in Q=0, this dt of melt will
            % go into migration module
            dtm=Q(i)/(Cpm*M(i)*f(i)); % excess energy to dt that goes into melt --  delQ=cpm*Mm*dtm of melt      
            Q2(i)=Q(i)-Cpm*M(i)*f(i);
            if (Q2(i)~=0) && (Q2(i)/Q(i) < 1e-15)   % rounding errors in 16th decimal place common in matlab (what this looks like)
                Q2(i)=0;                       % sets likely rounding error to 0; Q < ~ 1e27, results in less than 1K change in T
            end

            display ('Energy is greater than 0, can no longer neglect temperature change of parcel (e.g. include cp*M*dT + L*M)')
            display ('Fault is in solidus approximation: assume either f increases, or melt takes away excess heat')
            counts=counts+1;
        end
        Tt(i)=Tsol(i);
    end
end

%Find change in temperarture due to melting effects  (note: will force Temps to solidus if melt is present)
dT=Ti-Tt;

end
