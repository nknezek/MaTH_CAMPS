function [f,dT]=melt_ini(Ti,pm)

%Finds initial melt f and dT to lid from Ti choice; 

%Note: Tacitly assumes heating up from solidus/sub-solidus, not cooling 
%down from super-solidus/liquidus; can add assumption of crystallization 
%from liquid (true magma ocean) later as needed.

%Note: As of now, this is just the phase change energy loss, no energy
%is being expended in heating these parcels up. Will be added functionality
%as needed

% ADD change in T if latent heat does not remove all energy : Q=mcpdT

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

layer=n-1;   % neglect melt in core --> solidii and liquidii are wrong here 
%layer=n;     % Allow for melting in Core, should not 'remove' heat however, retool abit


%Cpm=1500; %Melt Cp from di genova et al 2014


for i=1:layer
    for j=1:1
        if (Ti(j,i) - Tliq(i) >= 0) && (Ti(j,i) - Tliq(i) >= 0)                                % find Ti from solver above the liquidus 
            f0(j,i) = (Ti(j,i)-Tsol(i))/(Tliq(i)-Tsol(i));          % melt fraction (not real)
            %Q0(j,i)=(Ti(j,i) - Tliq(i))*Cp(i)*M(i);                 % CURRENTLY UNDERESTIMATES MELT energy with no phase changes (not real), complete melt above liquidus
            Q0(j,i)=(Ti(j,i) - Tliq(i))*Cp(i)*M(i)*f0(j,i);         % energy with no phase changes (not real), complete melt above liquidus
            Q(j,i)=Q0(j,i)-M(i)*L(i);                              % System assumes Q0 is energy max, Q accounts for phase change (IS THIS CORRECT?) 
            if Q(j,i) < 0                                          % if energy is negative, not enough T increase to melt entire layer
                f(j,i)=Q0(j,i)/(M(i)*L(i));                       % Find f that will give Q=0
                Q1(j,i)=Q0(j,i)-M(i)*L(i)*f(j,i);
                if (Q1(j,i)~=0) && (Q1(j,i)/Q(j,i) < 1e-15)        % rounding errors in 16th decimal place common in matlab (what this looks like)
                    Q1(j,i)=0;                                     % sets likely rounding error to 0; Q < ~ 1e27, results in less than 1K change in T
                end
            end
            if Q(j,i) > 0
%                 f(j,i)=f(j,i)+ Q0(j,i)/(M(i)*L(i));   %cpMdt
%                  display ([num2str(Q(j,i)) ' ' num2str(Tt(j,i)-Tsol(i))])
                  display ('Energy is greater than 0, can no longer neglect temperature change of parcel (e.g. include cp*M*dT + L*M)')
            end
            
            % this approach causes see-saw problems, think on.
%             if (j < length(Ti(:,i))) && (Tt(j+1,i) > Tliq(i)) 
%                 Tt(j+1,i)=(Ti(j,i) - Ti(j+1,i)) +Tliq(i);
%                 
%             end

            Tt(j,i)=Tliq(i);                                       % NOTE: this is currently incorrect, assumes energy carried away brings T to liquidius/or solidus... think on
           
            if f(j,i) > 1
                display ('Energy balance requires the melt fraction to exceed 1, violates conditions - melt module is invalid')
                display ('Fault is in liquidus approximation')
                pause
            end
        end


         if (Tt(j,i) <= Tliq(i)) && (Tt(j,i) >= Tsol(i)) % find Ti from solver above the solidus
             f0(j,i) = (Tt(j,i)-Tsol(i))/(Tliq(i)-Tsol(i));          % melt fraction(not real)
             Q0(j,i)=(Tt(j,i) - Tsol(i))*Cp(i)*M(i)*f0(j,i);         % energy with no phase changes (not real), fractional melt above solidus
             Q(j,i)=Q0(j,i)-M(i)*L(i)*f0(j,i);                              % System assumes Q0 is energy max, Q accounts for phase change (IS THIS CORRECT?)
            if Q(j,i) < 0                                          % if energy is negative, not enough T increase to melt entire layer
                f(j,i)=f(j,i)+ Q0(j,i)/(M(i)*L(i));                       % Find f that will give Q=0
                Q1(j,i)=Q0(j,i)-M(i)*L(i)*f(j,i);
                if (Q1(j,i)~=0) && (Q1(j,i)/Q(j,i) < 1e-15)        % rounding errors in 16th decimal place common in matlab (what this looks like)
                    Q1(j,i)=0;                                     % sets likely rounding error to 0; Q < ~ 1e27, results in less than 1K change in T
                end
            end
            if Q(j,i) > 0
%                 f(j,i)=f(j,i)+ Q0(j,i)/(M(i)*L(i));   %cpMdt
%                  display ([num2str(Q(j,i)) ' ' num2str(Tt(j,i)-Tsol(i))])
                  display ('Energy is greater than 0, can no longer neglect temperature change of parcel (e.g. include cp*M*dT + L*M)')
            end
            
            
             Tt(j,i)=Tsol(i);
         end
                                              

    end
 
end
    %Find change in temperarture due to melting effects  (note: will force Temps to solidus if melt is present)
dT=Ti-Tt;
end
