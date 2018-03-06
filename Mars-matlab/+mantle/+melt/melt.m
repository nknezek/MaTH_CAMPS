function [f,dT]=melt(Ti,pm)

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
    for j=2:length(Ti(:,i))
        if (Ti(j,i) - Tliq(i) >= 0) && (Ti(j,i) - Tliq(i) >= 0) && (Ti(j,i)>Ti(j-1,i)) % find Ti from solver above the liquidus 
            Tt(j,i) = (Ti(j,i)-Ti(j-1,i))+ Tliq(i);                  % temperature value change over the solidus from preceding step - thermal history
            f0(j,i) = (Tt(j,i)-Tsol(i))/(Tliq(i)-Tsol(i));          % melt fraction (not real)
            %Q0(j,i)=(Ti(j,i) - Tliq(i))*Cp(i)*M(i);                 % CURRENTLY UNDERESTIMATES MELT energy with no phase changes (not real), complete melt above liquidus
            Q0(j,i)=(Tt(j,i) - Tliq(i))*Cp(i)*M(i)*f0(j,i);         % energy with no phase changes (not real), complete melt above liquidus
            Q(j,i)=Q0(j,i)-M(i)*L(i);                              % System assumes Q0 is energy max, Q accounts for phase change (IS THIS CORRECT?) 
            if Q(j,i) < 0                                         % if energy is negative, not enough T increase to melt entire layer
                f(j,i)=Q0(j,i)/(M(i)*L(i));                       % Find f that will give Q=0
                Q1(j,i)=Q0(j,i)-M(i)*L(i)*f(j,i);
                if (Q1(j,i)~=0) && (Q1(j,i)/Q(j,i) < 1e-15)        % rounding errors in 16th decimal place common in matlab (what this looks like)
                    Q1(j,i)=0;                                     % sets likely rounding error to 0; Q < ~ 1e27, results in less than 1K change in T
                end
            end
            if Q(j,i) > 0 && countl == 0
%                 f(j,i)=f(j,i)+ Q0(j,i)/(M(i)*L(i));   %cpMdt
%                  display ([num2str(Q(j,i)) ' ' num2str(Tt(j,i)-Tsol(i))])
                  display ('Energy is greater than 0, can no longer neglect temperature change of parcel (e.g. include cp*M*dT + L*M)')
                  display ('Fault is in liquidus approximation: assume either f increases, or melt takes away excess heat')
                  countl=countl+1;
            end
            

            Tt(j,i)=Tliq(i);                                       % NOTE: this is currently incorrect, assumes energy carried away brings T to liquidius/or solidus... think on
           
            if f(j,i) > 1 
                display ('Energy balance requires the melt fraction to exceed 1, violates conditions - melt module is invalid')
                display ('Fault is in liquidus approximation')
                pause
            end
        end


         if (Tt(j,i) <= Tliq(i)) && (Tt(j,i) >= Tsol(i)) % find Ti from solver above the solidus
             Tt(j,i) = (Ti(j,i)-Ti(j-1,i))+ Tsol(i);                 % temperature value change over the solidus from preceding step - thermal history
             if Tt(j,i) > Tsol(i)
                 f0(j,i) = (Tt(j,i)-Tsol(i))/(Tliq(i)-Tsol(i));          % melt fraction(not real)
                 %Q0(j,i)=(Tt(j,i) - Tsol(i))*Cp(i)*M(i)*f0(j,i);         % energy with no phase changes (not real), fractional melt above solidus
                 Q0(j,i)=(Tt(j,i) - Tsol(i))*Cp(i)*M(i);
                 Q(j,i)=Q0(j,i)-M(i)*L(i)*f0(j,i);                              % System assumes Q0 is energy max, Q accounts for phase change (IS THIS CORRECT?)
                 if Q(j,i) < 0                                          % if energy is negative, not enough T increase to melt entire layer
                     f(j,i)=f(j,i)+Q0(j,i)/(M(i)*L(i));                       % Find f that will give Q=0
                     Q1(j,i)=Q0(j,i)-M(i)*L(i)*f(j,i);
                     if (Q1(j,i)~=0) && (Q1(j,i)/Q(j,i) < 1e-15)        % rounding errors in 16th decimal place common in matlab (what this looks like)
                         Q1(j,i)=0;                                     % sets likely rounding error to 0; Q < ~ 1e27, results in less than 1K change in T
                     end
                 end
             end
             
            if Q(j,i) > 0 && counts == 0
                  % follow example above, calcualte a dt, plug into energy,
                  % find a new dt that results in Q=0, this dt of melt will
                  % go into migration module
                  dtm=Q(j,i)/(Cpm*M(i)*f(j,i));                          % excess energy to dt that goes into melt --  delQ=cpm*Mm*dtm of melt      %Q2(j,i)=Q(j,i)-cpm*Mm*dtm
                  Q2(j,i)=Q(j,i)-Cpm*M(i)*f(j,i);
                  if (Q2(j,i)~=0) && (Q2(j,i)/Q(j,i) < 1e-15)        % rounding errors in 16th decimal place common in matlab (what this looks like)
                    Q2(j,i)=0;                                     % sets likely rounding error to 0; Q < ~ 1e27, results in less than 1K change in T
                  end
%                   if Q2 ~= 0
%                       display ([num2str(Q2(j,i)) '  ' num2str(dtm) '  '  'Here there is a problem'])
%                   else 
%                       display ([num2str(Q2(j,i)) '  ' num2str(dtm) '  '   'It works?'])
%                   end
                  %display ([num2str(Q(j,i)) ' ' num2str(Tt(j,i)-Tsol(i))])
                  display ('Energy is greater than 0, can no longer neglect temperature change of parcel (e.g. include cp*M*dT + L*M)')
                  display ('Fault is in solidus approximation: assume either f increases, or melt takes away excess heat')
                  counts=counts+1;
            end
            
            
             Tt(j,i)=Tsol(i);
         end
                                              

    end
end


% fixes jump in temperatures from melting interval to non-melting inteval
% --> sets a "thermal history due to melting
% 
% for i=1:layer;
%     for j=2:length(Ti(:,i));
%         if (Tt(j,i) > Tsol(i)) && (f(j,i) == 0)
%             if count == 0
%                 Toff = (Tt(j,i)- Tsol(i));                                % determine initial offest - Note: will not work if systems warms again from this point
%             end
%             count = count+1;
%             if count > 0 & Toff >1
%                 %Tt(j,i) = ((Tt(j,i)- Tt(j-1,i)) + Tsol(i);           
%                 Tt(j,i) = Tt(j,i) - Toff ;                                % fix remaining Temps to melt modified Temps
%             end
%         end
%     end
% end



%Find change in temperarture due to melting effects  (note: will force Temps to solidus if melt is present)
dT=Ti-Tt;

end
