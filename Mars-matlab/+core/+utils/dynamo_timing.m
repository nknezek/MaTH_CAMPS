function [DTimeh,DTimel] = dynamo_timing(Qu)
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
end