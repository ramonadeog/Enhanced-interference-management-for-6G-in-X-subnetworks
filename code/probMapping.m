function [PP] = probMapping(Pmap,Smap,SINR)
%Outage probability mapping over multiple channels
%Inputs:
%SINR=downSINR;

%Ramoni Adeogun (June 2019)
[N,M,P] = size(SINR);
SINR = round(SINR*10)/10; %round to 1 decimal place
Smap = round(Smap*10)/10;
PP = zeros(N,M);
for n = 1:N
    for m = 1:M
        Pnm = 1;
        for p = 1:P
            S1 = SINR(n,m,p);
            if S1 < min(Smap)
                Indx = 1;
                Pnm = Pnm*Pmap(Indx);
            elseif S1 > max(Smap)
                Indx = length(Smap);
                Pnm = Pnm*Pmap(Indx);
                 %Pnm = Pnm*0;
            else
              Indx = find(Smap == S1);
                if isempty(Indx)
                    Indx = find(round(Smap)==round(S1));
                end
                %% 
                Pnm = Pnm*Pmap(Indx(1));
%                 if length(Indx) > 1
%                     Indx = Indx(1);
%                     Pnm = Pnm*Pmap(Indx);
%                 else
%                    Pnm = Pnm*Pmap(Indx);
%                 end
            end
           % [Indx Pmap(Indx)]
            %Pnm = Pnm*Pmap(Indx);
          %pause  
        end
        %pause
        PP(n,m) = Pnm;
    end
end
end