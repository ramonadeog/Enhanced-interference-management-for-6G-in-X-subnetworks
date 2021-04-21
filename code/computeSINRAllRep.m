function [upSINR,downSINR,intPow,IndxS] = computeSINRAllRep(WIRT,action,t,IndxI,option)
%COMPUTESINR Summary of this function goes here
%Detailed explanation goes here

%%%
[B,L] = size(WIRT.envConstant.frameS);
switch option
    case 'random'
        intPow = 0;
        N = WIRT.envConstant.numCell;
        M = WIRT.envConstant.numDev;
        P = WIRT.envConstant.numRep;
        F = WIRT.envConstant.frameS;
        IndxC = zeros(N,M,P);
        frameC = zeros(N,B,L);
        for n = 1:N
            rng(seed);
            frN = F(randperm(size(F,1)),:);
            for m = 1:M
                IndxC(n,m,:) = find(frN==m);
            end
            frN = ismember(frN,1:M); 
            frameC(n,:,:) = frN;
        end
        
          %%%
        rxPow = WIRT.rxPower;
        for n = 1:N
            Indxn = squeeze(IndxC(n,:,:));
            pown = repmat(db2pow(squeeze(rxPow(n,N+(n-1)*M+1:N+(n)*M,t))),P,1);
            upIn = zeros(M,P);
            downIn = zeros(M,P);
            for m = 1:N
                if n ==m
                    
                else
                    %Indxm = squeeze(IndxC(m,:,:));
                    frameCm = squeeze(frameC(m,:,:));
                 
                    upI = db2pow(squeeze(rxPow(n,N+(m-1)*M+1:N+(m)*M,t)));
                    downI = db2pow(squeeze(rxPow(N+(n-1)*M+1:N+(n)*M,m,t)));
                    for p = 1:P
                        Indx = frameCm(squeeze(Indxn(:,p)));
                        upIn(:,p) = upIn(:,p)+upI'.*Indx;
                        downIn(:,p) = downIn(:,p)+downI.*Indx;
                    end
                end
              
            end
            upSINR(n,:,:) = pown'./(upIn+WIRT.envConstant.noisePow);
            downSINR(n,:,:) = pown'./(downIn+WIRT.envConstant.noisePow);
            %SINR(n) = min([min(upSINR(n,:)) min(downSINR(n,:))]);
        end
         
    otherwise 
        
      
        
        frameSub = WIRT.envConstant.frameS(1:B/WIRT.envConstant.numGroups,:); 

        %Create randomized frame for all cells
        frameC = zeros(WIRT.envConstant.numCell,B/WIRT.envConstant.numGroups,L);
        for n = 1:WIRT.envConstant.numCell
            rng(1);
            frN = frameSub(randperm(size(frameSub,1)),:);
            [~,frN] = ismember(frN,[1:WIRT.envConstant.numDev]);
            frameC(n,:,:) = frN;
        end


        %%%
        rxPow = WIRT.rxPower;
        N = WIRT.envConstant.numCell;
        M = WIRT.envConstant.numDev;
        cIndx = action;
        intPow = zeros(N,WIRT.envConstant.numGroups);
        domRatio = zeros(N,WIRT.envConstant.numGroups);
        P = WIRT.envConstant.numRep;
        for n = 1:N
            framen = squeeze(frameC(n,:,:));
            pown = repmat(db2pow(squeeze(rxPow(n,N+(n-1)*M+1:N+(n)*M,t))),P,1);
            upIn = zeros(M,P);
            downIn = zeros(M,P);
            for m = 1:N
                framem = squeeze(frameC(m,:,:));
                if n~=m && cIndx(n)==cIndx(m)
                    upI = db2pow(squeeze(rxPow(n,N+(m-1)*M+1:N+(m)*M,t)));
                    downI = db2pow(squeeze(rxPow(N+(n-1)*M+1:N+(n)*M,m,t)));
                    for p = 1:P
                        framenn = squeeze(framen(p,:));
                        framemm = squeeze(framem(p,:));
                        [~,Indx] = ismember(framenn(framenn~=0),framenn(framemm~=0));
                        if isempty(Indx)
                           upIn(:,p) = upIn(:,p)+0;
                           downIn(:,p) = downIn(:,p)+0;
                            
                        elseif p < P
                            upIn(:,p) = upIn(:,p)+upI(Indx)';
                            downIn(:,p) = downIn(:,p)+downI(Indx);
                        else 
                           
                            IndxII = (~ismember(Indx,IndxI{m}));
                            upIn(IndxII,p) = upIn(IndxII,p)+upI(~ismember(Indx,IndxI{m}))';
                            downIn(IndxII,p) = downIn(IndxII,p)+downI(~ismember(Indx,IndxI{m}));
                        end
                    end
                else

                end
            end
            upSINR(n,:,:) = pown'./(upIn+WIRT.envConstant.noisePow);
            downSINR(n,:,:) = pown'./(downIn+WIRT.envConstant.noisePow);
            upSINR(n,IndxI{n},P) = 1e-2;
            downSINR(n,IndxI{n},P) = 1e-2;
            IndxS{n} = find((squeeze(upSINR(n,:,1))) >= WIRT.envConstant.SINRvec(3));
            %SINR(n) = min([min(upSINR(n,:)) min(downSINR(n,:))]);
        end

        for m = 1:WIRT.envConstant.numGroups
            [indGwm1,~] = find(action'==m);
            for n = 1:N
                rxP = db2pow(squeeze(rxPow(n,indGwm1(indGwm1~=n),t)));
                intPow(n,m) = sum(rxP);
                if isempty(rxP)
                    rxP = 0;
                end
                domRatio(n,m) = max(rxP)/sum(rxP(rxP~= max(rxP)));
            end
        end
        intPow = pow2db(intPow+WIRT.envConstant.noisePow);
        %intPow(intPow == -inf) = -100;
end
        

end

