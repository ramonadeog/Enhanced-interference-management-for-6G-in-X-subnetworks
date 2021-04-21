function [upSINR,downSINR,intPow] = computeSINR(WIRT,t)
%=========================================================================
frameS = WIRT.envConstant.frameS;   %maximum number of devices
N = WIRT.envConstant.numCell;
M = WIRT.envConstant.numDev;
P = WIRT.envConstant.numRep;
maxDev = max(max(frameS));    %maximum number of devices
if M > maxDev
    error('Maximum number of devices exceeded!!')
end
load_fract = M/maxDev; 
rxPow = WIRT.rxPower;
for n = 1:N
    pown = repmat(db2pow(squeeze(rxPow(n,N+(n-1)*M+1:N+(n)*M,t))),P,1);
    upIn = zeros(M,P);
    downIn = zeros(M,P);
    for m = 1:N
        if n ~= m
            upI = db2pow(squeeze(rxPow(n,N+(m-1)*M+1:N+(m)*M,t)));
            downI = db2pow(squeeze(rxPow(N+(n-1)*M+1:N+(n)*M,m,t)));
            for p = 1:P
                Indx = rand(M,1) <= load_fract; 
                upIn(:,p) = upIn(:,p)+upI'.*Indx;
                downIn(:,p) = downIn(:,p)+downI.*Indx;
            end
        end

    end
    upSINR(n,:,:) = pown'./(upIn+WIRT.envConstant.noisePow);
    downSINR(n,:,:) = pown'./(downIn+WIRT.envConstant.noisePow);
end

for m = 1:WIRT.envConstant.numGroups
    [indGwm1,~] = find(action'==m);
    for n = 1:N
        intPow(n,m) = sum(db2pow(rxPow(n,indGwm1(indGwm1~=n),t)));
    end
end
intPow = pow2db(intPow+WIRT.envConstant.noisePow);

end

