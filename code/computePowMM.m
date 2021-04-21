function [rxPow] = computePowMM(S,dist,sigmaVec,txPow,gammaP,freq,BW,B,N,M)
d0 = 1;
fc = (freq+BW/2)+(0:B-1)*BW;
[S1,S2] = size(S);
S = S*sigmaVec(2);
gamma = ones(S1,S2)*gammaP(2);
for n = 1:N
    S(n,N+(n-1)*M+1:N+n*M) = S(n,N+(n-1)*M+1:N+n*M)*sigmaVec(1)/sigmaVec(2);
    S(N+(n-1)*M+1:N+n*M,n) = S(N+(n-1)*M+1:N+n*M,n)*sigmaVec(1)/sigmaVec(2);
    gamma(n,N+(n-1)*M+1:N+n*M) = gamma(n,N+(n-1)*M+1:N+n*M)*gammaP(1)/gammaP(2);   
    gamma(N+(n-1)*M+1:N+n*M,n) = gamma(N+(n-1)*M+1:N+n*M,n)*gammaP(1)/gammaP(2);
end
rxPow = zeros(S1,S2,B);
for b = 1:B
    rxPow(:,:,b) = (txPow- (20*log10(4*pi*fc(b)*d0/3e8)...
        +10*gamma.*log10(dist/d0)+S));
end 

rxPow(isinf(rxPow)) = 0;
end

    