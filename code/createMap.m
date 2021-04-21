function [map] = createMap(xPoints,yPoints,sigma,dc)%%
%clear
N1 = length(xPoints);
N2 = length(yPoints);
G = zeros(N1,N2);
for n=1:N1 % first row of cov. matrix, arranged in a matrix
for m=1:N2
     G(n,m)= 1*exp(-1*sqrt(min(abs(xPoints(1)-xPoints(n)), ...
     max(xPoints)-abs(xPoints(1)-xPoints(n)))^2 + min(abs(yPoints(1)-yPoints(m)), ...
      max(yPoints)-abs(yPoints(1)-yPoints(m)))^2)/dc); %Euclidean distance on the 
   %G(n,m) = sigma*exp(-1*sqrt(abs(xPoints(1)-xPoints(n))^2+abs(yPoints(1)-yPoints(m))^2)/dc);
end
end

Gamma = fft2(G); % the eigenvalue matrix n*fft2(G/n)
Z = randn(N1,N2) + sqrt(-1)*randn(N1,N2);
map = real(fft2(sqrt(Gamma).*Z/sqrt(N1*N2)))*sqrt(2);



end