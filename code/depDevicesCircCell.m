function [X]=depDevicesCircCell(centerLoc,rc,numDevices)
[L, ~] = size(centerLoc);
X = zeros(2,L,numDevices);
for c = 1:L
    a=2*pi*rand(numDevices,1);
    r=sqrt(rand(numDevices,1));
    X(1,c,:)=(rc*r).*cos(a)+centerLoc(c,1);
    X(2,c,:)=(rc*r).*sin(a)+centerLoc(c,2);
end
end
