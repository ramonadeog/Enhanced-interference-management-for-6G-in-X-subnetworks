function [Loc,angD,Xway,Yway] = commonWaypointMM(XBound,YBound,sampTime,simTime,Speed,maxAng,nCell,minDist,cellRad)
%% ================ CommonWayPoint Mobility in a Rectangular Grip ===========
%  ===================== Ramoni Adeogun [2020] ============================
%% ========================================================================
%rng('shuffle');
%initial location
XBound1 = XBound-2*cellRad;
YBound1 = YBound-2*cellRad;

%Random start
% Xtemp = cellRad+rand(nCell,1)*XBound;
% Ytemp = cellRad+rand(nCell,1)*YBound;

%%Enforce minimum distance
[gwLoc] = dropPoints(XBound1,YBound1,nCell,minDist);
Xtemp = gwLoc(:,1)+(XBound)/2;
Ytemp = gwLoc(:,2)+(YBound)/2;
Xway = XBound/2; Yway = YBound/2;
N = round(simTime/sampTime);
%D = atand((Yway-Ytemp)./(Xway-Xtemp));
D = atan2((Yway-Ytemp),(Xway-Xtemp));
D = D(:);
S = Speed;

Loc = zeros(2,nCell,N);
angD = zeros(nCell,N);
Loc(:,:,1) = [Xtemp(:) Ytemp(:)]';
angD(:,1) = D;
for n = 2:N
    Xtemp = Xtemp+cos(D).*S*sampTime;
    Ytemp = Ytemp+ sin(D).*S*sampTime;
    distGw = pdist2([Xtemp(:) Ytemp(:)],[Xtemp(:) Ytemp(:)]); 
    distGw = distGw+diag(190*ones(nCell,1));
    [Indx1,~] = find(distGw <= minDist);
    D(unique(Indx1)) = rand(numel(unique(Indx1)),1)*maxAng;
    Xtemp((unique(Indx1))) = Xtemp((unique(Indx1)))+cos(D(unique(Indx1))).*S*sampTime;
    Ytemp((unique(Indx1))) = Ytemp((unique(Indx1)))+ sin(D(unique(Indx1))).*S*sampTime;
     Xtemp(Xtemp < cellRad) = cellRad;
     Xtemp(Xtemp > XBound1) = XBound1;
     Ytemp(Ytemp < cellRad) = cellRad;
     Ytemp(Ytemp > XBound1) = XBound1;
     D(Xtemp==cellRad) = rand(numel(Xtemp(Xtemp==cellRad)),1)*maxAng;
     D(Xtemp==XBound1) = rand(numel(Xtemp(Xtemp==XBound1)),1)*maxAng;
     D(Ytemp==cellRad) = rand(numel(Ytemp(Ytemp==cellRad)),1)*maxAng;
     D(Ytemp==YBound1) = rand(numel(Ytemp(Ytemp==YBound1)),1)*maxAng;
     Loc(:,:,n) = [Xtemp(:) Ytemp(:)]';
     angD(:,n) = D; 
end

end