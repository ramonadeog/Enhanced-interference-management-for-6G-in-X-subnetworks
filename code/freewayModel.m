function [Loc,angD] = freewayModel(XBound,YBound,sampTime,simTime,Speed,nCell,cellRad)
%Freeway Mobility model 
XBound1 = XBound-2*cellRad;
YBound1 = YBound-2*cellRad;
minDist = cellRad;
% cellRad = cellRad;
num_lanes = 6;
num_per_lane = nCell/num_lanes;
lane_width = (YBound)/num_lanes;

lane_y = (0:num_lanes-1)*lane_width + lane_width/2; 
lane_y(4:6) = lane_y(4:6); %1m separate between roads
%starting x location

%gwLoc_x = XBound1*rand(nCell,1);
%gwLoc_x =  [2 4 6 4 2 6 3 5 7 5 9 7 4 2 6 7 9 5]'+4; %+randi(3,18,1);
%gwLoc_x =  [15 10 16 48 18 25 12 16 20 32 11 22 19 13 5 7 15 22 14 54 34 42 21 52]';
% gwLoc_x = rand(nCell,1)*(XBound-4*cellRad) + 2*cellRad;
% gwLoc_x(1:nCell/2,1) = randmin(nCell/2,2*cellRad/XBound)*(XBound/2)+2*cellRad;
% gwLoc_x(nCell/2+1:nCell,1) = randmin(nCell/2,2*cellRad/XBound)*(XBound/2)+(XBound/2)-2*cellRad;
gwLoc_x1 = randmin(nCell/3,3.5/60)*(XBound-4*cellRad) + 2*cellRad;
gwLoc_x2 = randmin(nCell/3,3.5/60)*(XBound-4*cellRad) + 2*cellRad;
gwLoc_x3 = randmin(nCell/3,3.5/60)*(XBound-4*cellRad) + 2*cellRad;
gwLoc_x = [gwLoc_x1(:); gwLoc_x2(:); gwLoc_x3(:)];
gwLoc_y = reshape(repmat(lane_y,num_per_lane,1),[],1);
lane_index = reshape(repmat(1:num_lanes,num_per_lane,1),[],1);
%gwLoc = [gwLoc_x(:) gwLoc_y(:)]
S = (Speed(2)-Speed(1))*ones(nCell,1)+Speed(1); 
%%Enforce minimum distance
N = round(simTime/sampTime);
D = zeros(nCell,1);
D(lane_index <= num_lanes/2) = pi;
D(lane_index > num_lanes/2) = 2*pi;

Loc = zeros(2,nCell,N);
angD = zeros(nCell,N);
Loc(:,:,1) = [gwLoc_x(:) gwLoc_y(:)]';
angD(:,1) = D;
Xtemp = gwLoc_x(:);
Ytemp = gwLoc_y;
for n = 2:N
      
    D((lane_index <= 3 & Xtemp <= 2*cellRad)) = 2*pi;
    D((lane_index > 3 & Xtemp >= XBound-2*cellRad)) = pi; 
    new_index = randi(3,length(Ytemp((lane_index <= 3 & Xtemp <= 2*cellRad))),1)+3;
    Ytemp((lane_index <= 3 & Xtemp <= 2*cellRad)) = lane_y(new_index);
    lane_index((lane_index <= 3 & Xtemp <= 2*cellRad)) = new_index;
    new_index1 = randi(3,length(Ytemp((lane_index > 3 & Xtemp >= XBound-2*cellRad))),1);
    Ytemp((lane_index > 3 & Xtemp >= XBound-2*cellRad)) = lane_y(new_index1);
    lane_index(lane_index > 3 & Xtemp >= XBound-2*cellRad) = new_index1;
    Indexx = (lane_index <= 3 & Xtemp <= 2*cellRad);
    Xtemp((lane_index <= 3 & Xtemp <= 2*cellRad)) = 2*cellRad*rand(length(find(Indexx ==1)),1);
    Indexx1 = (lane_index > 3 & Xtemp >= XBound-cellRad);
    Xtemp((lane_index > 3 & Xtemp >= XBound-cellRad)) = XBound-2*cellRad*rand(length(find(Indexx1==1)),1);
    Xtemp = Xtemp+cos(D).*S*sampTime;  
  
    Loc(:,:,n) = [Xtemp(:) Ytemp(:)]';
    angD(:,n) = D;
    % Minimum distance constraint
%     distGw = pdist2([Xtemp(:) Ytemp(:)],[Xtemp(:) Ytemp(:)]); 
%     distGw = distGw+diag(190*ones(nCell,1));
%     [Indx1,~] = find(distGw <= minDist);
%     D(unique(Indx1)) = randi(2,size(unique(Indx1)))*pi;
%     size(unique(Indx1))
%     Xtemp((unique(Indx1))) = Xtemp((unique(Indx1)))+cos(D(unique(Indx1))).*S*sampTime;
%     Ytemp((unique(Indx1))) = Ytemp((unique(Indx1)))+ sin(D(unique(Indx1))).*S*sampTime;
    
%     distGw = pdist2([Xtemp(:) Ytemp(:)],[Xtemp(:) Ytemp(:)]); 
%     distGw = distGw+diag(190*ones(nCell,1));
%     [Indx1,~] = find(distGw <= 3.5);
%     S(unique(Indx1)) = (Speed(2)-Speed(1))*rand(length(unique(Indx1)),1)+Speed(1);
end

end