function [C,G] = centralizedColoring(interMat,numChannel,option)
%% ====================================================================
% ======================= Centralized Graph Coloring====================
%% ======================================================================
%distance
%Create neigbouhood matrix:

N = size(interMat,1);
S = []; T = [];
for ii = 1:N
    [~,Indx] = sort(interMat(ii,:));
     S = [S ii*ones(1,numChannel-1)];
     T = [T Indx(2:numChannel)];
end

GG = graph(S,T);
G = simplify(GG);  
CC = color(G,option);   %See other options in color
[C,~] = find(full(CC.array)==1);
while max(C) > numChannel
    [~,Indxmin] = min(min(interMat));
    Indxx = find(S==Indxmin);
    [~,Indx] = sort(interMat(Indxmin,:));
    interMat(Indxmin,Indx) = 100;
    S(Indxx(end)) = [];
    T(Indxx(end)) = [];
    GG = graph(S,T);
    G = simplify(GG);  
    CC = color(G,option);   %See other options in color
    [C,~] = find(full(CC.array)==1);
end

% % [C1,~] = find(full(CC1.array)==1);

%Uncomment if plotting of coloring is desired 
% if plotoption == 1
%     figure(1);
%     H = plot(G,'XData',Loc(:,1),'YData',Loc(:,2));
%     H.NodeColor = C2(C,:);
%     grid on; 
%     xlabel('Length [m]');
%     ylabel('Width [m]');
%     figure(2);
%     H = plot(G,'XData',Loc(:,1),'YData',Loc(:,2));
%     H.NodeColor = C2(C1,:);
%     grid on; 
%     xlabel('Length [m]');
%     ylabel('Width [m]');
%     %title(['Maximum MC size:' num2str(maxClique)]);
% end

C = C';


end