function [S,dist] = computeShadowing(map,X,Y,dc,Loc)
%%Description=======================
%%Calculate the shadowing from Locations and map
[~,idx] = ismember(round(squeeze(Loc(1,:)),1),round(X,1));
[~,idy] = ismember(round(squeeze(Loc(2,:)),1),round(Y,1));
idxx = sub2ind(size(map),idx,idy);
f = map(idxx);
fAB = f+f';
dist = pdist2(Loc',Loc');
S = ((1-exp(-1*dist/dc))./(sqrt(2)*sqrt((1+exp(-1*dist/dc))))).*fAB;

end