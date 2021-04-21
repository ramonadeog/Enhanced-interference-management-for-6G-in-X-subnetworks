function [WIRT] = GenerateWIRTNetwork(envConstant,numSim,option)
%% =============== GENERATEWIRTNETWORK ===================================
% == Create mobility, compute shadowing and pairwise receive/interference
% power
%==================Ramoni Adeogun [AAU, 2020] ============================
%% ========================================================================

%rng('shuffle');
N = envConstant.numCell+envConstant.numCell*envConstant.numDev;
if strcmp(option,'freeway')
    [gwLoc,angD] = freewayModel(envConstant.XBound,envConstant.YBound,envConstant.transTime,envConstant.simT,[20 20],envConstant.numCell,envConstant.cellRad);
   
else
    [gwLoc,angD,~,~] = commonWaypointMM(envConstant.XBound,envConstant.YBound,...
        envConstant.transTime,envConstant.simT,envConstant.speed,2*pi,envConstant.numCell,...
        1.5*envConstant.cellRad,envConstant.cellRad);
end

map = createMap(envConstant.mapXPoints,envConstant.mapYPoints,....
    envConstant.sigmaS,envConstant.dc);

gwLocInit = squeeze(gwLoc(:,:,1)); 
devLocInit = depDevicesCircCell(gwLocInit',envConstant.cellRad,...
    envConstant.numDev); 

% %Initial device locations
% devLocX = squeeze(devLocInit(1,:,:))';
% devLocX = devLocX(:);
% devLocY = squeeze(devLocInit(2,:,:))';
% devLocY = devLocY(:);

 %diffLoc = squeeze(gwLoc(:,:,:))-gwLocInit;
% LocInit = [squeeze(gwLocInit(1,:))' squeeze(gwLocInit(2,:))'; devLocX devLocY];
%Locn = zeros(N,N,numSim); 
rxPow = zeros(N,N,numSim);
switch envConstant.shadowFadingModel
    case 'map'
    for n = 1:numSim
        gwLocn = squeeze(gwLoc(:,:,n));
        devLocn = devLocInit+reshape(repmat((squeeze(gwLoc(:,:,n))-gwLocInit)...
            ,1,envConstant.numDev),2,envConstant.numCell,envConstant.numDev); %squeeze(diffLoc(:,:,n));
        devLocXn = squeeze(devLocn(1,:,:))';
        devLocXn = devLocXn(:);
        devLocYn = squeeze(devLocn(2,:,:))';
        devLocYn = devLocYn(:);
        Loc = [squeeze(gwLocn(1,:))' squeeze(gwLocn(2,:))'; devLocXn devLocYn];
        %Locn(:,:,n) = Loc;
        [S,dist] = computeShadowing(map,envConstant.mapXPoints,envConstant.mapYPoints,envConstant.dc,Loc');
        rxPow(:,:,n) = computePowMM(S,dist,envConstant.sigmaVec,envConstant.txPow,envConstant.gammaP...
            ,envConstant.carrierFreq,envConstant.bandW,1,envConstant.numCell,envConstant.numDev);
        if mod(n,ceil(envConstant.updateT/envConstant.sampTime)) == 0
            map = createMap(envConstant.mapXPoints,envConstant.mapYPoints,....
        envConstant.sigmaS,envConstant.dc);
        end

    end
    otherwise
        error('unsurpoted option');
%     for n = 1:numSim
%         gwLocn = squeeze(gwLoc(:,:,n));
%         devLocn = devLocInit+squeeze(diffLoc(:,:,n));
%         devLocXn = squeeze(devLocn(1,:,:))';
%         devLocXn = devLocXn(:);
%         devLocYn = squeeze(devLocn(2,:,:))';
%         devLocYn = devLocYn(:);
%         Loc = [squeeze(gwLocn(1,:))' squeeze(gwLocn(2,:))'; devLocXn devLocYn];
%         %Locn(:,:,n) = Loc;
%         if n == 1
%            %L = size(Loc,1);
%            dist = pdist2(Loc,Loc);
%            C_shadow = exp(-1*dist/envConstant.dc);
% %            S = mvnrnd(zeros(1,L),eye(size(C_shadow))*envConstant.sigmaS,L);
% %            S = triu(S)+triu(S ,1)';
%            S = randn(size(C_shadow))*envConstant.sigmaS;
%            rxPow(:,:,n) = computePowMM(S,dist,envConstant.txPow,envConstant.gammaP...
%             ,envConstant.carrierFreq,envConstant.bandW,1);
%         else
%            [S,dist] = computeShadowTempCorr(S,Loc,dist,envConstant.dc,envConstant.sigmaS);
%            rxPow(:,:,n) = computePowMM(S,dist,envConstant.txPow,envConstant.gammaP...
%             ,envConstant.carrierFreq,envConstant.bandW,1);
%         end
% 
%     end
end
        
        
WIRT.gwLoc = gwLoc;
WIRT.angD = angD;
%WIRT.location = Locn;
WIRT.rxPower = rxPow;
% Loc = WIRT.gwLoc(:,:,1).'; 
% pdistance = pdist2(Loc,Loc);
WIRT.channel = randi(envConstant.numGroups,1,envConstant.numCell); %centralizedColoring(pdistance,envConstant.numGroups);%
%end
return 