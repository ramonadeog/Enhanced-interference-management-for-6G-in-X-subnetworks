function [] = main_run(X,Y,option)

%% Simulate and generate performance results for distributed channel
%% Selection algorithms with a single (per config.) switching decision
%% threshold
envConstant = envConstants();
envConstant.XBound = X;
envConstant.YBound = Y;
envConstant.mapXPoints = 0:0.1:envConstant.XBound; %linspace(0,envConstant.XBound,envConstant.XBound*10);
envConstant.mapYPoints = 0:0.1:envConstant.YBound;
numSim = ceil(envConstant.simT/envConstant.transTime);


WIRT = GenerateWIRTNetwork(envConstant,numSim,'nfreeway');
WIRT.envConstant = envConstant;
numE = envConstant.numCell*envConstant.numDev;
outageTarget = envConstant.outageTarget;
%numWarmUp = envConstant.warmUp/envConstant.sampTime;
initOption = 'nrandom';
numSensingSteps = envConstant.sampTime/envConstant.transTime;
%numSignallingSteps = envConstant.sigTime/envConstant.transTime;
nCDIndx = 2;
N = envConstant.numCell;
M = envConstant.numGroups;
SINRth = envConstant.SINRvec(nCDIndx);
%numAlgo = length(envConstant.algo);
CSF = 0;

interMat = (squeeze(WIRT.rxPower(1:envConstant.numCell,1:envConstant.numCell,1)));
action.CGC = centralizedColoring(abs(interMat),envConstant.numGroups,'greedy'); 

if strcmp(initOption, 'random')
    action.NNCA = WIRT.channel;
    action.Greedy = WIRT.channel;
    action.GreedyR = WIRT.channel;
    action.minSINR = WIRT.channel;
    action.NNCAS = WIRT.channel;
    action.RDM = WIRT.channel;
else
    action.NNCA = action.CGC;
    action.Greedy = action.CGC;
    action.GreedyR = action.CGC;
    action.minSINR = action.CGC;
    action.NNCAS = action.CGC;
    action.RDM = action.CGC;
end
for indxAlgo = 1:length(envConstant.algo)
    for sIndx = 1:envConstant.numCell
        algo(indxAlgo).IndxI{sIndx} = [];
    end
    TBC{indxAlgo} = [];
end
if strcmp(option,'single')
    disp('Single SINR threshold option selected')
    action_hist = zeros(6,envConstant.numCell,numSim);
    switch_indicator = 1;
    mapCt = 1;
    for stepCt = 1:envConstant.numMap*numSim
        if mod(stepCt,numSim) == 1 && stepCt ~= 1
            mapCt = 1;
            clear WIRT

            for plIndx = 1:length(envConstant.algo)
                CSFn(plIndx) = numel(find(diff(squeeze(action_hist(plIndx,:,:)),1,2)~=0))/numel(squeeze(action_hist(plIndx,:,:)));
            end
            CSF = CSF+CSFn;
            for plIndx = 1:length(envConstant.algo)
                CC = [];
                for sIndx = 1:envConstant.numCell
                    CC= [CC diff(find(diff(squeeze(action_hist(plIndx,sIndx,:)))~=0))'];
                end
                TBC{plIndx} = [TBC{plIndx}; CC(:)*envConstant.transTime];
            end
            save(horzcat('resultsn1',num2str(X),num2str(Y)),'perfResult','CSF','TBC','numE','envConstant','-v7')
            WIRT = GenerateWIRTNetwork(envConstant,numSim,'nfreeway');
            WIRT.envConstant = envConstant;
            interMat = (squeeze(WIRT.rxPower(1:envConstant.numCell,1:envConstant.numCell,1)));
            action.CGC = centralizedColoring(abs(interMat),envConstant.numGroups,'greedy'); 
            if strcmp(initOption, 'random')
                action.NNCA = WIRT.channel;
                action.Greedy = WIRT.channel;
                action.GreedyR = WIRT.channel;
                action.minSINR = WIRT.channel;
                action.NNCAS = WIRT.channel;
                action.RDM = WIRT.channel;
            else
                action.NNCA = action.CGC;
                action.Greedy = action.CGC;
                action.GreedyR = action.CGC;
                action.minSINR = action.CGC;
                action.NNCAS = action.CGC;
                action.RDM = action.CGC;
            end
            disp('Environment updated')
        end
        switch_delay = randi(envConstant.maxDelay,envConstant.numCell,1);
        % Calculate uplink/downlink SINR and aggregate interference power

       for indxAlgo = 1:length(envConstant.algo)
            if indxAlgo == 1
                [algo(indxAlgo).upSINR,algo(indxAlgo).downSINR,algo(indxAlgo).intPow,~] = computeSINRAllRep(WIRT,...
                    action.RDM,mapCt,algo(indxAlgo).IndxI,'nonrandom');
                algo(indxAlgo).upSINRCC = cumsum(algo(indxAlgo).upSINR,3);
                algo(indxAlgo).downSINRCC = cumsum(algo(indxAlgo).downSINR,3);
                algo(indxAlgo).minSINR = min([(min(algo(indxAlgo).upSINRCC(:,:,1),[],2)) (min(algo(indxAlgo).downSINRCC(:,:,1),[],2))],[],2);
            elseif indxAlgo == 2
                [algo(indxAlgo).upSINR,algo(indxAlgo).downSINR,algo(indxAlgo).intPow,~] = computeSINRAllRep(WIRT,...
                    action.Greedy,mapCt,algo(indxAlgo).IndxI,'nonrandom');
                algo(indxAlgo).upSINRCC = cumsum(algo(indxAlgo).upSINR,3);
                algo(indxAlgo).downSINRCC = cumsum(algo(indxAlgo).downSINR,3);
                algo(indxAlgo).minSINR = min([(min(algo(indxAlgo).upSINRCC(:,:,1),[],2)) (min(algo(indxAlgo).downSINRCC(:,:,1),[],2))],[],2);
            elseif indxAlgo == 3
                [algo(indxAlgo).upSINR,algo(indxAlgo).downSINR,algo(indxAlgo).intPow,algo(indxAlgo).IndxI] = computeSINRAllRep(WIRT,...
                    action.NNCA,mapCt,algo(indxAlgo).IndxI,'nonrandom');
                algo(indxAlgo).upSINRCC = cumsum(algo(indxAlgo).upSINR,3);
                algo(indxAlgo).downSINRCC = cumsum(algo(indxAlgo).downSINR,3);
                algo(indxAlgo).minSINR = min([(min(algo(indxAlgo).upSINRCC(:,:,1),[],2)) (min(algo(indxAlgo).downSINRCC(:,:,1),[],2))],[],2);
            elseif indxAlgo == 4
                [algo(indxAlgo).upSINR,algo(indxAlgo).downSINR,algo(indxAlgo).intPow,algo(indxAlgo).IndxI] = computeSINRAllRep(WIRT,...
                    action.GreedyR,mapCt,algo(indxAlgo).IndxI,'nonrandom');
                algo(indxAlgo).upSINRCC = cumsum(algo(indxAlgo).upSINR,3);
                algo(indxAlgo).downSINRCC = cumsum(algo(indxAlgo).downSINR,3);
                algo(indxAlgo).minSINR = min([(min(algo(indxAlgo).upSINRCC(:,:,1),[],2)) (min(algo(indxAlgo).downSINRCC(:,:,1),[],2))],[],2);
            elseif indxAlgo == 5
                [algo(indxAlgo).upSINR,algo(indxAlgo).downSINR,algo(indxAlgo).intPow,~] = computeSINRAllRep(WIRT,...
                    action.NNCAS,mapCt,algo(indxAlgo).IndxI,'nonrandom');
                algo(indxAlgo).upSINRCC = cumsum(algo(indxAlgo).upSINR,3);
                algo(indxAlgo).downSINRCC = cumsum(algo(indxAlgo).downSINR,3);
                algo(indxAlgo).minSINR = min([(min(algo(indxAlgo).upSINRCC(:,:,1),[],2)) (min(algo(indxAlgo).downSINRCC(:,:,1),[],2))],[],2);
            elseif indxAlgo == 6
                [algo(indxAlgo).upSINR,algo(indxAlgo).downSINR,algo(indxAlgo).intPow,~] = computeSINRAllRep(WIRT,...
                    action.CGC,mapCt,algo(indxAlgo).IndxI,'nonrandom');
                algo(indxAlgo).upSINRCC = cumsum(algo(indxAlgo).upSINR,3);
                algo(indxAlgo).downSINRCC = cumsum(algo(indxAlgo).downSINR,3);
                algo(indxAlgo).minSINR = min([(min(algo(indxAlgo).upSINRCC(:,:,1),[],2)) (min(algo(indxAlgo).downSINRCC(:,:,1),[],2))],[],2);
            end
        end
        % Compute estimate of the number of failed loops
        for nC = 1:envConstant.numBWConf
            Pmap = WIRT.envConstant.outageMap.Pout{16-(nC)+1}; 
            Smap = WIRT.envConstant.outageMap.SINR{16-(nC)+1};
            parfor indxAlgo = 1:length(envConstant.algo)
                uplinkPoutCC = probMapping(Pmap,Smap,10*log10(algo(indxAlgo).upSINRCC));
                downlinkPoutCC = probMapping(Pmap,Smap,10*log10(algo(indxAlgo).downSINRCC));
                upIndxCC = find(uplinkPoutCC>outageTarget);
                downIndxCC = find(downlinkPoutCC>outageTarget);
                nSDevCC(indxAlgo,nC) = numel(ismember(upIndxCC,downIndxCC));
            end
        end
        perfResult(:,:,stepCt) = nSDevCC;
        action_hist(1,:,mapCt) = action.RDM;
        action_hist(2,:,mapCt) = action.Greedy;
        action_hist(5,:,mapCt) = action.NNCAS;
        action_hist(3,:,mapCt) = action.NNCA;
        action_hist(4,:,mapCt) = action.GreedyR;
        action_hist(6,:,mapCt) = action.CGC;
        % Check sensing interval condition and begin channel selection
        % operations
        if switch_indicator >= envConstant.maxDelay
            switch_indicator = 1;
        end
        if mapCt >= numSensingSteps
    %         gwloc = squeeze(WIRT.gwLoc(:,:,mapCt)).';
    %         gwdist = pdist2(gwloc,gwloc);
    %         D = processRank(gwdist,algo(5).intPow);
            interMat = (squeeze(WIRT.rxPower(1:envConstant.numCell,1:envConstant.numCell,mapCt)));
            for indxAlgo = 1:length(envConstant.algo)
                for indxCell = 1:N
                    if algo(indxAlgo).minSINR(indxCell) <= SINRth && switch_delay(indxCell) == switch_indicator
                        if indxAlgo == 1
                            action.RDM(indxCell) = randi(M,1);
                        elseif indxAlgo == 2
                            [~,action.Greedy(indxCell)] = min(algo(indxAlgo).intPow(indxCell,:));
                        elseif indxAlgo == 3
                                action.NNCA(indxCell) = nearestNeighbourAvoidance(squeeze(interMat(indxCell,:)),...
                                    action.NNCA,envConstant.numGroups,'nonaggregate');   
                        elseif indxAlgo == 4
    %                         action.GreedyR(indxCell) = greedyRankExchange(squeeze(D(indxCell,:,:)));
                              [~,action.GreedyR(indxCell)] = min(algo(indxAlgo).intPow(indxCell,:));
                        elseif indxAlgo == 5
                            action.NNCAS(indxCell) = nearestNeighbourAvoidance(squeeze(interMat(indxCell,:)),...
                                action.NNCAS,envConstant.numGroups,'nonaggregate');

                        end

                    end

                end
            end              
            action.CGC = centralizedColoring(abs(interMat),envConstant.numGroups,'greedy'); 
            switch_indicator = switch_indicator + 1;
        end
        mapCt = mapCt +1;
        if mod(stepCt,100) == 0
            disp(['Number of steps completed:', num2str(stepCt)]);  
        end

    end 

    save(horzcat('resultsn1',num2str(X),num2str(Y)),'perfResult','CSF','TBC','numE','envConstant','-v7')
else
    disp('Per BW SINR threshold option selected')
    action_hist = zeros(6,envConstant.numBWConf,envConstant.numCell,numSim);
    CSF = zeros(6,envConstant.numBWConf);
    for nC = 1:envConstant.numBWConf
        switch_indicator = 1;
        mapCt = 1;
        SINRth = envConstant.SINRvec(nC);
        for stepCt = 1:envConstant.numMap*numSim
            if mod(stepCt,numSim) == 1 && stepCt ~= 1
                mapCt = 1;
                %clear WIRT
                for plIndx = 1:length(envConstant.algo)
                    CSFn(plIndx,nC) = numel(find(diff(squeeze(action_hist(plIndx,nC,:,:)),1,2)~=0))/numel(squeeze(action_hist(plIndx,nC,:,:)));
                end
                CSF(:,nC) = CSF(:,nC)+CSFn(:,nC);
                for plIndx = 1:length(envConstant.algo)
                    CC = [];
                    for sIndx = 1:envConstant.numCell
                        CC= [CC diff(find(diff(squeeze(action_hist(plIndx,nC,sIndx,:)))~=0))'];
                    end
                    TBC{plIndx,nC} = [TBC{plIndx}; CC(:)*envConstant.transTime];
                end
                save(horzcat('resultsMBW',num2str(X),num2str(Y)),'perfResult','CSF','TBC','numE','envConstant','-v7')
              
                WIRT = GenerateWIRTNetwork(envConstant,numSim,'nfreeway');
                WIRT.envConstant = envConstant;
                interMat = (squeeze(WIRT.rxPower(1:envConstant.numCell,1:envConstant.numCell,1)));
                action.CGC = centralizedColoring(abs(interMat),envConstant.numGroups,'greedy'); 
                if strcmp(initOption, 'random')
                    action.NNCA = WIRT.channel;
                    action.Greedy = WIRT.channel;
                    action.GreedyR = WIRT.channel;
                    action.minSINR = WIRT.channel;
                    action.NNCAS = WIRT.channel;
                    action.RDM = WIRT.channel;
                else
                    action.NNCA = action.CGC;
                    action.Greedy = action.CGC;
                    action.GreedyR = action.CGC;
                    action.minSINR = action.CGC;
                    action.NNCAS = action.CGC;
                    action.RDM = action.CGC;
                end
                disp('Environment updated')
            end
            switch_delay = randi(envConstant.maxDelay,envConstant.numCell,1);
            % Calculate uplink/downlink SINR and aggregate interference power

           for indxAlgo = 1:length(envConstant.algo)
                if indxAlgo == 1
                    [algo(indxAlgo).upSINR,algo(indxAlgo).downSINR,algo(indxAlgo).intPow,~] = computeSINRAllRep(WIRT,...
                        action.RDM,mapCt,algo(indxAlgo).IndxI,'nonrandom');
                    algo(indxAlgo).upSINRCC = cumsum(algo(indxAlgo).upSINR,3);
                    algo(indxAlgo).downSINRCC = cumsum(algo(indxAlgo).downSINR,3);
                    algo(indxAlgo).minSINR = min([(min(algo(indxAlgo).upSINRCC(:,:,1),[],2)) (min(algo(indxAlgo).downSINRCC(:,:,1),[],2))],[],2);
                elseif indxAlgo == 2
                    [algo(indxAlgo).upSINR,algo(indxAlgo).downSINR,algo(indxAlgo).intPow,~] = computeSINRAllRep(WIRT,...
                        action.Greedy,mapCt,algo(indxAlgo).IndxI,'nonrandom');
                    algo(indxAlgo).upSINRCC = cumsum(algo(indxAlgo).upSINR,3);
                    algo(indxAlgo).downSINRCC = cumsum(algo(indxAlgo).downSINR,3);
                    algo(indxAlgo).minSINR = min([(min(algo(indxAlgo).upSINRCC(:,:,1),[],2)) (min(algo(indxAlgo).downSINRCC(:,:,1),[],2))],[],2);
                elseif indxAlgo == 3
                    [algo(indxAlgo).upSINR,algo(indxAlgo).downSINR,algo(indxAlgo).intPow,algo(indxAlgo).IndxI] = computeSINRAllRep(WIRT,...
                        action.NNCA,mapCt,algo(indxAlgo).IndxI,'nonrandom');
                    algo(indxAlgo).upSINRCC = cumsum(algo(indxAlgo).upSINR,3);
                    algo(indxAlgo).downSINRCC = cumsum(algo(indxAlgo).downSINR,3);
                    algo(indxAlgo).minSINR = min([(min(algo(indxAlgo).upSINRCC(:,:,1),[],2)) (min(algo(indxAlgo).downSINRCC(:,:,1),[],2))],[],2);
                elseif indxAlgo == 4
                    [algo(indxAlgo).upSINR,algo(indxAlgo).downSINR,algo(indxAlgo).intPow,algo(indxAlgo).IndxI] = computeSINRAllRep(WIRT,...
                        action.GreedyR,mapCt,algo(indxAlgo).IndxI,'nonrandom');
                    algo(indxAlgo).upSINRCC = cumsum(algo(indxAlgo).upSINR,3);
                    algo(indxAlgo).downSINRCC = cumsum(algo(indxAlgo).downSINR,3);
                    algo(indxAlgo).minSINR = min([(min(algo(indxAlgo).upSINRCC(:,:,1),[],2)) (min(algo(indxAlgo).downSINRCC(:,:,1),[],2))],[],2);
                elseif indxAlgo == 5
                    [algo(indxAlgo).upSINR,algo(indxAlgo).downSINR,algo(indxAlgo).intPow,~] = computeSINRAllRep(WIRT,...
                        action.NNCAS,mapCt,algo(indxAlgo).IndxI,'nonrandom');
                    algo(indxAlgo).upSINRCC = cumsum(algo(indxAlgo).upSINR,3);
                    algo(indxAlgo).downSINRCC = cumsum(algo(indxAlgo).downSINR,3);
                    algo(indxAlgo).minSINR = min([(min(algo(indxAlgo).upSINRCC(:,:,1),[],2)) (min(algo(indxAlgo).downSINRCC(:,:,1),[],2))],[],2);
                elseif indxAlgo == 6
                    [algo(indxAlgo).upSINR,algo(indxAlgo).downSINR,algo(indxAlgo).intPow,~] = computeSINRAllRep(WIRT,...
                        action.CGC,mapCt,algo(indxAlgo).IndxI,'nonrandom');
                    algo(indxAlgo).upSINRCC = cumsum(algo(indxAlgo).upSINR,3);
                    algo(indxAlgo).downSINRCC = cumsum(algo(indxAlgo).downSINR,3);
                    algo(indxAlgo).minSINR = min([(min(algo(indxAlgo).upSINRCC(:,:,1),[],2)) (min(algo(indxAlgo).downSINRCC(:,:,1),[],2))],[],2);
                end
            end
            % Compute estimate of the number of failed loops

            Pmap = WIRT.envConstant.outageMap.Pout{16-(nC)+1}; 
            Smap = WIRT.envConstant.outageMap.SINR{16-(nC)+1};
            for indxAlgo = 1:length(envConstant.algo)
                uplinkPoutCC = probMapping(Pmap,Smap,10*log10(algo(indxAlgo).upSINRCC));
                downlinkPoutCC = probMapping(Pmap,Smap,10*log10(algo(indxAlgo).downSINRCC));
                upIndxCC = find(uplinkPoutCC>outageTarget);
                downIndxCC = find(downlinkPoutCC>outageTarget);
                nSDevCC(indxAlgo) = numel(ismember(upIndxCC,downIndxCC));
            end
            
            perfResult(:,nC,stepCt) = nSDevCC;
            action_hist(1,nC,:,mapCt) = action.RDM;
            action_hist(2,nC,:,mapCt) = action.Greedy;
            action_hist(5,nC,:,mapCt) = action.NNCAS;
            action_hist(3,nC,:,mapCt) = action.NNCA;
            action_hist(4,nC,:,mapCt) = action.GreedyR;
            action_hist(6,nC,:,mapCt) = action.CGC;
            % Check sensing interval condition and begin channel selection
            % operations
            if switch_indicator >= envConstant.maxDelay
                switch_indicator = 1;
            end
            if mapCt >= numSensingSteps
        %         gwloc = squeeze(WIRT.gwLoc(:,:,mapCt)).';
        %         gwdist = pdist2(gwloc,gwloc);
        %         D = processRank(gwdist,algo(5).intPow);
                interMat = (squeeze(WIRT.rxPower(1:envConstant.numCell,1:envConstant.numCell,mapCt)));
                for indxAlgo = 1:length(envConstant.algo)
                    for indxCell = 1:N
                        if algo(indxAlgo).minSINR(indxCell) <= SINRth && switch_delay(indxCell) == switch_indicator
                            if indxAlgo == 1
                                action.RDM(indxCell) = randi(M,1);
                            elseif indxAlgo == 2
                                [~,action.Greedy(indxCell)] = min(algo(indxAlgo).intPow(indxCell,:));
                            elseif indxAlgo == 3
                                    action.NNCA(indxCell) = nearestNeighbourAvoidance(squeeze(interMat(indxCell,:)),...
                                        action.NNCA,envConstant.numGroups,'nonaggregate');   
                            elseif indxAlgo == 4
        %                         action.GreedyR(indxCell) = greedyRankExchange(squeeze(D(indxCell,:,:)));
                                  [~,action.GreedyR(indxCell)] = min(algo(indxAlgo).intPow(indxCell,:));
                            elseif indxAlgo == 5
                                action.NNCAS(indxCell) = nearestNeighbourAvoidance(squeeze(interMat(indxCell,:)),...
                                    action.NNCAS,envConstant.numGroups,'nonaggregate');

                            end

                        end

                    end
                end              
                action.CGC = centralizedColoring(abs(interMat),envConstant.numGroups,'greedy'); 
                switch_indicator = switch_indicator + 1;
            end
            mapCt = mapCt +1;
            if mod(stepCt,100) == 0
                disp(['Number of steps completed:', num2str(stepCt),' for BW Index :' num2str(nC)]);  
            end

        end 
    end

    save(horzcat('resultsMBW',num2str(X),num2str(Y)),'perfResult','CSF','TBC','numE','envConstant','-v7')
    
end
 
end

