function [envConstant] = envConstants()
%Specify environment and other paramaters
envConstant.outageMap = load('DataPoutInterp480kHz.mat');

envConstant.XBound = 40;
envConstant.YBound = 40;
envConstant.speed = 2; 
envConstant.transTime = 0.005;
envConstant.sampTime = 0.050;
envConstant.sigTime = 0.05;
envConstant.maxDelay = 8;
envConstant.TTA = [1,8]; 
envConstant.simT = 20;
envConstant.updateT = 20; 
envConstant.numBWConf = 5; 
envConstant.warmUp = 1;
envConstant.numMap = 100;
envConstant.numCell = 16;
envConstant.numDev = 18;
envConstant.cellRad = 2.5;
envConstant.sigmaS = 3; 
envConstant.dc = 4;
envConstant.algo = {'Random','Greedy','NNCA - Rep','Greedy - Rep','NNCA','CGC'};
envConstant.shadowFadingModel = 'map';
envConstant.gammaP = [2.2 2.2];
envConstant.sigmaVec = [3 3];
envConstant.noiseFactor = 10;
envConstant.txPow = -10;
envConstant.numChannel = 12;
envConstant.numGroups = 6;
envConstant.numRep = 2;
envConstant.currentLoc = 0;
envConstant.mapXPoints = 0:0.1:envConstant.XBound; %linspace(0,envConstant.XBound,envConstant.XBound*10);
envConstant.mapYPoints = 0:0.1:envConstant.YBound; % linspace(0,envConstant.YBound,envConstant.YBound*10);
envConstant.carrierFreq = 6e9;
envConstant.bandW = 20e6;
envConstant.outageTarget = 1e-6;
envConstant.noisePow = 10^((-174+envConstant.noiseFactor...
    +10*log10(envConstant.bandW))/10);
envConstant.SINRth = db2pow(10);
envConstant.SINRvec = db2pow([34 17 10.5 7 4.5 2.5 1.2 0.5]+3);
envConstant.Sth = 5;
envConstant.Count = 0;
if envConstant.numGroups == 4 && envConstant.numRep == 3
%%12  channels 3 repetition
    envConstant.frames = [1     4     3     6     2     5;
       2     5     1     4     3     6;
        3     6     2     5     1     4];
    envConstant.frameS = [envConstant.frames envConstant.frames+6 envConstant.frames+12;...
        envConstant.frames+18 envConstant.frames+24 envConstant.frames+30;...
        envConstant.frames+36 envConstant.frames+42 envConstant.frames+48; envConstant.frames+54 ...
        envConstant.frames+60 envConstant.frames+66];
elseif envConstant.numGroups == 3 && envConstant.numRep == 4
%%12  channels 3 groups 3 repetition
    envConstant.frames=[1 5 4 8 2 6 3 7;2 6 3 7 1 5 4 8; 3 7 1 5 4 8 2 6; 4 8 2 6 3 7 1 5];
    envConstant.frameS = [envConstant.frames envConstant.frames+8 envConstant.frames+16;.... 
        envConstant.frames+24 envConstant.frames+32 envConstant.frames+40;...
        envConstant.frames+48 envConstant.frames+56 envConstant.frames+64];
elseif envConstant.numGroups == 3 && envConstant.numRep == 3
    envConstant.frames=[1 5 4 8 2 6 3 7;2 6 3 7 1 5 4 8; 3 7 1 5 4 8 2 6; 400 800 200 600 300 700 100 500];
    envConstant.frameS = [envConstant.frames envConstant.frames+8 envConstant.frames+16;.... 
        envConstant.frames+24 envConstant.frames+32 envConstant.frames+40;...
        envConstant.frames+48 envConstant.frames+56 envConstant.frames+64];

elseif envConstant.numGroups == 6 && envConstant.numRep == 2
    envConstant.frames = [1     4     3     6     2     5;
       2     5     1     4     3     6];
    envConstant.frameS = [envConstant.frames envConstant.frames+6 envConstant.frames+12;...
        envConstant.frames+18 envConstant.frames+24 envConstant.frames+30;...
        envConstant.frames+36 envConstant.frames+42 envConstant.frames+48; envConstant.frames+54 ...
        envConstant.frames+60 envConstant.frames+66; envConstant.frames+72 envConstant.frames+78 ...
        envConstant.frames+84; envConstant.frames+90 envConstant.frames+96 envConstant.frames+102];
   
elseif envConstant.numGroups == 12 && envConstant.numRep == 1
    envConstant.frames = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18];
    envConstant.frameS = [envConstant.frames;envConstant.frames+18;envConstant.frames+36;...
        envConstant.frames+54; envConstant.frames+72; envConstant.frames+90;...
        envConstant.frames+108;envConstant.frames+126;envConstant.frames+144;...
        envConstant.frames+162; envConstant.frames+180; envConstant.frames+198];

else
    error('Specified number of groups and repetitions not supported');
end


end



% function [envConstant] = envConstants()
% %Specify environment and other paramaters
% envConstant.outageMap = load('DataPoutInterp480kHz.mat');
% %%12  channels 3 repetition
% % envConstant.frames = [1     4     3     6     2     5;
% %    2     5     1     4     3     6;
% %     3     6     2     5     1     4];
% % envConstant.frameS = [envConstant.frames envConstant.frames+6 envConstant.frames+12;...
% %     envConstant.frames+18 envConstant.frames+24 envConstant.frames+30;...
% %     envConstant.frames+36 envConstant.frames+42 envConstant.frames+48; envConstant.frames+54 ...
% %     envConstant.frames+60 envConstant.frames+66];
% %%12  channels 3 groups 3 repetition
% envConstant.frames=[1 5 4 8 2 6 3 7;2 6 3 7 1 5 4 8; 3 7 1 5 4 8 2 6; 4 8 2 6 3 7 1 5];
% %envConstant.frames=[1 5 4 8 2 6 3 7;2 6 3 7 1 5 4 8; 3 7 1 5 4 8 2 6; 400 800 200 600 300 700 100 500];
% envConstant.frameS = [envConstant.frames envConstant.frames+8 envConstant.frames+16;.... 
%     envConstant.frames+24 envConstant.frames+32 envConstant.frames+40;...
%     envConstant.frames+48 envConstant.frames+56 envConstant.frames+64];
% 
% %12 channels 6 groups 2 repetition
% % envConstant.frames = [1     4     3     6     2     5;
% %    2     5     1     4     3     6];
% % envConstant.frameS = [envConstant.frames envConstant.frames+6 envConstant.frames+12;...
% %     envConstant.frames+18 envConstant.frames+24 envConstant.frames+30;...
% %     envConstant.frames+36 envConstant.frames+42 envConstant.frames+48; envConstant.frames+54 ...
% %     envConstant.frames+60 envConstant.frames+66; envConstant.frames+72 envConstant.frames+78 ...
% %     envConstant.frames+84; envConstant.frames+90 envConstant.frames+96 envConstant.frames+102];
% envConstant.XBound = 30;
% envConstant.YBound = 30;
% envConstant.speed = 2; 
% envConstant.sampTime = 0.05;
% envConstant.simT = 2;
% envConstant.updateT = 1; 
% envConstant.warmUp = 0;
% envConstant.numCell = 16;
% envConstant.numDev = 18;
% envConstant.cellRad = 2.5;
% envConstant.sigmaS = 3;
% envConstant.dc = 4;
% envConstant.shadowFadingModel = 'map';
% envConstant.gammaP = 2.2;
% envConstant.noiseFactor = 10;
% envConstant.txPow = -10;
% envConstant.numChannel = 12;
% envConstant.numGroups = 3;
% envConstant.numRep = 4;
% envConstant.currentLoc = 0;
% envConstant.mapXPoints = 0:0.1:envConstant.XBound; %linspace(0,envConstant.XBound,envConstant.XBound*10);
% envConstant.mapYPoints = 0:0.1:envConstant.YBound; % linspace(0,envConstant.YBound,envConstant.YBound*10);
% envConstant.carrierFreq = 6e9;
% envConstant.bandW = 20e6;
% envConstant.outageTarget = 1e-6;
% envConstant.noisePow = 10^((-174+envConstant.noiseFactor...
%     +10*log10(envConstant.bandW))/10);
% envConstant.SINRth = db2pow(5);
% envConstant.Count = 0;
% end
% 
