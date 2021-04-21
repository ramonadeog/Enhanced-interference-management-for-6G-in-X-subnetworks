%% Plotting performance results
clear all; clc;
X = 30;
Y = 30;
load(horzcat('resultsn1',num2str(X),num2str(Y),'.mat'));
% load(horzcat('results1',num2str(X),num2str(Y),'.mat'));
%%
% perfResult(:,:,1000001:1640000) = perfResult1;
% CSF = CSF1+CSF;
% 
% for i = 1:6
%     TBC{i} = [TBC{i}; TBC1{i}];
% end
%%
%perfResult(:,:,[1:10 20000:20020 40000:40020 60000 60020 80000 80020 100000 100020  ]) = [];
PLF = mean(perfResult,3)/numE;
linestyle = {'-bo','-rs','--gd','--rs','-gd','-kd'};
indBW = 1:envConstant.numBWConf;
tt= figure(); 
for plIndx = 2:length(envConstant.algo)
    semilogy(indBW*40,PLF(plIndx,:),linestyle{plIndx},'linewidth',2,'DisplayName',envConstant.algo{plIndx});  hold on
end
xlabel('Bandwidth per channel [MHz]')
ylabel('PLF');
ylim([7e-7 1]);
xlim([40 envConstant.numBWConf*40]);
legend show
hold off; box off
title(horzcat('Deployment area: ', num2str(X), ' m', ' by ', num2str(Y),' m'));
grid on
savefig(tt,horzcat('figures/PLFResult5',num2str(X),num2str(X)));
exportgraphics(tt,horzcat('figures/PLFResult5',num2str(X),num2str(X),'.pdf'));
%%
for plIndx = 1:length(envConstant.algo)
    [f{plIndx} , x{plIndx} ] = ecdf(TBC{plIndx});
end

cc= figure(); 
for plIndx = 2:length(envConstant.algo)
    semilogy(f{plIndx},x{plIndx},'linewidth',2,'DisplayName',envConstant.algo{plIndx});  hold on
end
hold off; grid off
legend show
xlabel('CDF')
ylabel('Time between channel switching');
grid on; box off;
title(horzcat('Deployment area: ', num2str(X), ' m', ' by ', num2str(Y),' m'));
savefig(cc,horzcat('figures/TBCResult5',num2str(X),num2str(X)));
exportgraphics(cc,horzcat('figures/TBCResult5',num2str(X),num2str(X),'.pdf'));

%%
dd = figure();
algoList = categorical({'Greedy','NNCA - Rep','Greedy - Rep','NNCA','CGC'});
bar(algoList, CSF(2:6)*100/50); 
ylabel('Channel Switching Frequency (%)');
grid on
box off
title(horzcat('Deployment area: ', num2str(X), ' m', ' by', num2str(Y),' m'));
savefig(dd,horzcat('figures/CSFResult5',num2str(X),num2str(X)));
exportgraphics(dd,horzcat('figures/CSFResult5',num2str(X),num2str(X),'.pdf'));
