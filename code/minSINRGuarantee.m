function [action] = minSINRGuarantee(int_pow,des_pow,numChan,noise_pow,sinr_th)
%% =========================================================================
% Minimum SINR Guarantee dynamic channel allocation
% Select the worst channel that means SINR requirements
%OUTPUTS:
%action: selected channel
%succ_flag: a binary indicator showing whether a channel meeting the 
%SINR constraint is found (too frequent 0's indicate too high SINR
%threshold)
%INPUTS:
%int_pow: aggregate interference power on all channel groups
%des_pow: worst case desired signal power
%channel: currently occupied channel (or channel group)
%noise_pow: noise power
%sinr_th: SINR threshold 
%% =======================================================================
% Step 1: Compute possible SINR on all channel
sinr_ach = des_pow./(int_pow+noise_pow); 
% Step 2: Sort sinr is increasing order of magnitude
[sinr_sorted,indx_sorted] = sort(sinr_ach);
% Step 3: Select channels satisfying sinr requirement
chann_avai = indx_sorted(sinr_sorted > sinr_th);
% Step 4: Check if there are channels meeting the condition and choose the worst
if isempty(chann_avai)
    chann_list = 1:numChan;
    action = chann_list(sinr_ach == max(sinr_ach));
    if length(action) > 1
        action = action(randi(length(action),1));
    end
    %succ_flag = 0;
else
    action = chann_avai(1);
    %succ_flag = 1; 
    
end

return