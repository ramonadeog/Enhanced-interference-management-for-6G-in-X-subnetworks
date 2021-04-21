function [action] = selectChannel(intPow,epsilon,channel,numG,option)
%% ======================================================================
% Epsilon-greedy channel selection 
chann_list = 1:numG;
%rng(seed);
switch option
    %% Greedy selection with 1-step prediction
    case 'filter'
       intPow_pre = squeeze(intPow(:,2));
       intPow_next = intPow_pre+1*(intPow(:,2)-intPow(:,1));
       [~,Indx] = sort(intPow_pre); 
       nextPow_sort = intPow_next(Indx(2:end));
       Indxx = find(nextPow_sort < intPow_pre(Indx(1)));
       if isempty(Indxx)
           action = Indx(1);
       else
           action = Indxx(nextPow_sort(Indxx)==min(nextPow_sort(Indxx)));
       end
       %% Epsilon-greedy channel selection
    otherwise 
        [~,Indx] = sort(intPow);
        if rand() <= epsilon
            action = Indx(1);
            %display(action)
        else
            av_chann = chann_list(chann_list ~= channel);
            action = av_chann(randi(numG-1,1));
        end 
end 
return