function [action] = nearestNeighbourAvoidance(pow,channel,numChann,methods)
%% =======================Nearest neighbour conflict avoidance ==========
%===================== Ramoni Adeogun [AAU, 2019] =======================
%% ======================================================================
%
chan_list = 1:numChann;
%rng(seed); 
switch methods
    case 'aggregate'
      
        [~,Indx] = sort(pow,'descend');
        action = Indx(1);  
        %display(action)   %debugging
    otherwise       
        [~,Indx] = sort(pow,'descend');
        blocked_chan = channel(Indx(2:numChann));
        av_chan = chan_list(~ismember(chan_list,blocked_chan));
        if length(av_chan) > 1
            action = av_chan(randi(length(av_chan)));
        else
            action = av_chan;
        end       
end
action = action(:)';

return
