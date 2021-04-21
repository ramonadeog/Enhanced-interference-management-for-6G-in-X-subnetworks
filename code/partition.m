classdef partition
    properties
        param
        array
    end
    methods
        function p = partition(param)
 
             if nargin==0
                 param = 0;
             end

             if ~isa(param,'cell')
                 n = param;
                 p.array = logical(speye(n,n));  % singleton parts
                 p.array = sortrows(p.array);
                 p = class(p,'partition');
                 return
             end

             % otherwise, called with a cell array containing the parts of p

             % special case: param is an empty cell array
             if isempty(param)
                 p.array = speye(0);
                 p = class(p,'partition');
                 return
             end
             % Find nv (n) and np (m)
             maxv = 0;
             for i=1:length(param)
                 m = max(param{i});
                 maxv = max([m,maxv]);
             end
             n = maxv;
             m = length(param);
             p.array = logical(sparse([],[],[],m,n,n));
             % load the entries
             for i=1:m
                 row = zeros(1,n);
                 row(param{i}) = 1;
                 p.array(i,:) = logical(row);
             end
             p.array = sortrows(p.array);

             %p = class(p,'partition');

             if ~check(p)
                 p.array = speye(0);
                 error('The cell array does not define a valid partition');
             end
        end
    end
end
 
 function yn = check(p)
 % check(p) --- check that the datastructure holding p is a valid partition
 
 s = full(sum(p.array,1));  % col sums should all be 1
 t = full(sum(p.array,2));  % row sums should all be positive
 yn = all(s == 1)  & all(t>0);
 end