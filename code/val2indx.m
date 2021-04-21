
function [i,j] = val2indx(ix, a, b)
% val2indexes provides the indexes i,j from index ix
%   OUTPUT: 
%       * i - index of the first element
%       * j - index of the second element
%   INPUT: 
%       * ix - index of the (i,j) combination
%       * a - size of elements containing "i"
%       * b - size of elements containing "j"

i = mod(ix, a); 
if i == 0, i = a; end    
y = mod(ix, (a * b));
j = ceil(y/a);
if j == 0, j = b; end 
    
end

