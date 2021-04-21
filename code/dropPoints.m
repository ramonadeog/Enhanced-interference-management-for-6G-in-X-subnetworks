function [Loc]=dropPoints(W,H,Num,R0)

%function [X,Y,Nmax,Dmatrix,delay]=scatter_points71(void)

% given H and W size of 2D rectangle scatter_points7.m does the following:
% 1. calculate Nmax, the maximum amount of circles radius R0 that orderly fit in rectangle HxW
% 2. scatter Ap random points with spatial resolution dx =- 0.1 and dy=0.1
%           and spaced at least distance R0
%           if requested amount Ap>Nmax rand2Dpoints breaks because there
%           is no space to fit in so many random points complying with min
%           distance R0
% 3. calculate distance matrix Dmatrix among all points
% 4. return the coordinates of all generated points in X and Y, along with Dmatrix, Nmax, and delay in seconds tic toc around main loop
% 5. plot points
% 6. plot safety circles to visually verify 
% 
% the amount of generated points is numel(X) and cannot be larger than Nmax 
% because of the request for the points to be randomly generated.
% 
% in this initial version, only manual input through message box.
% 
% call examples:
% 1.
% [X,Y,Dmatrix,Nmax]=scatter_points7
% 
% 2.
% [X,Y,Nmax]=scatter_points7
% 2 examples how to verify distance values and check distances meet requirement > R0
% 1.
% L=combinator(Ap,2,'c'); 
% relD2=((X(L(:,2))-X(L(:,1))).^2+(Y(L(:,2))-Y(L(:,1))).^2).^.5
% find(relD2<R0)
% 
% 2.
% L2=combinator(Ap,2);  
% Dmatrix=reshape(((X(L2(:,2))-X(L2(:,1))).^2+(Y(L2(:,2))-Y(L2(:,1))).^2).^.5,[Ap Ap]);
% Dmatrix([1:21:end])=NaN ;
% Dmatrix(Dmatrix<R0)
% 
% February 10th 2017  version: 1.0
% February 16th 2017 version: 1.01 improved speed for large (area>900^2) rectangles with double rand
%        now delay ~1sec regardless of rectangle size, delay can be    returned in additional output.
%        if A<900^2 dx=0.1 dy=0.1  
%        if A>=900^2 circles calculated with R+1
%             and then shake each point [-0.5 0.5]
%  
% author: John Bofarull Guix, any feedback to build next version is welcome at
% jgb2012@sky.com or through the Mathworks website
% this script was inspired by Mr Marwen Tarhoumi marwentarhoumi@rocketmail.com
% and Matt Fig's popkenai@yahoo.com mighty function combinator.m available
% from Mathworks File exchange

%clear all;clc;close all;
%format bank;
rng('Shuffle');
delay=0;

%prompt = {'rectangle width W: ','rectange height H: ','amount points to scatter: ','safety distance Radius: '};
%dlg_title = 'Input';
%num_lines = 1;
%defaultans = {'100','100','20','3'};
%input_answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

Diam0=2*R0;
W=W+Diam0;
H=H+Diam0;
Ap=Num;
Ap=floor(Ap);
Nmax=calc_amount_circles(W-Diam0,H-Diam0,Diam0/2);
%fprintf('\nRectangle %f x %f has\n max capacity: %i circles radius %f\n',W-Diam0,H-Diam0,Nmax,Diam0/2);
if Nmax<Ap 
%     tic;
    fprintf('\nCannot fit in more than %i circles.\n',Nmax);
%     X=0;Y=0;Nmax=0;Dmatrix=0;delay=toc;
    return; 
end                                                         % error message in case requested amount points above Nmax
 
As=21;                                                   % amount sides polygon approximating circles
a=linspace(0,2*pi,As);                       % angle for  circles

dx=1;
dy=1;
refine=0;
X=zeros(1,Ap);Y=zeros(1,Ap);
if W*H>=900^2
dx=1;dy=1;
refine=1;
else
    dx=.1;dy=.1;
    refine=0;
end

x_grid=[(-W+Diam0)/2:dx:(W-Diam0)/2];y_grid=[(-H+Diam0)/2:dy:(H-Diam0)/2];       % avoid circles hitting frame
[X_grid,Y_grid]=meshgrid(x_grid,y_grid);
P=[X_grid(:)';Y_grid(:)'];
[sz1P,sz2P]=size(P);

if refine==0
xc2_base=Diam0/2*cos(a);yc2_base=Diam0/2*sin(a); 
elseif refine==1
    xc2_base=(Diam0/2+1)*cos(a);yc2_base=(Diam0/2+1)*sin(a); 
end

% figure;ax=gca;ax.DataAspectRatio=[1 1 1];
% ax.XLim=[-W/2 W/2];ax.YLim=[-H/2 H/2];
% ax.XTick=[(-W+Diam0)/2:10:(W-Diam0)/2];ax.YTick=[(-H+Diam0)/2:10:(H-Diam0)/2];grid on;hold all;
%perimeter=[(-W+Diam0)/2-.25 (W-Diam0)/2+.25 (W-Diam0)/2+.25 (-W+Diam0)/2-.25 (-W+Diam0)/2-.25;
                     %(-H+Diam0)/2-.25 (-H+Diam0)/2-.25 (H-Diam0)/2+.25 (H-Diam0)/2+.25 (-H+Diam0)/2-.25];
%plot(perimeter(1,:),perimeter(2,:),'Color',[.3 .3 .3]);
% tic;
for k=1:1:Ap
   [sz1P,sz2P]=size(P);
   if sz2P>0
       nP = randi(sz2P,1,1);
       else
            break
   end
   
   if refine==0
       X(k)=P(1,nP);Y(k)=P(2,nP);
       elseif refine==1
            dec_xP=randi([0 499],1,1);dec_xP=dec_xP/1e3;  % worst case added jitter bringing 2 points on crash course
            dec_yP=randi([0 499],1,1);dec_yP=dec_yP/1e3;
            X(k)=P(1,nP)+dec_xP;Y(k)=P(2,nP);Y(k)=P(2,nP)+dec_yP;
   end
   
   xc2=xc2_base+X(k);yc2=yc2_base+Y(k);
   in=inpolygon(P(1,:),P(2,:),xc2,yc2);
   P(:,find(in>0))=[];                                                                          % exclude area already busy
%     figure(1);plot(X(k),Y(k),'r*');                                                         % centre circles
%     figure(1);plot(xc2,yc2,'Color',[0.8 0.8 1]);                                % circles radius Diam0/2
end

L2=combinator(Ap,2);  % test 2
Dmatrix=reshape(((X(L2(:,2))-X(L2(:,1))).^2+(Y(L2(:,2))-Y(L2(:,1))).^2).^.5,[Ap Ap]);
Dmatrix([1:Ap+1:end])=NaN ;
%hold off
% delay=toc;
Loc = [X(:) Y(:)];
end


function Anxny=calc_amount_circles(H_,L_,D_)
% calculates 1. amount of circles in hex pattern that fit within 2D rectangle L (columns) x H (tall, lines or rows)
% used graph from http://www.engineeringtoolbox.com/circles-within-rectangle-d_1905.html to calibrate

R_=D_/2;

if L_<(2*R_) 
    Nx=0; 
else
    Nx=floor(L_/(2*R_));
end
    
if H_<(2*R_) 
    Ny=0;
else
    s=1;
    while H_/(2*R_+s*R_*3^.5)>1
        s=s+1; 
    end
    Ny=s;
end

if rem(L_,2*R_)>=R_
   min1_evenlines=0;
   else min1_evenlines=1; 
end

Anxny=Nx*Ny-floor(Ny/2)*min1_evenlines;
end

% 
% function handling_input_errors
% % input checks
% narginchk(5,5);nargoutchk(4,4);
% 
% err_message={'error input type W';'error input type H';'error input type R0';'error input type R0';'error input saturate'};
% 
% if(~isreal(H) || ~isscalar(W) || W<=0 )    
%     error(err_message{1}); 
% end
% if(~isreal(H) || ~isscalar(H) || H<=0 )    
%     error(err_message{2}); 
% end
% if(~isreal(R0) || ~isscalar(R0) || R0<=0 )    
%     error(err_message{3});
% end
% 
% if(~isreal(Ap) || ~isscalar(Ap) || Ap<=0 )    
%     error(err_message{3}); 
% end
% if(~isreal(saturate) || ~isscalar(saturate) || saturate<0 || saturate>1)    
%     error(err_message{3}); 
% end
% end


function [A] = combinator(N,K,s1,s2)
%COMBINATOR  Perform basic permutation and combination samplings.
% COMBINATOR will return one of 4 different samplings on the set 1:N,  
% taken K at a time.  These samplings are given as follows:
%    
% PERMUTATIONS WITH REPETITION/REPLACEMENT
%   COMBINATOR(N,K,'p','r')  --  N >= 1, K >= 0
% PERMUTATIONS WITHOUT REPETITION/REPLACEMENT
%   COMBINATOR(N,K,'p')  --  N >= 1, N >= K >= 0
% COMBINATIONS WITH REPETITION/REPLACEMENT
%   COMBINATOR(N,K,'c','r')  --  N >= 1, K >= 0
% COMBINATIONS WITHOUT REPETITION/REPLACEMENT
%   COMBINATOR(N,K,'c')  --  N >= 1, N >= K >= 0
%
% Example:
%
% To see the subset relationships, do this:  
%     combinator(4,2,'p','r')  % Permutations with repetition
%     combinator(4,2,'p')      % Permutations without repetition
%     combinator(4,2,'c','r')  % Combinations with repetition
%     combinator(4,2,'c')      % Combinations without repetition
%
%
% If it is desired to use a set other than 1:N, simply use the output from 
% COMBINATOR as an index into the set of interest.  For example:
% 
%    MySet = ['a' 'b' 'c' 'd'];
%    MySetperms = combinator(length(MySet),3,'p','r'); % Take 3 at a time.
%    MySetperms = MySet(MySetperms)
%    
%   
%    Class support for input N:
%       float: double, single 
%       integers: int8,int16,int32
%
%
% Notes: 
% All of these algorithms have the potential to create VERY large outputs.
% In each subfunction there is an anonymous function which can be used to
% calculate the number of row which will appear in the output.  If a rather
% large output is expected, consider using an integer class to conserve
% memory.  For example: 
%
%          M = combinator(int8(30),3,'p','r');  % NOT uint8(30)
%
% will take up 1/8 the memory as passing the 30 as a double.  See the note
% below on using the MEX-File.
%
% To make your own code easier to read, the fourth argument can be any 
% string.  If the string begins with an 'r' (or 'R'), the function
% will be called with the replacement/repetition algorithm.  If not, the
% string will be ignored.  
% For instance, you could use:  'No replacement', or 'Repetition allowed'
% If only two inputs are used, the function will assume 'p','r'.
% The third argument must begin with either a 'p' or a 'c' but can be any
% string beyond that.
%
% The permutations with repetitions algorithm uses cumsum.  So does the
% combinations without repetition algorithm for the special case of K=2.
% Unfortunately, MATLAB does not allow cumsum to work with integer classes.
% Thus a subfunction has been placed at the end for the case when these
% classes are passed.  The subfunction will automatically pass the
% necessary matrix to the built-in cumsum when a single or double is used.
% When an integer class is used, the subfunction first looks to see if the
% accompanying MEX-File (cumsumall.cpp) has been compiled.  If not, 
% then a MATLAB For loop is used to perform the cumsumming.  This is 
% VERY slow!  Therefore it is recommended to compile the MEX-File when 
% using integer classes. 
% The MEX-File was tested by the author using the Borland 5.5 C++ compiler.
% 
% See also, perms, nchoosek, npermutek (on the FEX)
%             
% Author:   Matt Fig
% Contact:  popkenai@yahoo.com
% Date:     5/30/2009
%
% Reference:  http://mathworld.wolfram.com/BallPicking.html

% N=vpa(N,27)
% K=vpa(K,27)

ng = nargin;

if ng == 2
    s1 = 'p';
    s2 = 'r';
elseif ng == 3 
    s2 = 'n';
elseif ng ~= 4
    error('Only 2, 3 or 4 inputs are allowed.  See help.')
end

if isempty(N) || K == 0
   A = [];  
   return
elseif numel(N)~=1 || N<=0 || ~isreal(N) || floor(N) ~= N 
    error('N should be one real, positive integer. See help.')
elseif numel(K)~=1 || K<0 || ~isreal(K) || floor(K) ~= K
    error('K should be one real non-negative integer. See help.')
end

STR = lower(s1(1)); % We are only interested in the first letter.

if ~strcmpi(s2(1),'r')
    STR = [STR,'n'];
else
   STR = [STR,'r']; 
end


switch STR
    case 'pr'
        A = perms_rep(N,K);     % strings
    case 'pn'
        A = perms_no_rep(N,K);  % permutations
    case 'cr'
        A = combs_rep(N,K);     % multichoose
    case 'cn'
        A = combs_no_rep(N,K);  % choose
    otherwise
        error('Unknown option passed.  See help')
end

end



function PR = perms_rep(N,K)
% This is (basically) the same as npermutek found on the FEX.  It is the  
% fastest way to calculate these (in MATLAB) that I know.  
% pr = @(N,K) N^K;  Number of rows.
% A speed comparison could be made with COMBN.m, found on the FEX.  This
% is an excellent code which uses ndgrid.  COMBN is written by Jos.
%
%      % All timings represent the best of 4 consecutive runs.
%      % All timings shown in subfunction notes used this configuration:
%      % 2007a 64-bit, Intel Xeon, win xp 64, 16 GB RAM  
%      tic,Tc = combinator(single(9),7,'p','r');toc  
%      %Elapsed time is 0.199397 seconds.  Allow Ctrl+T+C+R on block
%      tic,Tj = combn(single(1:9),7);toc  
%      %Elapsed time is 0.934780 seconds.
%      isequal(Tc,Tj)  % Yes

if N==1
   PR = ones(1,K,class(N)); 
   return
elseif K==1
    PR = (1:N).';
    return
end

CN = class(N);
M = double(N);  % Single will give us trouble on indexing.
L = M^K;  % This is the number of rows the outputs will have.
PR = zeros(L,K,CN);  % Preallocation.
D = ones(1,N-1,CN);  % Use this for cumsumming later.
LD = M-1;  % See comment on N. 
VL = [-(N-1) D].';  % These values will be put into PR.
% Now start building the matrix.
TMP = VL(:,ones(L/M,1,CN));  % Instead of repmatting.
PR(:,K) = TMP(:);  % We don't need to do two these in loop.
PR(1:M^(K-1):L,1) = VL;  % The first column is the simplest.
% Here we have to build the cols of PR the rest of the way.
for ii = K-1:-1:2
    ROWS = 1:M^(ii-1):L;  % Indices into the rows for this col.
    TMP = VL(:,ones(length(ROWS)/(LD+1),1,CN));  % Match dimension.
    PR(ROWS,K-ii+1) = TMP(:);  % Build it up, insert values.
end

PR(1,:) = 1;  % For proper cumsumming.
PR = cumsum2(PR);  % This is the time hog.
end



function PN = perms_no_rep(N,K)
% Subfunction: permutations without replacement.
% Uses the algorithm in combs_no_rep as a basis, then permutes each row.
% pn = @(N,K) prod(1:N)/(prod(1:(N-K)));  Number of rows.

if N==K
    PN = perms_loop(N);  % Call helper function.
%     [id,id] = sort(PN(:,1));  %#ok  Not nec., uncomment for nice order.
%     PN = PN(id,:);  % Return values.
    return
elseif K==1
    PN = (1:N).';  % Easy case.
    return
end

if K>N  % Since there is no replacement, this cannot happen.
    error(['When no repetitions are allowed, '...
           'K must be less than or equal to N'])
end

M = double(N);  % Single will give us trouble on indexing.
WV = 1:K;  % Working vector.
lim = K;   % Sets the limit for working index.
inc = 1;   % Controls which element of WV is being worked on.
BC = prod(M-K+1:M);  % Pre-allocation of return arg.
BC1 = BC / ( prod(1:K)); % Number of comb blocks.
PN = zeros(round(BC),K,class(N));
L = prod(1:K) ;  % To get the size of the blocks.
cnt = 1+L;
P = perms_loop(K);  % Only need to use this once.
PN(1:(1+L-1),:) = WV(P);  % The first row.

for ii = 2:(BC1 - 1);
    if logical((inc+lim)-N)  % The logical is nec. for class single(?)
        stp = inc;  % This is where the for loop below stops.
        flg = 0;  % Used for resetting inc.
    else
        stp = 1;
        flg = 1;
    end
    
    for jj = 1:stp
        WV(K  + jj - inc) = lim + jj;  % Faster than a vector assignment!
    end
                                                                             
    PN(cnt:(cnt+L-1),:) = WV(P);  % Assign block.
    cnt = cnt + L;  % Increment base index.    
    inc = inc*flg + 1;  % Increment the counter.
    lim = WV(K - inc + 1 );  % lim for next run.
end

V = (N-K+1):N;  % Final vector.
PN(cnt:(cnt+L-1),:) = V(P);  % Fill final block.
% The sorting below is NOT necessary.  If you prefer this nice
% order, the next two lines can be un-commented.
% [id,id] = sort(PN(:,1));  %#ok  This is not necessary!
% PN = PN(id,:);  % Return values.
end



function P = perms_loop(N)
% Helper function to perms_no_rep.  This is basically the same as the
% MATLAB function perms.  It has been un-recursed for a runtime of around  
% half the recursive version found in perms.m  For example:
%
%      tic,Tp = perms(1:9);toc
%      %Elapsed time is 0.222111 seconds.  Allow Ctrl+T+C+R on block
%      tic,Tc = combinator(9,9,'p');toc  
%      %Elapsed time is 0.143219 seconds.
%      isequal(Tc,Tp)  % Yes

M = double(N); % Single will give us trouble on indexing.
P = 1;  % Initializer.
G = cumprod(1:(M-1));  % Holds the sizes of P.
CN = class(N);

for n = 2:M
    q = P;
    m = G(n-1);
    P = zeros(n*m,n,CN);
    P(1:m, 1) = n;
    P(1:m, 2:n) = q;
    a = m + 1;
    
    for ii = n-1:-1:1,
        t = q;
        t(t == ii) = n;
        b = a + m - 1;
        P(a:b, 1) = ii;
        P(a:b, 2:n) = t;
        a = b + 1;
    end 
end
end




function CR = combs_rep(N,K)
% Subfunction multichoose:  combinations with replacement.
% cr = @(N,K) prod((N):(N+K-1))/(prod(1:K)); Number of rows.

M = double(N);  % Single will give us trouble on indexing.
WV = ones(1,K,class(N));  % This is the working vector.
mch = prod((M:(M+K-1)) ./ (1:K));  % Pre-allocation.
CR = ones(round(mch),K,class(N));

for ii = 2:mch
    if WV(K) == N
        cnt = K-1;  % Work backwards in WV.
        
        while WV(cnt) == N
            cnt = cnt-1;  % Work backwards in WV.
        end

        WV(cnt:K) = WV(cnt) + 1;  % Fill forward.
    else
        WV(K) = WV(K)+1;   % Keep working in this group.
    end

    CR(ii,:) = WV;
end
end



function CN = combs_no_rep(N,K)
% Subfunction choose:  combinations w/o replacement.
% cn = @(N,K) prod(N-K+1:N)/(prod(1:K));  Number of rows.
% Same output as the MATLAB function nchoosek(1:N,K), but often faster for
% larger N.
% For example: 
%
%      tic,Tn = nchoosek(1:17,8);toc
%      %Elapsed time is 0.430216 seconds.  Allow Ctrl+T+C+R on block
%      tic,Tc = combinator(17,8,'c');toc  
%      %Elapsed time is 0.024438 seconds.
%      isequal(Tc,Tn)  % Yes

if K>N
    error(['When no repetitions are allowed, '...
           'K must be less than or equal to N'])
end

M = double(N);  % Single will give us trouble on indexing.

if K == 1
   CN =(1:N).';  % These are simple cases.
   return
elseif K == N
    CN = (1:N);
    return
elseif K==2 && N>2  % This is an easy case to do quickly.
    BC = (M-1)*M / 2;
    id1 = cumsum2((M-1):-1:2)+1;
    CN = zeros(BC,2,class(N));
    CN(:,2) = 1;
    CN(1,:) = [1 2];
    CN(id1,1) = 1;
    CN(id1,2) = -((N-3):-1:0);
    CN = cumsum2(CN);
    return
end

WV = 1:K;  % Working vector.
lim = K;   % Sets the limit for working index.
inc = 1;   % Controls which element of WV is being worked on.
BC = prod(M-K+1:M) / (prod(1:K));  % Pre-allocation.
CN = zeros(round(BC),K,class(N));
CN(1,:) = WV;  % The first row.

for ii = 2:(BC - 1);   
    if logical((inc+lim)-N) % The logical is nec. for class single(?)
        stp = inc;  % This is where the for loop below stops.
        flg = 0;  % Used for resetting inc.
    else
        stp = 1;
        flg = 1;
    end
    
    for jj = 1:stp
        WV(K  + jj - inc) = lim + jj;  % Faster than a vector assignment.
    end
    
    CN(ii,:) = WV;  % Make assignment.
    inc = inc*flg + 1;  % Increment the counter.
    lim = WV(K - inc + 1 );  % lim for next run. 
end
  
CN(ii+1,:) = (N-K+1):N;
end


function A = cumsum2(A)
%CUMSUM2, works with integer classes. 
% Duplicates the action of cumsum, but for integer classes.
% If Matlab ever allows cumsum to work for integer classes, we can remove 
% this.

if isfloat(A)
    A = cumsum(A);  % For single and double, use built-in.
    return
else 
    for ii = 2:size(A,1)
        A(ii,:) = A(ii,:) + A(ii-1,:); % User likes it slow.
    end  
end
end



