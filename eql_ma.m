



function [] = eql_ma()

global c wmax beta delta a
global newvalue newx newp oldvalue oldx oldp smax 

wmax = c.WMAX;
beta = c.BETA;
delta = c.DELTA;
a = c.INV_MULT;

tol = 0.1;  % Tolerance for convergence
newvalue = []; newx = []; newp = []; 
oldvalue = []; oldx = []; oldp = [];


    % Update values, or define starting values.
    
    update;
    disp(sprintf('Contraction ...'));
    ix = 1;
    norm = tol + 1;
    avgnorm = norm;
    while (norm > tol) & (avgnorm > 0.0001*tol);
        contract;
        norm = max(max(abs(oldvalue - newvalue)));
        avgnorm = mean(mean(abs(oldvalue-newvalue)));
        
        disp(sprintf('  %2d    Sup norm: %8.4f      Mean norm: %8.4f', ...
            ix, norm, avgnorm));
        ix = ix+1;
        
        iterdiff = abs(oldvalue-newvalue);
        normind1 = maxind(max(iterdiff'));
        normind2 = maxind(max(iterdiff));
        normcode = (qdecode(normind1))';
        
        % "Max. elt is: " normind2 "," normcode "; Old value: "
        % oldvalue(normind1,normind2) "; New value: "
        % newvalue(normind1,normind2) "";
        
        oldx = newx; oldp = nwep; oldvalue = newvalue;
    end
    
    % d2 = date;
    % Now find if there is any investment at the highest level.
    
    w=kmax;
    if nfirms > 1;
        w = [w;zeros(nfirms-1,1)];
    end
    if max(newx(qencode(w):wmax,1)) > 0;
        disp('Warning: Positive investment recorded at highest efficiency level.')
        disp('Please consider increasing the maximum efficiency level (kmax).')
    end
    
    % Store data in file for inspection
    % Store data in file for comparative statics program to read
    

    save(['a.' c.PREFIX '_markov' '.mat'], ...
        'newvalue', 'newx')
    
    % disp(sprintf('\n'))
    % disp('Value Function (wmax x nfirms)')
    % disp([dtable' newvalue])
    % disp('Investment (wmax x nfirms)')
    % disp([dtable' newx])
    % disp('Probability of p rising (wmax x nfirms)'),
    % disp([dtable' prising])
    % disp('Probability of entry (wmax x nfirms)')
    % disp([dtable' isentry])
      
end

c.EQL_DONE = 1;



function [] = contract()
% This procedure does one iterative step on investment and the value fn
% Implicit parameters are oldx, oldvalue (passed in), and newx, newvalue,
%  which are returned

global newvalue newx  newp oldvalue oldx old p smax
s = 1;
while s <= smax;
    [newx(s,:), newp(s,:), newvalue(s,:)] = optimize(w);
    w=w+1;
end

% Implicit returned parameters: newx, newvalue



function [] = update()
% This procedure takes the solved newx, newvalue matrix for the nfirms - 1
% problem, and puts them into the nfirms matrices oldx, oldvalue, for use
% as starting values
% local w,i,n,tuple;

global isentry nfirms wmax newvalue newx oldvalue oldx
oldx = zeros(wmax,nfirms);
oldvalue = zeros(wmax,nfirms);
if nfirms == 1;
    i=1;
    while i <= wmax;
        oldvalue(i,:) = 1 + 0.1*i;
        i=i+1;
    end
else;
    w=1;
    while w <= wmax;
        tuple = qdecode(w);
        nfirms = nfirms - 1;
        n = encode(tuple(1:nfirms));
        oldx(w,1:nfirms) = newx(n,1:nfirms);
        oldvalue(w,1:nfirms) = newvalue(n,1:nfirms);
        nfirms = nfirms + 1;
        tuple(nfirms-1) = tuple(nfirms);
        tuple(nfirms) = 0;
        oldvalue(w,nfirms) = oldvalue(encode(tuple),nfirms-1);
        oldx(w,nfirms) = oldx(encode(tuple),nfirms-1);
        w=w+1;
    end
end
isentry = zeros(wmax,1);
newx = zeros(wmax,nfirms);
newvalue = zeros(wmax,nfirms);

% Implicit returned value: oldx, oldvalue



function [out1,out2,out3] = optimize(s)
% This procedure calculates optimal investment, and value fn., for a
% given industry structure w. 

global a beta oldvalue oldx phi profit
s_tuple = qdecode(s);
oval = oldvalue(s,:);
ox = oldx(s,:);
op = oldp(s,:);
nval = 0;
nx = 0;
np = 0;


% Now calculate the optimal policies for this industry structure

 



out1 = nx';
out2 = nval';





function [out1,out2] = calcval(place,w,x,k)
% This procedure calculates val = EEEV(.,.,.,.)p(.)p(.)p(.), where E
% represents sums, and this is the calculation of the 4-firm problem
% Vars: place = place of own omega, for calculating value function (v)
%       w = the vector of omegas; already decoded
%       x = the vector of investments (nfirms of them)
% Implicit parameter: oldvalue
% For efficiency reasons, it outputs the following vector:
% [ calcval(k_v+1,w,x), calcval(k_v,w,x) ]
% local i,valA,valB,d,e,probmask,z1,z2,locmask,
%   p_up,  % p_down, p of going up/down for all other firms
%   temp,
%   pl1,justone;

global a delta kmax mask nfirms oldvalue two_n
z1 = zeros(nfirms,1);
z2 = kmax*ones(nfirms,1);

% Expand mask to allow for the non-inclusion of the ith plant

if nfirms > 1;
    if place == 1;
        locmask = [zeros(1,two_n);mask];
    elseif place == nfirms;
        locmask = [mask;zeros(1,two_n)];
    else; locmask = [mask(1:place-1,:);zeros(1,two_n);mask(place:nfirms-1,:)];
    end
else; locmask = zeros(1,1);
end
x(place) = 0;
w(place) = k;
justone = zeros(nfirms,1);
justone(place) = 1;
p_up = (a .* x) ./ (1 + a .* x);
% p_down = 1 - p_up;
valA=0; valB=0;
i=1;

while i <= two_n;
    % probmask = prod(mask(:,i) .* p_up + (1 - mask(:,i)) .* p_down);
    probmask = prod(2 .* locmask(:,i) .* p_up + 1 - locmask(:,i) - p_up);
    d = w+locmask(:,i); 
    temp = flipud(sortrows(flipud([d,justone]),1));%里面的flipud多余
    d = temp(:,1);
    e = d - 1;
    
    % Check for evaluation of value fn. at -1
    e = max(([e,z1])')';%max比较每列的最大
    % Check for evaluation of value fn. at kmax+1
    d = min(([d,z2])')';
    pl1 = maxind(temp(:,2));% sum(d(1:place)>=k) + sum(d(place:nfirms)>k);
    
    valB = valB + ((1-delta)*oldvalue(qencode(d),pl1) ...
        + delta*oldvalue(qencode(e),pl1))*probmask;
    
    d = w+locmask(:,i)+justone;
    temp = flipud(sortrows(flipud([d,justone]),1));
    d = temp(:,1);
    e = d - 1;
    
    % Check for evaluation of value fn. at -1
    e = max(([e,z1])')';
    % Check for evaluation of value fn. at kmax+1
    d = min(([d,z2])')';
    pl1 = maxind(temp(:,2)); %sum(e(1:place)>=k) + sum(e(place:nfirms)>k);
    
    valA = valA + ((1-delta)*oldvalue(qencode(d),pl1) ...
        + delta*oldvalue(qencode(e),pl1))*probmask;
    i=i+1;
end

out1 = valA;
out2 = valB;


function [out1] = encode(ntuple)
% This procedure takes a weakly descending n-tuple (n = nfirms), with
% min. elt. 0, max. elt. kmax, and encodes it into an integer
% local code,digit,i;
%将n维状态映射为自然数
global binom nfirms
code = 1; % Coding is from 1 to wmax
i = 1;
while i <= nfirms;
    digit = ntuple(i);
    code = code + binom(digit+nfirms+1-i,digit+1);
    i=i+1;
end

out1 = code;


function [out1] = qencode(ntuple)
% This procedure does a quick encode of any n-tuple given in weakly
% descending order. Encoding is done using a table lookup. Each
% column of the table consists of an n-tuple; the ith column is the ith
% n-tuple to be decoded. The table is stored in the variable "etable".

global encfirm etable1 etable2 multfac1 multfac2 nfirms
if nfirms <= encfirm;
    out1 = etable1(sum(ntuple.*multfac1)+1);
else;
    out1 = etable1(sum(ntuple.*multfac1)+1) ...
        + etable2(sum(ntuple.*multfac2)+1);
end



function [out1] = qdecode(code)
% This procedure does a quick decode of a previously encoded number into
% a weakly descending n-tuple. Decoding is done using a table lookup. Each
% column of the table consists of an n-tuple; the ith column is the ith
% n-tuple to be decoded. The table is stored in the variable "dtable".

global dtable

out1 = dtable(:,code);


function [out1] = decode(code)
% This procedure takes a previously encoded number, and decodes it into
% a weakly descending n-tuple (n = nfirms)
% local ntuple,digit,i;

global binom nfirms
code = code-1;
ntuple = zeros(nfirms,1);
i = 1;
while i <= nfirms;
    digit = 0;
    while binom(digit+nfirms-i+2,digit+2) <= code;
        digit=digit+1;
    end
    ntuple(i) = digit;
    code = code-binom(digit+nfirms-i+1,digit+1);
    i = i+1;
end

out1 = ntuple;
