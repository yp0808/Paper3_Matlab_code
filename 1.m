s_n_max = (wmax+1)^2*(1/kappa+1)^2;
wmax = 10;
kappa = 0.05;

%将状态编号化
s=1;
s_table = cell(1000,2);
while s <= 1000;
    w1 = fix((s-1)/((wmax+1)*(1/kappa+1)^2));
    w2 = fix(mod(s-1,(wmax+1)*(1/kappa+1)^2)/(1/kappa+1)^2);
    sigma_s = fix((s-1-w1*(wmax+1)*(1/kappa+1)^2-w2*(1/kappa+1)^2)/(1/kappa+1));
    sigma_b = s-1-w1*(wmax+1)*(1/kappa+1)^2-w2*(1/kappa+1)^2-sigma_s*(1/kappa+1);
    s_tuple = [w1,w2,sigma_s,sigma_b];
    s_table{s,1} = s;
    s_table{s,2} = s_tuple;
    s = s+1;
end

tol = 0.1;  %收敛误差
newvalue = []; newx = []; newp = []; 
oldvalue = []; oldx = []; oldp = [];
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
       
        oldx = newx; oldp = nwep; oldvalue = newvalue;
    end
%-----------------------------    
function [] = contract()
%更新决策变量和值函数
global newvalue newx oldvalue oldx s_max

w = 1;
while s <= s_n_max;
    [newx(s), newp(s)newvalue(s)] = optimize(s);
    w=w+1;
end


%--------------------------------------
function [out1,out2] = optimize(w)
% This procedure calculates optimal investment, and value fn., for a
% given industry structure s.

global a beta entry_k isentry nfirms oldvalue oldx phi profit
locw = qdecode(w);
locwx = locw;
oval = oldvalue(w,:)';
ox = oldx(w,:)';
nval = zeros(nfirms,1);
nx = zeros(nfirms,1);

% Find out which firms want to exit

i = (min(oval) == phi)*(minind(oval)-1) + (min(oval) > phi)*nfirms;
%有误,但在本例中，oval一定会>=phi。
%i = nfirms-sum(oval <= phi)会更好。

% Replace efficiency levels of exitors with zero

if i < nfirms;
    locwx(i+1:nfirms) = zeros(nfirms-i,1);
    %locwx：退出后的状态
end

% Figure out the probability of entry

entered = isentry(qencode(flipud(sortrows(flipud(locwx),1)))); %由大到小重新排列
locwe = locwx;
locwe(nfirms) = entry_k;%退出和进入后的状态

% Now calculate the optimal policies for this industry structure, given that
% entry and exit are as specified.
 
j=1;
while j <= nfirms;
    if locw(j) == 0;
        nval(j:nfirms) = phi*ones(nfirms-j+1,1);
        break;
    end
    v1=0; v2=0;
    if entered < 1;
        
        % First: Calculate v, without entry
        
        [v1, v2] = calcval(j,locwx,ox,locw(j));
    end
    
    if entered > 0;
        
        % A firm wants to enter with positive probability
        
        [tempv1, tempv2] = calcval(j,locwe,ox,locw(j));
        v1 = entered*tempv1 + (1-entered)*v1;
        v2 = entered*tempv2 + (1-entered)*v2;
    end
    
    % Calculate values for firm, given that it is not leaving
    
    if v1 <= v2; % Avoid division by zeros
        r = 1.0;%代表p = 0
    else; r = 1.0/(beta*a*(v1-v2));
    end
    
    % r now contains the value r = (1 - p)^2. => p = 1 - sqrt(r)),
    % where p is the optimal prob. of having k rise, cond. on world
    
    r = min([max([r;0.0000000000001]);1]);
    p = 1.0 - sqrt(r);
    nx(j) = p/(a - a * p);
    
    % Now calculate the value from staying in
    % Ask: given this optimal investment level, will there be exit?
    
    nval(j) = profit(w,j) - nx(j) + beta*(v1*p + v2*(1-p));
    if nval(j) <= phi;
        nval(j) = phi;
        nx(j) = 0;
    end
    if (j < nfirms) & (nval(j) == phi);
        nval(j+1:nfirms) = ones(nfirms-j,1) * phi;
        break;
    end
    ox(j) = nx(j);
    locwx(j) = (nval(j) > phi)*locw(j);
    locwe(j) = locwx(j);
    j=j+1;
end

out1 = nx';
out2 = nval';

