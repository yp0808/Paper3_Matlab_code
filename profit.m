% oct2mat
% gauss2oct


function [] = profit()




global c binom p mc sigma wstar M
global ggamma D f
global profit agprof csurplus share pmcmarg concent

nfmax = c.MAX_FIRMS; % max # of active firms最大的活跃公司数量
kkmax = c.KMAX;      % max efficiency level attainable能达到的最高效率水平
it = c.IND_TYPE;     % investment type (quality/mc/capacity)投资类型（质量/mc/产量）
et = c.EQL_TYPE;     % equilibrium type (Nash/monopoly/social planner)均衡类型

% "*** Computing profit function ***";




nfirms=1;
while nfirms <= nfmax;
    %     Number of descending n-tuples
    
    disp(sprintf('\nFirms: %d', nfirms))
    descn = binom(nfirms+kmax+1,kmax+2);
    disp(sprintf('Industry structures to compute: %d', descn))
    if strcmp(et, 'COMPETITION') nagents = nfirms;%nagent决策主体
    else nagents = 1; end
    
    profit = zeros(descn,nagents);
    % The following variables are used for comparative statics
    agprof = zeros(descn,nfirms); % Aggregate profits per industry structure
    csurplus = zeros(descn,1);    % Consumer surplus
    share = zeros(descn,nfirms);  % Market share of each firm
    pmcmarg = zeros(descn,1);     % Price/mc margin, average by sales
    concent = zeros(descn,1);     % One-firm concentration ratio
    
    
    %      now, call appropriate profit function
    
    if strcmp(it, 'QUALITY');
        %      mc = c.MC;
        %      M = c.MKT_SIZE;
        %      wstar = c.WSTAR;
        %      w = []; egw = []; egwp = []; p = []; profstar = [];
        %      sigma=zeros(nfirms,1);
        %      p=5.5*ones(nfirms,1);
        %      if strcmp(et, 'COMPETITION');
        %         cqprofit(nfirms, descn);
        %      elseif strcmp(et, 'MONOPOLY');
        %         mqprofit(nfirms,descn);
        %      elseif strcmp(et, 'PLANNER');
        %         sqprofit(nfirms,descn);
        %      end
    elseif strcmp(it, 'COST');
        D = c.INTERCEPT;
        f = c.FIXED_COST;
        ggamma = c.GAMMA;
        quan = []; profstar = []; w = []; theta = []; pstar = [];
        if strcmp(et, 'COMPETITION');
            ccprofit(nfirms, descn);
        elseif strcmp(et, 'MONOPOLY');
            mcprofit(nfirms,descn);
        elseif strcmp(et, 'PLANNER');
            scprofit(nfirms,descn);
        end
    elseif strcmp(it, 'CAPACITY');
        %      D = c.INTERCEPT;
        %      mc = c.MC;
        %      tau = c.TAU;
        %      quan = []; profstar = []; w = []; pstar = [];
        %      if strcmp(et, 'COMPETITION');
        %         cpprofit(nfirms, descn);
        %      elseif strcmp(et, 'MONOPOLY');
        %         mpprofit(nfirms,descn);
        %      elseif strcmp(et, 'PLANNER');
        %         spprofit(nfirms,descn);
        %      end
    end
    
    %   write output
    %     "Generating output --> ";;file1;;", ";;file2;
    
    s = int2str(nfirms);
    save(['a.' c.PREFIX '_pr' s '.mat'], 'profit');
    save(['a.' c.PREFIX '_cons' s '.mat'], ...
        'agprof', 'csurplus', 'share', 'pmcmarg', 'concent')
    
    nfirms=nfirms+1;
end
c.PROFIT_DONE=1;



function [] = progress(i)
% report progress
if mod(i, 50) == 0;
    % "Computed: ";;compact(i);;"\r";;
    disp(sprintf('  Computed: %d', i));
end


function [] = ccprofit(nfirms,descn)
% competition
% local i, n, q, p;

global ggamma D f
global profit agprof csurplus share pmcmarg concent
i = 1;
while i <= descn;
    progress(i);
    w = cdecode(i,nfirms+1);
    theta = ggamma * exp(-w);  % marginal cost
    
    %   quan = solveeq(theta); inlined
    %   Solve for equilibrium with n firms; reduce n until all firms
    %   want to produce quantity > 0
    
    n=nfirms;
    p = (D + sum(theta(1:n)))/(n+1);
    while ~((p - theta(n) >= 0) | (n==1)); %两个条件同时不成立
        n=n-1;
        p = (D + sum(theta(1:n)))/(n+1);
    end
    q = zeros(nfirms,1);
    if p - theta(n) > 0;
        q(1:n) = p - theta(1:n);
    end
    quan=q;%quantity
    
    pstar = D - sum(quan);   % Equilibrium price
    profstar = (pstar>theta).*(pstar-theta).*quan - f; % Equilibrium profits
    
    profit(i,:) = profstar';
    csurplus(i) = 0.5*sum(quan)*sum(quan);
    agprof(i,:) = profstar';
    share(i,:) = quan';
    if sum(quan) > 0;
        pmcmarg(i) = pstar / (sum(theta.*quan)) * sum(quan);
        concent(i) = max(quan)/sum(quan);
    else;
        pmcmarg(i) = 1;
    end
    i = i+1;
end



function [] = mcprofit(nfirms,descn)
% monopolist
% local i, numin;

global ggamma D f
global profit agprof csurplus share pmcmarg concent
i = 1;
while i <= descn;
    progress(i);
    
    %   The monopolist will always choose to produce everything from the lowest
    %   priced firm, so it acts like a 1-plant firm for the static profits.
    
    w = cdecode(i,nfirms+1);
    numin = sum(w > -4); % No. of firms in, for fixed-fee computation
    agprof(i,:) = -f * ((w > -4)');
    w = max(([(w-1),(-4*ones(nfirms,1))])');
    
    %   As zero here represents being out, move everything down by one, except for
    %   zero, as there is no negative efficiency level
    
    theta = ggamma * exp(-w(1));  % marginal cost
    pstar = 0.5*(D + theta);   % One-plant monopolist price
    quan = (pstar>theta)*(pstar-theta);
    profstar = quan*(pstar-theta) - f*numin; % Monopolist profits
    
    profit(i) = profstar;
    agprof(i,1) = agprof(i,1) + profstar;
    csurplus(i) = 0.5*quan*quan;
    pmcmarg(i) = pstar / theta;
    concent(i) = 1;
    if nfirms > 1;
        quan = [quan;zeros(nfirms-1,1)];
    end
    share(i,:) = quan';
    i = i+1;
end



function [] = scprofit(nfirms, descn)
% local i, numin;

global ggamma D f
global profit agprof csurplus share pmcmarg concent
i = 1;
while i <= descn;
    progress(i);
    
    %   The social planner will always choose to produce everything from the lowest
    %   priced firm, so it acts like a 1-plant firm for the static profits.
    
    w = cdecode(i,nfirms+1);
    agprof(i,:) = -f*((w>-4)');
    numin = sum(w > -4); % No. of firms in, for fixed-fee computation
    w = max(([(w-1),(-4*ones(nfirms,1))])');
    
    %   As zero here represents being out, move everything down by one, except for
    %   zero, as there is no negative efficiency level
    
    theta = ggamma * exp(-w(1));  % marginal cost
    pstar = theta;  % Set price = mc, for social planner solution
    quan = (D>theta)*(D-pstar);
    profstar = 0.5*(D-theta)*quan; % Consumer surplus
    
    profit(i) = profstar - f*numin;  % Producer surplus
    csurplus(i) = profstar;
    pmcmarg(i) = 1;
    concent(i) = 1;
    if nfirms > 1;
        quan = [quan;zeros(nfirms-1,1)];
    end
    share(i,:) = quan';
    i = i+1;
end



function [out1] = cdecode(code,nfirms)
% Cournot
% This procedure takes a previously encoded number, and decodes it into
% a weakly descending n-tuple (n = nfirms - 1)
% local ntuple,digit,i;

global binom
code = code-1;
ntuple = zeros(nfirms-1,1);
i = 1;
while i <= nfirms - 1;
    digit = 0;
    while binom(digit+nfirms-i+1,digit+2) <= code;
        digit=digit+1;
    end
    ntuple(i) = digit;
    code = code-binom(digit+nfirms-i,digit+1);
    i = i+1;
end

% Now convert to format of starting at -4, and jumping by 1's

ntuple = ntuple-4;
out1 = ntuple;


function [out1] = pdecode(code,nfirms)
% Capacity
% This procedure takes a previously encoded number, and decodes it into
% a weakly descending n-tuple (n = nfirms - 1)
% local ntuple,digit,i;

global binom
code = code-1;
ntuple = zeros(nfirms-1,1);
i = 1;
while i <= nfirms - 1;
    digit = 0;
    while binom(digit+nfirms-i+1,digit+2) <= code;
        digit=digit+1;
    end
    ntuple(i) = digit;
    code = code-binom(digit+nfirms-i,digit+1);
    i = i+1;
end
out1 = ntuple;
