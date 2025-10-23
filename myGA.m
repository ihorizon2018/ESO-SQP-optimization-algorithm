function [bestGene_Position,best_Gene,bestfitGene] = myGA(ObjFun,LB,UB,nVar,nPop,MaxIt)
% Simple real-coded Genetic Algorithm (minimization)
% Interface: myGA(ObjFun, LB, UB, nVar, nPop, MaxIt)
% Returns:
%   bestGene_Position  - best solution found (1 x nVar)
%   best_Gene          - best objective value
%   bestfitGene        - best-so-far fitness per iteration (1 x MaxIt)

% ---------- standardize inputs ----------
dim = nVar;
N   = nPop;

% row-vector bounds of length dim
if isscalar(LB), lbv = LB*ones(1,dim); else, lbv = LB(:)'; end
if isscalar(UB), ubv = UB*ones(1,dim); else, ubv = UB(:)'; end
span = ubv - lbv;

% ---------- GA hyperparameters ----------
pc   = 0.9;             % crossover probability
pm   = 1/dim;           % per-gene mutation probability
eliteRate = 0.05;       % elitism ratio
eliteCount = max(1, round(eliteRate*N));
sigma = 0.1*span;       % mutation std (per gene)

% ---------- initialization ----------
P = repmat(lbv,N,1) + rand(N,dim).*repmat(span,N,1);  % population
f = arrayfun(@(i) ObjFun(P(i,:)), (1:N)');            % fitness (column)

% track best
[f_gbest, idx] = min(f);
gbest = P(idx,:);

bestfitGene = nan(1, MaxIt);

% ---------- main loop ----------
for t = 1:MaxIt
    % sort by fitness (ascending)
    [f, order] = sort(f);
    P = P(order,:);

    % elitism pool
    elites  = P(1:eliteCount, :);
    fElite  = f(1:eliteCount);

    % offspring container
    offCount = N - eliteCount;
    Off = zeros(offCount, dim);
    fOff = zeros(offCount, 1);

    % offspring generation (tournament selection + uniform crossover + Gaussian mutation)
    k = 1;
    while k <= offCount
        % --- tournament selection (size 2) ---
        p1 = tournamentSelect(f, 2);
        p2 = tournamentSelect(f, 2);
        x1 = P(p1,:); x2 = P(p2,:);

        % --- crossover ---
        if rand < pc
            mask = rand(1,dim) < 0.5;
            c1 = x1; c2 = x2;
            c1(mask) = x2(mask);
            c2(mask) = x1(mask);
        else
            c1 = x1; c2 = x2;
        end

        % --- mutation (Gaussian) ---
        mutMask1 = rand(1,dim) < pm;
        mutMask2 = rand(1,dim) < pm;
        if any(mutMask1)
            c1(mutMask1) = c1(mutMask1) + sigma(mutMask1).*randn(1, sum(mutMask1));
        end
        if any(mutMask2)
            c2(mutMask2) = c2(mutMask2) + sigma(mutMask2).*randn(1, sum(mutMask2));
        end

        % --- bounds handling ---
        c1 = min(max(c1, lbv), ubv);
        c2 = min(max(c2, lbv), ubv);

        % assign to Off (ensure not exceeding needed count)
        Off(k,:) = c1;  fOff(k) = ObjFun(c1);
        if k+1 <= offCount
            Off(k+1,:) = c2; fOff(k+1) = ObjFun(c2);
        end
        k = k + 2;
    end

    % next generation: elites + offspring
    P = [elites; Off];
    f = [fElite; fOff];

    % update global best
    [f_curBest, idx] = min(f);
    if f_curBest < f_gbest
        f_gbest = f_curBest;
        gbest   = P(idx, :);
    end

    % record curve (best-so-far)
    bestfitGene(t) = f_gbest;
end

% ---------- outputs ----------
best_Gene          = f_gbest;
bestGene_Position  = gbest;
bestfitGene        = cummin(bestfitGene(:))';   % ensure non-increasing and row vector
end


% ===== local helper =====
function idx = tournamentSelect(f, tourSize)
% Select index by tournament of size tourSize (minimization)
n = numel(f);
cands = randi(n, [tourSize, 1]);
[~, k] = min(f(cands));
idx = cands(k);
end
