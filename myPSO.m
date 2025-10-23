%% PARTICLE SWARM OPTIMIZATION (PSO)
function [bestParticle_Position,best_Particle,bestfitIter] = myPSO(ObjFun,LB,UB,nVar,nPop,MaxIt)

%% Problem Parameters
prob   = ObjFun;      % fitness function
D      = nVar;        % number of variables
cPop   = nPop;        % population size
Maxiter= MaxIt;       % max no of iteration

% --- 将上下界规范为 1×D 行向量 ---
if isscalar(LB), lbv = LB*ones(1,D); else, lbv = LB(:)'; end
if isscalar(UB), ubv = UB*ones(1,D); else, ubv = UB(:)'; end
span = ubv - lbv;

%% PSO Parameters
w  = 0.45;            % inertia
c1 = 1;               % cognitive
c2 = 1;               % social

%% Initialization
% 位置 P：均匀采样到 [lbv,ubv]
P = repmat(lbv,cPop,1) + rand(cPop,D).*repmat(span,cPop,1);
% 速度 V：给个与范围同量级的初值（可选）
V = 0.2*repmat(span,cPop,1).*(rand(cPop,D)-0.5);

f = NaN(cPop,1);
for q = 1:cPop
    f(q) = prob(P(q,:));
end

pbest   = P;             % 个体最优位置
f_pbest = f;             % 个体最优值
[f_gbest,ind] = min(f_pbest);
gbest   = pbest(ind,:);  % 全局最优位置

bestfitIter = NaN(Maxiter,1);

%% Main Loop
for t = 1:Maxiter
    for q = 1:cPop
        % 速度更新
        V(q,:) = w*V(q,:) ...
               + c1*rand(1,D).*(pbest(q,:)-P(q,:)) ...
               + c2*rand(1,D).*(gbest      -P(q,:));
        % 位置更新
        P(q,:) = P(q,:) + V(q,:);

        % 边界处理（逐元素裁剪到 [lbv,ubv]）
        P(q,:) = max(P(q,:), lbv);
        P(q,:) = min(P(q,:), ubv);

        % 评估
        f(q) = prob(P(q,:));

        % 更新个体/全局最优
        if f(q) < f_pbest(q)
            f_pbest(q) = f(q);
            pbest(q,:) = P(q,:);
            if f_pbest(q) < f_gbest
                f_gbest = f_pbest(q);
                gbest   = pbest(q,:);
            end
        end
    end
    bestfitIter(t) = f_gbest;   % 逐代全局最优（已单调不增）
end

% Results
best_Particle          = f_gbest;
bestParticle_Position  = gbest;
end
