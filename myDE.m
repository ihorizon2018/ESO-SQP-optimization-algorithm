%% DIFFERENTIAL EVOLUTION (DE)
function [bestPop_Position,best,bestfitIter] = myDE(N,T,lb,ub,dim,fobj)
%% Problem Parameters
prob = fobj;                  % fitness function
Var = dim;                     % number of variables
lb = lb;                        % variable Lower Bound
ub = ub;                        % variable Upper Bound
cPop = N;                    % class population size
Maxiter = T;                % max no of iteration

%% PSO Parameters
pc=0.8;                         % cross over probability
sf=0.45;                        % scaling factor


%% Starting the Differential Evolution
f = NaN(cPop,1);                    % target vector for fitness fun value of the population
fu = NaN(cPop,1);                   % trial vector for fitness fun value of the new population
D = Var;                            % length of decision variables

U = NaN(cPop,D);                    % matrix to store the inital solutions

bestfitIter = NaN(Maxiter,1);     % vetcor to store best fitness fun value per iteration
if isscalar(lb), lbv = lb*ones(1,dim); else, lbv = lb(:)'; end
if isscalar(ub), ubv = ub*ones(1,dim); else, ubv = ub(:)'; end
X = lbv + rand(N,dim).*(ubv - lbv);
P=X;
if size(P,2)~=dim
    P=P';
end

for q = 1:cPop
    f(q) = prob(P(q,:));        % evaluating the fitness fun value of initial population
end

% bestfitIter(1) = min(f);            % vector to store best fitness value of the 0th iteration
%% Iteration Loop (Main Loop)
for t = 1:Maxiter
    for q = 1:cPop
        
        % Mutation
        Candidates = [1:q-1 q+1:cPop];          % ensuring that the current member is not a parter
        idx = Candidates(randperm(cPop-1,3));    % selection of three random parters
        
        S1 = P(idx(1),:);           % assign randomly selected solution 1
        S2 = P(idx(2),:);           % assign randomly selected solution 1
        S3 = P(idx(3),:);           % assign randomly selected solution 1
        
        v = S1 + sf*(S2 - S3);      % generating the donor vector
        
        %% Crossover
        del = randi(D,1);           % genertaing the random variable delta
        for j = 1:D
            
            if (rand <= pc) || del == j
                U(q,j) = v(j);      % check for donor vector or target vector
            else
                U(q,j) = P(q,j);    % check for donor vector or target vector
            end
        end
    end
    
    %% Bounding and Greedy Selection
    for k = 1:cPop
        U(k,:) = max(lb,U(k,:));    % bounding the violating to their lower bound
        U(k,:) = min(ub,U(k,:));    % bounding the violating to their upper bound
        
        fu(k) = prob(U(k,:));       % evaluating the fitness of thr trial solution
        
        if fu(k) < f(k)
            P(k,:) = U(k,:);
            f(k) = fu(k);
        end
    end
    bestfitIter(t) = min(f);



end

% Results
[best,ind] = min(f);
bestPop_Position = P(ind,:);
