function [Xfood, Ffood, gbest_t2, fminconX, fminconF, curve] = myESOSQP(N, T, lb, ub, dim, fobj)
% myESOSQP
% Runs your ESO first (via myESO), then a local SQP refinement with fmincon.
% Signature unchanged from your version.

% --- bounds as row vectors ---
if isscalar(lb), lbv = lb * ones(1,dim); else, lbv = lb(:)'; end
if isscalar(ub), ubv = ub * ones(1,dim); else, ubv = ub(:)'; end

% --- Phase 1: ESO global search (3/4 of the budget) ---
T_eso = max(1, floor(3*T/4));
[Xfood, Ffood, gbest_t] = myESO(N, T_eso, lb, ub, dim, fobj);  

% Ensure ESO output is within bounds before passing to SQP
Xfood = min(max(Xfood(:)'.*1, lbv), ubv);

% --- Phase 2: SQP local refinement (fmincon) ---
% Objective wrapper to ensure row-vector input to fobj
obj = @(x) fobj(x(:)');

% Budgeting
Mconnt = (T - T_eso) * max(N,1);              % rough FE budget for local phase
MaxIterLocal = max(1, floor(T/2));            % your original intent

% Capture per-iteration values via nested OutputFcn (no globals)
curve_rec = [];
function stop = outfun(~, optimValues, state)
    stop = false;
    if strcmp(state,'iter')
        curve_rec(end+1) = optimValues.fval;
    end
end

opts = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'off', ...
    'MaxIterations', MaxIterLocal, ...
    'MaxFunctionEvaluations', max(1000, 10*Mconnt), ...
    'StepTolerance', 1e-12, ...
    'OptimalityTolerance', 1e-12, ...
    'SpecifyObjectiveGradient', false, ...
    'OutputFcn', @outfun);

A = []; b = []; Aeq = []; beq = []; nonlcon = [];

[fminconX, fminconF] = fmincon(obj, Xfood, A, b, Aeq, beq, lbv, ubv, nonlcon, opts);

% SQP curve: best-so-far, row vector
if isempty(curve_rec)
    curve = fminconF;
else
    curve = cummin(curve_rec);
end
curve = curve(:)';

% --- Final outputs ---
% Concatenate ESO curve and SQP curve
gbest_t = gbest_t(:)';            % ensure row
gbest_t2 = [gbest_t, curve];      % combined convergence history

% Choose the better final solution
if fminconF < Ffood
    Xfood = fminconX;
    Ffood = fminconF;
end

% Ensure outputs are row-shaped as expected elsewhere
Xfood = Xfood(:)';

end
