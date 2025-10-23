
clc; clear; close all;

% ==========
N        = 50;
MaxIt    = 200;
func_name= 'sumproduct';
use_logy = true;   
[lb, ub, dim, fobj] = get_benchmark(func_name);
algos  = {'ESO','ESOSQP','SO','DE','PSO','GA','WOA'};
colors = lines(numel(algos));
legend_entries = cell(1, numel(algos));

figure('Position',[200 200 720 440]); hold on; grid on;
xlabel('Iteration'); ylabel('Objective');
if use_logy
    set(gca,'YScale','log');
else
    set(gca,'YScale','linear');
end

for ai = 1:numel(algos)
    aname = algos{ai};
    fprintf('Running %s...\n', aname);
    try
        switch upper(aname)
            case 'ESO'
                [~, ~, curve] = myESO(N, MaxIt, lb, ub, dim, fobj);
            case 'ESOSQP'
                [~, ~, curve] = myESOSQP(N, MaxIt, lb, ub, dim, fobj);
            case 'SO'
                [~, ~, curve] = mySO(N, MaxIt, lb, ub, dim, fobj);
            case 'DE'
                [~, ~, curve] = myDE(N, MaxIt, lb, ub, dim, fobj);
            case 'PSO'
                [~, ~, curve] = myPSO(fobj, lb*ones(1,dim), ub*ones(1,dim), dim, N, MaxIt);
            case 'GA'
                [~, ~, curve] = myGA(fobj, lb*ones(1,dim), ub*ones(1,dim), dim, N, MaxIt);
            case 'WOA'
                [~, ~, curve] = myWOA(N, MaxIt, lb, ub, dim, fobj);
            otherwise
                error('Unknown algorithm: %s', aname);
        end
    catch ME
        warning('Algorithm %s failed: %s', aname, ME.message);
      
        curve = nan(1, MaxIt);
    end

    curve = curve(:)';                
    curve(~isfinite(curve)) = NaN;    
    curve = max(curve, 1e-5);          
    % ==========
    if use_logy
        semilogy(curve, 'LineWidth', 1.5, 'Color', colors(ai,:));
    else
        plot(curve, 'LineWidth', 1.5, 'Color', colors(ai,:));
    end
    legend_entries{ai} = aname;
end

legend(legend_entries, 'Location','northeast');
box on;
