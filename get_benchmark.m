function [lb, ub, dim, fobj] = get_benchmark(name)
    key = lower(strtrim(name));
    dim = 30;

    switch key
        case 'sphere'
            lb = -100; ub = 100; fobj = @(x) sum(x.^2);

        case 'sumproduct'
            lb = -10;  ub = 10;  fobj = @(x) sum(abs(x)) + prod(abs(x));

        case 'maxvalue'
            lb = -10;  ub = 10;  fobj = @(x) max(abs(x));

        case 'rosenbrock'
            lb = -30;  ub = 30;  fobj = @(x) sum(100*(x(2:end)-x(1:end-1).^2).^2 + (x(1:end-1)-1).^2);

        case 'schwefel'
            lb = -500; ub = 500; fobj = @(x) 418.9829*numel(x) - sum(x.*sin(sqrt(abs(x))));

        case 'rastrigin'
            lb = -5.12; ub = 5.12; fobj = @(x) 10*numel(x) + sum(x.^2 - 10*cos(2*pi*x));

        case 'penalized1'
            lb = -50;  ub = 50;  fobj = @(x) penalized1(x);

        case 'penalized2'
            lb = -50;  ub = 50;  fobj = @(x) penalized2(x);

        otherwise
            error('Unknown function name: %s', name);
    end
end

function y = penalized1(x)
    x = x(:)';
    n = numel(x);
    term1 = pi/n*(10*sin(pi*(1+(x(1)+1)/4))^2);
    term2 = sum( ((x(1:n-1)+1)/4).^2 .* (1+10*(sin(pi*(1+(x(2:n)+1)/4))).^2) );
    term3 = ((x(n)+1)/4)^2;
    U = Ufun(x,10,100,4);
    y = term1 + term2 + term3 + sum(U);
end

function y = penalized2(x)
    x = x(:)';
    n = numel(x);
    term = (sin(3*pi*x(1)))^2 + sum((x(1:n-1)-1).^2 .* (1+(sin(3*pi*x(2:n))).^2)) + (x(n)-1)^2*(1+(sin(2*pi*x(n)))^2);
    U = Ufun(x,5,100,4);
    y = 0.1*term + sum(U);
end

function y = Ufun(x,a,k,m)
    y = k*((x>a).*(x-a).^m + (x<-a).*(-x-a).^m);
end
