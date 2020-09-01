function [y, t]= genFun(X, fun, noisestructure, noisevar, f, r)

% Generates noisy observation for the true function
%
% Description:
%   [Y,...] = GENFUN(X, FUN, NOISESTRUCTURE, NOISEVAR, R) generates noisy
%   observations Y
%       X - sampling location
%       FUN - true function
%       NOISESTRUCTURE - the structure of noise (normal/t with constant
%       degree of freedom/t with heteroskedastic degree of freedom)
%       NOISEVAR - the variance of the noise (small/large/heteroskedastic)
%       R - batch size

gen_t = tic;

if (nargin < 5)
    f = fun(X);
end

if (nargin < 6)
    r = 1;
end

n = size(f,1);

switch noisevar
    case 'small'
        sd = 0.1;
    case 'middle'
        sd = 0.2;
    case 'large'
        sd = 1;
    case 'hetero'
        sd = 0.4.*(4.*X(:,1)+1);
    case 'mixed'
        sd = [0.5 1];
end

if (isscalar(r) && r >= 1)
    % handles constant r %
    switch noisestructure
        case 'normal'
            ytmp = random('normal', repmat(f,1,r), sd, n, r);
        case 't_constdf'
            ytmp = random('tlocationscale', repmat(f,1,r), sd, 3, n, r);
        case 't_heterodf'
            ytmp = random('tlocationscale', repmat(f,1,r), repmat(sd, 1, r), repmat(6-4*X(:,1), 1, r), n, r);
        case 'mixed'
            f_grid = repelem(f, r);
            mu_grid = [f_grid f_grid]';
            s1_grid = sd(1)^2 .* diag(ones(r * n, 1));
            s2_grid = sd(2)^2 .* diag(ones(r * n, 1));
            ytmp = random(gmdistribution(mu_grid, cat(3, s1_grid, s2_grid), [0.5, 0.5]));
            ytmp = reshape(ytmp, r, n)';
        case 'null'
            ytmp = zeros(n,r);
            for j = 1:r
                ytmp(:,j) = fun(X);
            end
    end
    y = mean(ytmp,2);
else
    r_total = sum(r);
    id = (1:n)';
    id = repelem(id, r);
    f_grid = repelem(f, r);
    switch noisestructure
        case 'normal'
            ytmp = random('normal', f_grid, sd);
        case 't_constdf'
            ytmp = random('tlocationscale', f_grid, sd, 3);
        case 't_heterodf'
            sd_grid = repelem(sd, r);
            nu_grid = repelem(6-4*X(:, 1), r);
            ytmp = random('tlocationscale', f_grid, sd_grid, nu_grid);
        case 'mixed'
            mu_grid = [f_grid f_grid]';
            s1_grid = sd(1)^2 .* diag(ones(r_total, 1));
            s2_grid = sd(2)^2 .* diag(ones(r_total, 1));
            ytmp = random(gmdistribution(mu_grid, cat(3, s1_grid, s2_grid), [0.5, 0.5]));
        case 'null'
            ytmp = f_grid;
    end
    sums  = accumarray(id, ytmp);
    y = sums ./ r;
end
t = toc(gen_t);
end
