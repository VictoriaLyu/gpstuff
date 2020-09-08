function [x_seq, y_seq, ee, er, bias, metric, t, Ef, Varf, xtt, l, sigma2, sigman, t_optim, t_gen] = intgp(k, m0, d, budget, fun)

% Initializes the variable space for updatapar.

%%% construct testing points %%%

if ( d == 1 )
    xtt = linspace(0,1,m0)';
end

if ( d == 2 )
    if (~isequal(fun, @braninsc2))
        x1 = repmat(linspace(0,1,m0)',1,m0);
        x2 = repmat(linspace(0,1,m0)',1,m0)';
        xtt = [x1(:) x2(:)];
    else
        % Approximate empirical error rate %
        rng default  % For reproducibility
        p = sobolset(2,'Skip',1e3,'Leap',1e2);
        p = scramble(p,'MatousekAffineOwen');
        xtt = net(p,100000);
        f = fun(xtt);
        xt1 = xtt(abs(f) < 0.5,:);
        xt1 = xt1(1:400,:);
        xt2 = xtt(abs(f) >= 0.5,:);
        xt2 = xt2(1:100,:);
        xtt = [xt1;xt2];
    end
end

if ( d == 6 )
    % Approximate empirical error rate %
    rng default  % For reproducibility
    p = sobolset(6,'Skip',1e3,'Leap',1e2);
    p = scramble(p,'MatousekAffineOwen');
    xtt = net(p,10000000);
    f = fun(xtt);
    xt1 = xtt(abs(f) < 0.7,:);
    xt1 = xt1(1:800,:);
    xt2 = xtt(abs(f) >= 0.7,:);
    xt2 = xt2(1:200,:);
    xtt = [xt1;xt2];
end

%%% initialize variables %%%

m = size(xtt,1);
num_diff_samples = budget; % the variable size is set up to 5000
x_seq = zeros(num_diff_samples,d,k);
y_seq = zeros(num_diff_samples,k);

er = zeros(num_diff_samples,k);
ee = zeros(num_diff_samples,k);
bias = zeros(num_diff_samples,k);
metric = zeros(num_diff_samples,k);
t = zeros(num_diff_samples,k);
t_optim = zeros(num_diff_samples,k);
t_gen = zeros(num_diff_samples,k);

Ef = zeros(m,num_diff_samples,k);
Varf = zeros(m,num_diff_samples,k);

l = zeros(num_diff_samples,d,k);
sigma2 = zeros(num_diff_samples,k);
sigman = zeros(num_diff_samples,k);
end