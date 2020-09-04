function [x_seq, y_seq, ee, er, bias, metric, t, Ef, Varf, l, sigma2, sigman, t_optim, t_gen] = updategppar(I, k, m0, d, budget, fun, noisestructure, noisevar, model, design)

% Generates and records the synthetic experiments results in parallel
% computing. 
%
% Description:
%   [x_seq, y_seq,...] = UPDATEGPPAR(I, k, m0, d, budget, fun, 
%   noisestructure, noisevar, model, design, r) returns the
%   entire designs x_seq with observations y_seq for a
%   given true function and a given noise structure.

%   PARAMS:
%       I - the initial design
%       k - number of runs of experiments
%       m0 - size of test set 
%       d - dimension
%       budget - total budget (N)
%       fun - the true function
%       noisestructure - the distribution of noise: 'normal', 't_constdf',
%       't_heterodf', 'mixed'
%       noisevar - variance of noise: 'small' (0.1), 'large'(1), 'mixed', 
%       'hetero'
%       model - Gaussian Process ('gauss') or t-GP ('t') or Cl-GP
%       ('probit') or MGP ('mgauss') or MCl-GP ('mprobit')
%       design - 'cUCB', 'gSUR', 'tMSE' or 'SUR'

%   Outputs:
%       x_seq - designs 
%       y_seq - noisy observations
%       ee - empirical error 
%       er - error rate 
%       bias - bias 
%       metric - optimized acquisition function value 
%       t - total time 
%       Ef - posterior mean function 
%       Varf - posterior function variance 
%       l - lengthscale 
%       sigma2 - sigma2 in covariance function 
%       sigman - sigma2 in likelihood 
%       t_optim - time to optimize acquisition function 
%       t_gen - time to generate observation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize variables %
[x_seq, y_seq, ee, er, bias, metric, t, Ef, Varf, xtt, l, sigma2, sigman, t_optim, t_gen] = intgp(k, m0, d, budget, fun);


% critical probability pcr in eq. 5.2 %
pcr = 0.4;
lambda = 0.8;
ft = fun(xtt);

% do experiments in parellel %
parfor j = 1:k
        
    rng(j+2);

    Xint = lhsdesign(I,d);

    % update gp %
    [x_seq_tmp, y_seq_tmp, metric_tmp, t_tmp, ee_tmp, er_tmp, bias_tmp, Ef_tmp, Varf_tmp, l_tmp, sigma2_tmp, sigman_tmp, t_optim_tmp, t_gen_tmp] = updategp(fun, noisestructure, noisevar, Xint, xtt, model, design, budget, ft, pcr, lambda);

    
    % record performance %
    l(:,:,j) = l_tmp;
    sigma2(:,j) = sigma2_tmp;
    sigman(:,j) = sigman_tmp;
    metric(:, j) = metric_tmp;

    Ef(:,:,j) = Ef_tmp;
    Varf(:,:,j) = Varf_tmp;
    x_seq(:,:,j) = x_seq_tmp;
    y_seq(:,j) = y_seq_tmp;
    
    ee(:,j) = ee_tmp;
    er(:,j) = er_tmp;
    bias(:,j) = bias_tmp;
    t(:,j) = t_tmp;
    t_optim(:,j) = t_optim_tmp;
    t_gen(:,j) = t_gen_tmp;
end
end