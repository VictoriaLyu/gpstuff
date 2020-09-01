function [x_seq, y_seq, metric, time, ee, er, bias, Ef, Varf, l, sigma2, sigman, t_optim, t_gen] = updategp(fun, noisestructure, noisevar, Xint, xtt, model, design, budget, ft, pcr, lambda)

% Generates the synthetic function in each run
%
% Description:
%   [x_seq, y_seq, ...] = UPDATEGP(fun, noisestructure, noisevar, 
%   Xint, xt, model, design, budget, xtt, r, ft, pcr, lambda, batch, t0)
%   generates the designs, observations and batch size in each run of the
%   synthetic experiments for a given true function with a given noise
%   structure.

%   PARAMS:
%       budget - total budget (N)
%       fun - the true function
%       noisestructure - the distribution of noise: 'normal', 't_constdf',
%       't_heterodf', 'mixed'
%       noisevar - variance of noise: 'small' (0.1), 'large'(1), 'mixed', 
%       or 'hetero'
%       model - GP or t-GP
%       design - 'MCU', 'cSUR', 'ABSUR'
%       batch - 'FB', 'MLB', 'RB', 'ABSUR'
%       t0 - hyperparameter in ABSUR
%       lambda and pcr - hyperparameters in estimating the model 
%       performance for test set xtt

%   Outputs:
%       x_seq - designs 
%       y_seq - noisy observations
%       r_seq - batch size
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
%       steps - design size
%       t_optim - time to optimize acquisition function 
%       t_gen - time to generate observation
%       overhead - c_over in ABSUR/ the estimated overhead to optimize the 
%       acquisition function 
%       gamma - hyperparameter in MLB and RB 

if (nargin < 10)
    pcr = 0;
    lambda = 0;
end

%%%%%%% Initialization %%%%%%%
% n0 %
I = size(Xint,1);
% m %
m = size(xtt,1); 
num_diff_samples = budget; % max length of data 
% dimension %
d = size(Xint,2);
% Initialize observations y_{1:n0} %
yint = genFun(Xint, fun, noisestructure, noisevar);

if ( strcmp(model,'probit') || strcmp(model,'mprobit') )
    yint = 2*(yint>0)-1;
end

% initializes design x and observation y
x_seq = zeros(num_diff_samples,d);
y_seq = zeros(num_diff_samples,1);

x_seq(1:I,:) = Xint;
y_seq(1:I) = yint;

% initializes metric, empirical error, error rate, bias %
metric = zeros(num_diff_samples, 1);
ee = zeros(num_diff_samples, 1);
er = zeros(num_diff_samples, 1);
bias = zeros(num_diff_samples, 1);

% initializes hyperparemeters %
l = zeros(num_diff_samples, d);
sigma2 = zeros(num_diff_samples, 1);
sigman = zeros(num_diff_samples, 1);
time = zeros(num_diff_samples, 1);
t_optim = zeros(num_diff_samples, 1);
t_gen = zeros(num_diff_samples, 1);
opt=optimset('TolFun',1e-3,'TolX',1e-3);

% initializes estimated f, and var %
Ef = zeros(m, num_diff_samples);
Varf = zeros(m,num_diff_samples);

budget = budget - I;
step = 1;
while (budget > 0)
    total_t = tic;
    
    % updates the hyperparameters every 80 budgets and record the performance %
    if (step == 1 || mod(step, 10) == 0)
        % estimate hyperparameters %
        [gprocess, l(step+I,:), sigma2(step+I), sigman(step+I)] = gp_setup(model, x_seq(1:(I+step-1),:), y_seq(1:(I+step-1)), opt);
    end
    
    % makes predictions %
    [Ef(:,step+I), Varf(:,step+I)] = gp_pred(gprocess, x_seq(1:(I+step-1),:), y_seq(1:(I+step-1)), xtt);

    % calculates the performance metric %
    [er(I+step), ee(I+step), bias(I+step)] = gp_perf(Ef(:,step+I), Varf(:,step+I), xtt, ft, pcr, lambda);
    
    % optimizes the performance metric and chooses the next design %
    [x_seq(I+step,:), metric(I+step), t_optim(I+step)] = seq_design(gprocess, x_seq(1:(I+step-1),:), y_seq(1:(I+step-1)), design, xtt);
    
    % generates next observation y
    [y_seq(I+step), t_gen(I+step)] = genFun(x_seq(I+step,:), fun, noisestructure, noisevar);
    if ( strcmp(model,'probit') || strcmp(model,'mprobit') )
        y_seq(I+step) = 2*(y_seq(I+step)>0)-1;
    end
    
    time(step+I) = toc(total_t);
    step = step + 1;
    budget = budget - 1;
end
end

        
