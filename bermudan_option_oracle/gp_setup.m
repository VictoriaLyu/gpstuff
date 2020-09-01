function [fit, l, sigma2, sigman, nu] = gp_setup(model, x, y, opt)

% Sets up Gaussian Process model
% Description:
%   GP_SETUP constructs Gaussian Process simulator for American 
%   Put/Call Option function.
%   [FIT,...] = GP_SETUP(MODEL, X, Y, OPT) sets up the Gaussian
%   Process fit for American put/call options.
% 
%       MODEL - the likelihood function used in GP: GP or t-GP
%       X - current designs
%       Y - current observations

d = size(x,2);

x0 = x;

if (nargin == 3)
    opt=optimset('TolFun',1e-3,'TolX',1e-3);
end

gpcf = gpcf_sexp('lengthScale',repmat(10,1,d), 'magnSigma2', 1);

switch model
    case 'gauss'
        gprocess = gp_set('lik', lik_gaussian(), 'cf', gpcf, 'jitterSigma2', 1e-9);
        gprocess = gp_optim(gprocess, x, y, 'opt', opt);
    case 't'
        gprocess = gp_set('lik', lik_t(), 'cf', gpcf, 'jitterSigma2', 1e-9);
        gprocess = gp_optim(gprocess, x, y, 'opt', opt);
    case 'probit'
        gprocess = gp_set('lik', lik_probit(), 'cf', gpcf, 'jitterSigma2', 1e-9);
        gprocess = gp_optim(gprocess, x, y, 'opt', opt);
    case 'mprobit'
        gprocess = gp_set('lik', lik_probit(), 'cf', gpcf, 'jitterSigma2', 1e-9);
        gprocess = gp_optim(gprocess, x, y, 'opt', opt);
        gprocess = gp_monotonic(gprocess, x, y, 'optimize', 'on','opt', opt);
    case 'mgauss'
        gprocess = gp_set('lik', lik_gaussian(), 'cf', gpcf, 'jitterSigma2', 1e-9);
        gprocess = gp_optim(gprocess, x, y, 'opt', opt);
        gprocess = gp_monotonic(gprocess, x, y, 'optimize', 'on','opt', opt);
end

w = gp_pak(gprocess);
sigma2 = exp(w(1))^2;
l = exp(w(2:(2+d-1)));
if (~strcmp(model, 't'))
    sigman = exp(w(end))^2;
    nu = 0;
else
    sigman = exp(w(end - 1))^2;
    nu = exp(w(end));
end

fit.gp = gprocess;
fit.x = x0;
fit.y = y;
fit.t = 0;
end