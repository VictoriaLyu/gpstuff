function [gprocess, l, sigma2, sigman] = gp_setup(model, x, y, opt)

% Sets up/Initializes Gaussian Process simulator
%
% Description:
%   [GPROCESS,...] = GP_SETUP(MODEL, X, Y, OPT) sets up the
%   Gaussian Process simulator GPROCESS
%       MODEL - the likelihood for Gaussian Process (Gaussian/t)
%       X - inputs
%       Y - observations
%       OPT - parameters used in optimization of hyperparameters

d = size(x,2);

if (nargin == 4)
    opt=optimset('TolFun',1e-3,'TolX',1e-3);
end

gpcf = gpcf_sexp('lengthScale',repmat(0.5,1,d), 'magnSigma2', 1);
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
        if (d == 1)
            gprocess = gp_monotonic(gprocess, x, y, 'optimize', 'on','opt', opt);
        else
            gprocess = gp_monotonic(gprocess, x, y, 'optimize', 'on','opt', opt, 'nvd', -1);
        end
    case 'mgauss'
        gprocess = gp_set('lik', lik_gaussian(), 'cf', gpcf, 'jitterSigma2', 1e-9);
        gprocess = gp_optim(gprocess, x, y, 'opt', opt);
        if (d == 1)
            gprocess = gp_monotonic(gprocess, x, y, 'optimize', 'on','opt', opt);
        else
            gprocess = gp_monotonic(gprocess, x, y, 'optimize', 'on','opt', opt, 'nvd', -1);
        end
end

w = gp_pak(gprocess);
sigma2 = exp(w(1))^2;
l = exp(w(2:(2+d-1)));

if (ismember(model, {'gauss', 't', 'mgauss'}))
    sigman = exp(w(end));
else
    sigman = 0;
end
end