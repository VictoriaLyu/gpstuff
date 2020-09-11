function [fun, X0, Y0, Kmatinv] = initializeGPTP(d, N0, original_func, noisestructure, noisevar)

% Initializes true function for a Gaussian/T process
%
% Description:
%   FUN = initializeGPTP(D, N0, ORIGINAL_FUN, NOISESTRUCTURE, NOISEVAR) 
%   initializes the Gaussian/T Process as a true function for simulations.
%       D - number of dimentions
%       N0 - number of sampled observations
%       ORIGINAL_FUNC - the type of process (tp/gp)
%       NOISESTRUCTURE - the distribution of noise (Normal/t)
%       NOISEVAR - the variance of noise

% initialize gp function
rng(0);
X0 = lhsdesign(N0, d);
    
gpcf_init = gpcf_sexp('lengthScale',repmat(0.2,1,d), 'magnSigma2', 2);

gprocess_init = gp_set('lik', lik_gaussian('sigma2', 0), 'cf', gpcf_init, 'jitterSigma2', 1e-9);    

% covariance matrix
Kmat = gp_cov(gprocess_init, X0, X0);
Kmatinv = inv(Kmat);

% sample observations
if (isequal(original_func, @gp))
    f = mvnrnd(zeros(N0, 1), Kmat)';
else
    f = mvtrnd(Kmat, 3)';
end

% sample noisy observations
rng(0);
Y0 = genFun(X0, original_func, noisestructure, noisevar, f);

fun = @(x)process_mean(x, X0, Y0, Kmatinv, gprocess_init);
end