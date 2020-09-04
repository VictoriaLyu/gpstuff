function [xoptim, metric, t_optim] = seq_design(gprocess, x, y, design, xtest)

% Optimizes the acquisition function in sequential design and selects the
% next input
% Description:
%   SEQ_DESIGN optimizes the acquisition function in sequential design 
%   [XOPTIM, METRIC,...] = SEQ_DESIGN(GPROCESS, X, Y, DESIGN) returns the optimized 
%   sampling location in sequential design
%
%       GPROCESS - the GP fit used to calculate the acquisition function
%       X - inputs
%       Y - observations
%       DESIGN - sequential design algorithm (cUCB/gSUR/tMSE/SUR)

optim_t = tic;

[xoptim, metric] = optimAcquiFun(gprocess, x, y, design, xtest);        
metric = -metric;
t_optim = toc(optim_t);
end

function EI = acquiFun(gprocess, x, y, xt, design, xtest)
% Description:
%   ACQUIFUN calculates the acquisition function in sequential
%   design
%       GPROCESS - fitted Gaussian Process to calculate the EI
%       X - inputs
%       Y - observations
%       XT - candidate samples
%       DESIGN - design algorithm (cUCB/tMSE/gSUR/SUR)
%       XTEST - testing inputs to calculate SUR acquisition function

switch design
    case 'cUCB'
        [Ef, Varf] = gp_pred(gprocess, x, y, xt);
        gamma = (quantile(Ef, 0.75) - quantile(Ef, 0.25))/mean(sqrt(Varf));
        EI = abs(Ef) - gamma.*sqrt(Varf)./3;
    case 'gSUR'
        [Ef, Varf] = gp_pred(gprocess, x, y, xt);
        if(isfield(gprocess.lik, 'sigma2'))
            sigmanoise = gprocess.lik.sigma2;
        else
            sigmanoise = inverseterm_probit(Ef, Varf);
        end
        % approximate the new variance in eq. 4.3 %
        Varfnew = sigmanoise .* Varf./(sigmanoise + Varf);
        % calculate the reduced empirical error in eq. 3.2 %
        EI = -normcdf( -abs(Ef)./sqrt(abs(Varf)) ) + normcdf( -abs(Ef)./sqrt(abs(Varfnew)) );
    case 'tMSE'
        [Ef, Varf] = gp_pred(gprocess, x, y, xt);
        sigmaepsilon = 0;
        w = exp(-0.5.*(Ef./sqrt(Varf + sigmaepsilon.^2)).^2)./sqrt(2*pi*(Varf + sigmaepsilon.^2));
        % eq. 3.5 %
        EI = -Varf.*w;
    case 'SUR'
        pcr = 0.4;
        [Efall, Covfall] = gp_jpred(gprocess, x, y, [xtest;xt]);
        n = size(xtest,1);
        d = size(xtest,2);
        Varfall = diag(Covfall);
        m = size(xt,1);
        Varftest = repmat(Varfall(1:n),1,m);
        Varf = Varfall((n+1):end);
        Eftest = repmat(Efall(1:n,1),1,m);
        Ef = Efall((n+1):end,1);
        if(isfield(gprocess.lik, 'sigma2'))
            sigmanoise = gprocess.lik.sigma2;
        else
            sigmanoise = inverseterm_probit(Ef, Varf);
        end
        % approximate the new variance in eq. 4.5 %
        Varfupdate = Varftest - Covfall(1:n,(n+1):end).^2 ./repmat((sigmanoise + Varf)',n,1);
        ilocal = normcdf(-abs(Eftest)./sqrt(abs(Varfupdate)));

        % calculate the expected empirical error in eq. 3.4 %
        if (d > 2)
            EI = (pcr*sum(ilocal(1:(0.8*n),:))/(0.8*n)+(1-pcr)*sum(ilocal((0.8*n+1):end,:))/(0.2*n))';
        else
            EI = sum(ilocal)';
        end
       
    case 'SUR-C'
        pcr = 0.4;
        [Efall, Covfall] = gp_jpred(gprocess, x, y, [xtest;xt]);
        n = size(xtest,1);
        d = size(xtest,2);
        Varfall = diag(Covfall);
        m = size(xt,1);
        Varftest = repmat(Varfall(1:n),1,m);
        Varf = Varfall((n+1):end);
        Eftest = repmat(Efall(1:n,1),1,m);
        Ef = Efall((n+1):end,1);
        if(isfield(gprocess.lik, 'sigma2'))
            sigmanoise = gprocess.lik.sigma2;
        else
            sigmanoise = inverseterm_probit(Ef, Varf);
        end
        % approximate the new variance in eq. 4.5 %
        Varfupdate = Varftest - Covfall(1:n,(n+1):end).^2 ./repmat((sigmanoise + Varf)',n,1);
        a_x = Eftest ./ sqrt(abs(Varfupdate));
        c_x = Varftest ./ Varfupdate;
        p_x = normcdf(Eftest./sqrt(abs(Varftest)));
        phi_x = normcdf(Eftest./sqrt(abs(c_x - 1)));
        phi_2 = zeros(n, m);
        mu = zeros(1, 2);
        for i = 1:n
            for j = 1:m
                phi_2(i, j) = mvncdf([a_x(i, j), -a_x(i, j)], mu, [c_x(i, j) + 1e-8 1 - c_x(i, j); 1 - c_x(i, j) c_x(i, j) - 1 + 1e-8]);
            end
        end
        ilocal = 2 * phi_2 - p_x + phi_x; 

        % calculate the expected empirical error in eq. 3.4 %
        if (d > 1)
            EI = (pcr*sum(ilocal(1:(0.8*n),:))/(0.8*n)+(1-pcr)*sum(ilocal((0.8*n+1):end,:))/(0.2*n))';
        else
            EI = sum(ilocal)';
        end
        
    case 'gSUR-C'
        [Ef, Varf] = gp_pred(gprocess, x, y, xt);
        
        if(isfield(gprocess.lik, 'sigma2'))
            sigmanoise = gprocess.lik.sigma2;
        else
            sigmanoise = inverseterm_probit(Ef, Varf);
        end
        % approximate the new variance in eq. 4.3 %
        Varfnew = sigmanoise .* Varf./(sigmanoise + Varf);
        a_x = Ef ./ sqrt(abs(Varfnew));
        c_x = Varf ./ Varfnew;
        p_x = normcdf(Ef./sqrt(abs(Varf)));
        phi_x = normcdf(Ef./sqrt(abs(c_x - 1)));
        phi_2 = zeros(size(xt, 1), 1);
        mu = zeros(1, 2);
        for i = 1:size(xt, 1)
            phi_2(i) = mvncdf([a_x(i), -a_x(i)], mu, [c_x(i) 1 - c_x(i); 1 - c_x(i) c_x(i) - 1]);
        end
        one_step_ahead_prob = 2 * phi_2 - p_x + phi_x; 
        % calculate the reduced empirical error in eq. 3.2 %
        EI = -normcdf( -abs(Ef)./sqrt(abs(Varf)) ) + one_step_ahead_prob;
end
end

function [xoptim, metric] = optimAcquiFun(gprocess, x, y, design, xtest)
d = size(x, 2);
opts = gaoptimset('Display','off','TolFun', 1e-6, 'Generations', 500, 'Vectorized','on','UseParallel',true);
[xoptim, metric]= ga(@(xt) acquiFun(gprocess, x, y, xt, design, xtest), d, [],[], [], [], zeros(d,1), ones(d,1),[], opts);
end
