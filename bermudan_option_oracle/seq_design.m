function [x_optim, metric] = seq_design(fit, x, y, model, lhs_rect, i, xtest, xt_dens)

% Optimizes the acquisition function in sequential design
% Description:
%   SEQ_DESIGN optimizes the acquisition function in sequential design for
%   American put/call option
%   [X_OPTIM, METRIC] = SEQ_DESIGN(FIT, X, Y, MODEL, LHS_RECT) returns
%   the optimized sampling location and batch size in
%   sequential design
%       FIT - the GP fit used to calculate the acquisition function
%       X - inputs
%       Y - observations
%       MODEL - the list of model information (including GP metamodels and
%       designs)
%       LHS_RECT - the edge of design space

gprocess = fit.gp; % fitted GP
design = model.design; % sequantial design algorithm
d = size(x,2); % dimension

opts = gaoptimset('Display','off','TolFun', 1e-6, 'Generations', 500, 'Vectorized','on','UseParallel',true);
[x_optim, metric]= ga(@(xt) acquiFun(gprocess, x, y, model.km_batch, xt, model, design, i, xtest, xt_dens), d, [],[], [], [], lhs_rect(:,1), lhs_rect(:,2), [], opts);
metric = -metric;
end

function EI = acquiFun(gprocess, x, y, r0, xt, model, design, i, xtest, xt_dens)

% Acquisition function for MCU, tMSE, cSUR and ICU

% Description:
%   ACQUIFUN calculates the weighted acquisition function in sequential
%   design
%   EI = ACQUIFUN(GPROCESS, X, Y, R, XT, MODEL, DESIGN, i) returns the
%   weighted expected improvement at ith time stamp in put/call option
%       GPROCESS - fitted Gaussian Process to calculate the EI
%       X - inputs
%       Y - observations
%       R - batch size
%       XT - candidate samples
%       MODEL - GP model(GP/t-GP/Cl-GP/MCl-GP/MGP)
%       XTEST and XT_DENS - testing samples and their weights for ICU

d = size(x, 2);
switch design
    case 'cSUR'
        [Ef, Varf] = gp_pred(gprocess, x, y, xt);
        Ef = Ef(:,1);
        switch model.method
            case {'gauss', 'mgauss', 't'}
                sigmanoise = gprocess.lik.sigma2/r0;
            case {'probit', 'mprobit'}
                sigmanoise = inverseterm_probit(Ef, Varf);
        end
        % approximates the new variance %
        Varfnew = sigmanoise .* Varf./(sigmanoise + Varf);
        % calculates the reduced empirical error %
        EI = -normcdf( -abs(Ef)./sqrt(abs(Varf)) ) + normcdf( -abs(Ef)./sqrt(abs(Varfnew)));
        
    case 'MCU'
        [Ef, Varf] = gp_pred(gprocess, x, y, xt);
        Ef = Ef(:,1);
        % calculate the local empirical error %
        beta = (quantile(Ef, 0.75) - quantile(Ef, 0.25))/mean(sqrt(Varf));
        EI = abs(Ef) - beta.*sqrt(Varf);
        
    case 'tMSE'
        [Ef, Varf] = gp_pred(gprocess, x, y, xt);
        Ef = Ef(:,1);
        sigmaepsilon = 0;
        w = exp(-0.5.*(Ef./sqrt(Varf + sigmaepsilon^2)).^2)./sqrt(2*pi*(Varf + sigmaepsilon^2));
        % eq. 3.5 %
        EI = -Varf.*w;
    
    otherwise
        n = size(xtest,1);
        [Efall, Covfall] = gp_jpred(gprocess, x, y, [xtest;xt]);
        Varfall = diag(Covfall);
        m = size(xt,1);
        Varftest = repmat(Varfall(1:n),1,m);
        Varf = Varfall((n+1):end);
        Eftest = repmat(Efall(1:n,1),1,m);
        Ef = Efall((n+1):end,1);
        switch model.method
            case {'gauss', 'mgauss', 't'}
                sigmanoise = gprocess.lik.sigma2/r0;
            case {'probit', 'mprobit'}
                sigmanoise = inverseterm_probit(Ef, Varf);
        end
        % approximate the new variance in eq. 4.5 %
        Varfupdate = Varftest - Covfall(1:n,(n+1):end).^2 ./repmat((sigmanoise + Varf)',n,1);
        ilocal = normcdf(-abs(Eftest)./sqrt(abs(Varfupdate)));
        
        ilocal = ilocal.*repmat(xt_dens,1,m);
        % calculate the expected empirical error in eq. 3.4 %
        EI = sum(ilocal)';
end

% density of x %
if (~strcmp(design, 'ICU'))
    x_dens = lognpdf(xt(:,1), log(model.x0(1))+(model.r - model.div - model.sigma(1)^2/2)*i*model.dt, model.sigma(1)*sqrt(i*model.dt));

    if ( d >= 2 )
        for j = 2:d
            x_dens = x_dens.*lognpdf(xt(:,j), log(model.x0(j))+(model.r - model.div - model.sigma(j)^2/2)*i*model.dt, model.sigma(j)*sqrt(i*model.dt));
        end
    end
else
    x_dens = ones(size(xt,1), 1);
end
        
% excludes deep in-the-money designs %
if (d == 2  && model.K == 40)
    if (~strcmp(design, 'ICU'))
        x_dens(xt(:,1) + xt(:,2) > 80) = 0;
    else
        x_dens(xt(:,1) + xt(:,2) > 80) = Inf;
    end
end
if (d == 2  && model.K == 100)
    if (~strcmp(design, 'ICU'))
        x_dens(xt(:,1) < 100 & xt(:,2) < 100) = 0;
    else
        x_dens(xt(:,1) < 100 & xt(:,2) < 100) = Inf;
    end
end
if (d == 3  && model.K == 100)
    if (~strcmp(design, 'ICU'))
        x_dens(xt(:,1) < 100 & xt(:,2) < 100 & xt(:,3) < 100) = 0;
    else
        x_dens(xt(:,1) < 100 & xt(:,2) < 100 & xt(:,3) < 100) = Inf;
    end
end
        
EI = EI.*x_dens;

end

