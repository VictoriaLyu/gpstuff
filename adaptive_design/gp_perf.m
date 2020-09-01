function [er, ee, bias] = gp_perf(Ef, Varf, xt, ft, pcr, lambda)

% Calculates the performance metric for level set estimation.
%
% Description:
%   [ER, EE, BIAS] = GP_PERF(EF, VARF, XT, FT, PCR, LAMBDA) calculates
%   the error rate, integrated empirical error and
%   bias for level set estimation.
%       EF - posterior mean of function
%       VARF - posterior variance of function
%       XT - test samples
%       FT - true function value at test samples
%       pcr and lambda - parameters in 2d and 6d "trick" to approximate the
%       integrated error rate and empirical error

d = size(xt,2);
m = size(xt,1);
nt1 = lambda*m;

lee = metricmee(Ef, Varf);

if ( d == 1 || (pcr == 0 && lambda == 0))
    ee = sum(lee)/m;
    er = sum(Ef.*ft<0)/m;
    bias = (sum((Ef>0)&(ft<0))-sum((Ef<0)&(ft>0)))/m;
else
    er = pcr*sum(Ef(1:nt1).*ft(1:nt1)<0)/nt1 + (1-pcr)*sum(Ef((nt1+1):end).*ft((nt1+1):end)<0)/(m - nt1);
    ee = pcr*sum(lee(1:nt1))/nt1+(1-pcr)*sum(lee((nt1+1):end))/(m - nt1);
    bias = pcr*(sum((Ef(1:nt1)>0)&(ft(1:nt1)<0))-sum((Ef(1:nt1)<0)&(ft(1:nt1)>0)))/nt1 + (1-pcr)*(sum((Ef((nt1+1):end)>0)&(ft((nt1+1):end)<0))-sum((Ef((nt1+1):end)<0)&(ft((nt1+1):end)>0)))/(m - nt1);
end

end