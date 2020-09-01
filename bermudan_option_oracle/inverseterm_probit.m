function sigmanoise = inverseterm_probit(Ef, Varf)
vpos = normpdf(Ef).^2./normcdf(Ef).^2 + Ef.*normpdf(Ef)./normcdf(Ef);
vneg = normpdf(Ef).^2./normcdf(Ef).^2 - Ef.*normpdf(Ef)./normcdf(Ef);
ppos = normcdf(Ef./sqrt(1+Varf));
pneg = 1 - ppos;
v = vpos.*ppos + vneg.*pneg;
sigmanoise = 1./v;