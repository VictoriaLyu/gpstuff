function Imee = metricmee(Ef, Varf)
Imee = normcdf(-abs(Ef)./sqrt(abs(Varf)));
end