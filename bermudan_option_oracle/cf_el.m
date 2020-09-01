function el = cf_el(objMean,objSd)

% Calculates expected loss for contour finding
%
% Discription:
%   CF_EL Computes expected loss using the optimal stopping loss function.
%   EL = CF_EL(OBJMEAN, OBJSD) returns the expected loss for contour
%   finding calculated with the optimal stopping loss function.
%
%       OBJMEAN, OBJSD - a list with a mean and an sd field

el = objSd.*normpdf( -abs(objMean)./objSd ) - abs(objMean).*normcdf( -abs(objMean)./objSd);
el(el<0) = 0;

end