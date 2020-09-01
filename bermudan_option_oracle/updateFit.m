function fit = updateFit(fit, add_grid, add_mean)

% Updates the Gaussian Process fit
% Description:
%   UPDATEFIT updates the Gaussian Process fit by added new sample and
%   observation
%   FIT = UPDATEFIT(FIT, ADD_GRID, ADD_MEAN, ADD_VAR, R_CURR) returns the
%   updated GP fit
%       FIT - Original fit
%       ADD_GRID - new sample
%       ADD_MEAN - observation mean at new sample

fit.x = [fit.x; add_grid];
fit.y = [fit.y; add_mean];
end
