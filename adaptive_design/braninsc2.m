function [y] = braninsc2(xx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BRANIN-HOO FUNCTION, RESCALED
%
% The rescaled Branin-Hoo function is monotone along both dimensions and
% has mean 0 and variance 1 over sample space [0, 1]*[0, 1].

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
% xx = [x1, x2]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1 = xx(:,1);
x2 = xx(:,2);

x1bar = 15 .* x1;
x2bar = 15 .* x2 - 5;

term1 = x1bar - 5.1 .* x2bar .^ 2 / (4 * pi^2) + 5 .* x2bar./pi - 20;
term2 = (10 - 10/(8*pi)) .* cos(x2bar);

y = (term1.^2 + term2 - 181.47) / 178;

end
