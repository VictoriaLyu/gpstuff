function [y] = hart6(xx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HARTMANN 6-DIMENSIONAL FUNCTION, RESCALED
%
% The rescaled function has mean 0 and variance 1 in sample space [0, 1]^6.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
% xx = [x1, x2, x3, x4, x5, x6]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = [0.2, 0.22, 0.28, 0.3]';
A = [8, 3, 10, 3.50, 1.7, 6;
     0.5, 8, 10, 1, 6, 9;
     3, 3.5, 1.7, 8, 10, 6;
     10, 6, 0.5, 8, 1, 9];
P = 10^(-4) * [1312, 1696, 5569, 124, 8283, 5886;
               2329, 4135, 8307, 3736, 1004, 9991;
               2348, 1451, 3522, 2883, 3047, 6650;
               4047, 8828, 8732, 5743, 1091, 381];

n = size(xx,1);
outer = zeros(n,1);
for ii = 1:4
	inner = zeros(n,1);
	for jj = 1:6
		xj = xx(:,jj);
		Aij = A(ii, jj);
		Pij = P(ii, jj);
		inner = inner + Aij.*(xj-Pij).^2;
	end
	new = alpha(ii) .* exp(-inner);
	outer = outer + new;
end

y = -outer;

y = (y + 0.1)/0.1;
end
