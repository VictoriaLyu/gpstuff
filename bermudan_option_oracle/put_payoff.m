function y = put_payoff(x,K)

% Arithmetic basket Put for d-dim x
%
% Description:
%   PUT_PAYOFF calculates the arithmetic basket put option for d-dimention
%   x
%   Y = PUT_PAYOFF(X, K) returns the payoff at asset price x and strike K
% 
%       K - the strike price
%       X - a vector of asset prices

if (nargin == 1)
    K = 40;
end

y = K-mean(x,2);
y(y<0) = 0;

end

