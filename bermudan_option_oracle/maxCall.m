function y = maxCall(x, K)

% Max Call payoff
%
% Description:
%   MAXCALL calculates the max call payoff for basket of call option with
%   asset prices x at strike K.
%   Y = MAXCALL(X, K) returns the max call payoff.
%       X - asset prices
%       K - strike price

if (nargin == 1)
    K = 100;
end

y = max(x,[],2)-K;

y(y<0) = 0;

end