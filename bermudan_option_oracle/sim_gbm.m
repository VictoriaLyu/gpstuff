function newX = sim_gbm( x0, model, dt )

% Simulates paths of GBM
%
% Description:
%   SIM_GBM simulates the design weight from p(X_t|X_{t-1}) using 
%   log-normal transition density specified by the model
%   NEWX = SIM_GBM(X0, MODEL, DT) returns the simulated paths from
%   Geometric Brownian Motion.
%       x0 - the starting values
%       dt - the step size
%       model - contains all the other parameters in GBM

if (nargin == 2)
    dt = model.dt;
end

len = size(x0,1);

newX = x0;

for j = 1:model.dim
   newX(:,j) = x0(:,j).*exp( normrnd(0,1,[len,1])*model.sigma(j)*sqrt(dt) + (model.r- model.div- model.sigma(j)^2/2)*dt );
end

end