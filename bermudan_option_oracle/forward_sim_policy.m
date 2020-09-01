function fsim = forward_sim_policy( x, M, fit, model, offset, compact, use_qv )

% FORWARD_SIM_POLICY Calculates the payoff in Amerian put/call option

% Description:
%   FORWARD_SIM_POLICY calculates the payoff in American put/call option 
%   using forward stopping criteria given the forward boundaries.
%
%   FSIM = FORWARD_SIM_POLICY( X, M, FIT, MODEL, OFFSET, COMPACT, USE_QV )
%   returns the resulting payoff NPV from t=0 (payoff), list of resulting 
%   payoffs (on paths still not stopped) NPV from t=i (fvalue[i]), the
%   times when stopped and the forward x-values of paths at t=i (those not 
%   stopped yet) (sims[[i]]).

%       X - The inputs at t = 0
%       M - The current time stamp
%       FIT - The forward fitted models
%       MODEL - The model to fit the boundary


nsim = 0;

if( nargin == 4 )
    offset = 1;
    compact = true;
    use_qv = false;
end

if ( nargin == 5 )
    compact = true;
    use_qv = false;
end

if ( nargin == 6 )
    use_qv = false;
end

if (ismatrix(x) || isnumeric(x) ) 
    curX = model.sim_func( x, model, model.dt);
    nsim = size(x,1);
end

if ( size(x, 3) > 1 )
    curX = x(:,:,1);
end

payoff = zeros(size(curX,1), 1);
tau = zeros(size(curX,1), 1);
sims = zeros(size(curX,1), 1, model.look_ahead+1);
save_ndx = zeros(size(curX,1), 1, model.look_ahead+1);
fvalue = zeros(size(curX,1), 1, model.look_ahead+1);

contNdx = (1:size(curX,1))';
i = 1;
payoff(contNdx)  = exp(-(i)*model.dt*model.r)*model.option_payoff( curX(contNdx,:), model.K );

% main loop forward
while ( ( i < (M + use_qv) ) && size(contNdx,1) > 0 && ~isempty(contNdx)) 

    in_the_money =  find( payoff(contNdx) > 0); % 1:length(contNdx) %
    
    if ( size(in_the_money,1) >0 )
        
        rule = gp_pred(fit(i+1-offset).gp, fit(i+1-offset).x, fit(i+1-offset).y, curX(contNdx(in_the_money),:));

        if (use_qv == true && i== M) 
          payoff(contNdx) = payoffCont(contNdx) + rule;  % continuation value of paths that didn't stop yet)
          break
        end

        % stop if the expected gain is negative
        if (size(find(rule <0),1) > 0)
           contNdx(in_the_money(rule < 0)) = [];
        end

    end
    
    tau(contNdx) = (i)*model.dt;

    if (compact == false) 
      sims(:,:,min(i,model.look_ahead+1)) = curX(contNdx,:);
      save_ndx(:,:,min(i,model.look_ahead+1)) = contNdx;
    end
    
    % update the x values by taking a step of length dt
    
    i = i+1;

    if (ismatrix(x) || isnumeric(x) ) 
        curX(contNdx,:) = model.sim_func( curX(contNdx,:),model,model.dt);
        nsim = nsim + size(contNdx,1);
    end

    if (size(x,3) > 1)  % stored list of paths
        curX(contNdx,:) = x(contNdx,:,i);
    end
    
    % payoff for next timestep -- used for terminal payoff at i=M
    payoff(contNdx)  = exp(-(i)*model.dt*model.r)*model.option_payoff( curX(contNdx,:), model.K);

end

for i = 2:(model.look_ahead)   % payoff for a trajectory starting at x^n_{t+i} which was still alive then
    fvalue(:,:,i) = payoff( save_ndx(:,:,i) )*exp((i-1)*model.dt*model.r);
end

fsim.payoff = payoff;
fsim.fvalue = fvalue;
fsim.sims = sims;
fsim.tau = tau;
fsim.nsim = nsim;
end
