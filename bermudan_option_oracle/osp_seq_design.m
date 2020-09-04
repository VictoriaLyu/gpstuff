function[fit, timeElapsed, nsims, empLoss] = osp_seq_design(model)

% Optimal stopping criteria for American put/call option
% 
% OSP_SEQ_DESIGN calculates the optimal stopping criteria for American
% put/call option
% [FIT, TIMEELAPSED, NSIMS, EMPLOSS] = OSP_SEQ_DESIGN(MODEL) returns a list
% containing:
% fit - a list of fitted response surfaces. 
% timeElapsed - time cost in seconds
% nsims - total number of 1-step sim.func calls
% empLoss - vector of empirical losses

% PARAMS:
%       model - a list including all parameters for American put/call
%       option simulation, and the model (GP/t-GP) and designs
%       (cUCB/tMSE/gSUR/SUR) to estimate the stoppint criteria


M = model.T/model.dt;
offset = 0;
compact = true;

tic;

cur_sim = 0;

d = model.dim;

method = model.method;
design = model.design;
option_payoff = model.option_payoff;
budget = model.budget;

fits = repmat(struct('gp',[],'x',[],'y',[], 't', []), M, 1);   % list of emulator objects at each step

emp_loss = zeros(M,model.adaptive_grid_loop-model.init_size);

update_kernel_iters = (0:10:model.adaptive_grid_loop)';   % when to refit the whole GP
 
% % set-up initial grid to understand the distribution of X

init_grid = model.init_grid;

% % testing size for SUR 
switch d
    case 1
        nt = 100;
    case 2
        nt = 400;
    case 3
        nt = 1000;
end

%%%%%%%%%%%% step back in time
for i = (M-1):-1:1
    t_one_pass = tic;
    
    running_budget = 0;
    all_X = zeros(model.adaptive_grid_loop, d+3);

    % construct the input domain where candidates will be looked for
    if (isempty(model.lhs_rect))
      model.lhs_rect = 0.02;
    
      lhs_rect = zeros(d, 2);
      % create a box using empirical quantiles of the init.grid cloud
      for jj = 1:d
        lhs_rect(jj,:) = quantile( init_grid(:,jj), [model.lhs_rect, 1-model.lhs_rect] );
      end
    else   % already specified
      lhs_rect = model.lhs_rect;
    end
    
    if ( strcmp(design, 'lhs') )
        if ( d == 2 && model.K == 40 )
            X_int = lhsCons( model.adaptive_grid_loop*5, lhs_rect );
            X_int = X_int(X_int(:,1) + X_int(:,2) <= 80, :);
            all_X(:,1:d) = X_int(1:model.adaptive_grid_loop,:);
        else
            if (d == 3  && model.K == 100)
                X_int = lhsCons( model.adaptive_grid_loop*5, lhs_rect );
                X_int = X_int(max([X_int(:,1), X_int(:,2), X_int(:,3)], [], 2) >= 100, :);
                all_X(:,1:d) = X_int(1:model.adaptive_grid_loop,:);
            else
                all_X(:,1:d) = lhsCons( model.adaptive_grid_loop, lhs_rect );
            end
        end
        
        K0 = size(all_X,1);
        
        big_grid = repmat(all_X(:,1:d), model.km_batch, 1);

        fsim = forward_sim_policy( big_grid, M-i, fits(i:M), model, offset, compact);
        cur_sim = cur_sim + fsim.nsim;
    
        % payoff at t
        immPayoff = option_payoff( all_X(:,1:d), model.K );

        % batched mean and variance
        for jj = 1:K0 
          all_X(jj, d+1) = mean( fsim.payoff( jj + (0:K0:K0*(model.km_batch-1))' ) ) - immPayoff(jj);
          all_X(jj, d+2) = var( fsim.payoff( jj + (0:K0:K0*(model.km_batch-1))' ) ) + 0.00001;
        end

        if (strcmp(method, 'probit') || strcmp(method, 'mprobit') )
            all_X(:, d+1) = 2 * (all_X(:, d+1) >= 0) - 1;
        end
        
        all_X(:, end) = model.km_batch;
        % create the gp object
        fits(i) = gp_setup(model.method, all_X(:,1:d), all_X(:,d+1));
        
    else
        
        % initial design
        if ( isempty(model.init_grid) )
            init_grid = lhsCons( model.init_size, lhs_rect);
        else
            init_grid = model.init_grid;
        end

        K0 = size(init_grid,1);

       % initial conditions for all the forward paths: replicated design with km.batch
        big_grid = repmat(init_grid, model.km_batch, 1);

        fsim = forward_sim_policy( big_grid, M-i, fits(i:M), model, offset, compact);
        cur_sim = cur_sim + fsim.nsim;

        % payoff at t
        immPayoff = option_payoff( init_grid, model.K );

        % batched mean and variance
        for jj = 1:K0 
          all_X(jj, d+1) = mean( fsim.payoff( jj + (0:K0:K0*(model.km_batch-1))' ) ) - immPayoff(jj);
          all_X(jj, d+2) = var( fsim.payoff( jj + (0:K0:K0*(model.km_batch-1))' ) ) + 0.00001;
        end
        
        if (strcmp(method, 'probit') || strcmp(method, 'mprobit') )
            all_X(1:K0, d+1) = 2 * (all_X(1:K0, d+1) >= 0) - 1;
        end

        all_X(1:K0, 1:d) = init_grid;  % use  first dim+1 columns for batched GP regression
        all_X(1:K0, end) = model.km_batch;
        
        % create the gp object
        fits(i) = gp_setup(model.method, init_grid, all_X(1:K0,d+1));
        
        % active learning loop
        add_more_sites = true; 
        k = K0;
        running_budget = running_budget + model.km_batch * K0;
        step = 1;
        
        % test points for SUR
        if (strcmp(design, 'SUR'))
            xtest = lhsCons(nt, lhs_rect);
            xt_dens = lognpdf(xtest(:,1), log(model.x0(1))+(model.r - model.div - model.sigma(1)^2/2)*i*model.dt, model.sigma(1)*sqrt(i*model.dt));

            if ( d >= 2 )
                for j = 2:d
                    xt_dens = xt_dens.*lognpdf(xtest(:,j), log(model.x0(j))+(model.r - model.div - model.sigma(j)^2/2)*i*model.dt, model.sigma(j)*sqrt(i*model.dt));
                end
            end
        else
            xtest = 0;
            xt_dens = 0;
        end
        % to do it longer for later points to minimize error build-up use: *(1 + i/M)
        %  for (k in (K0+1):round(model$adaptive.grid.loop))
 
        while (add_more_sites)
    
            % sequential design to optimize the sampling location and
            % design batch
            
            [add_grid, emp_loss(i, k-model.init_size + 1)] = seq_design(fits(i), all_X(1:k,1:d), all_X(1:k,d+1), model, lhs_rect, i, xtest, xt_dens);
            add_grid = repmat(add_grid(1,:), model.km_batch, 1);  %batch 

            % compute corresponding y-values
            fsim = forward_sim_policy( add_grid,M-i,fits(i:M),model,offset );
            cur_sim = cur_sim + fsim.nsim;

            immPayoff = option_payoff(add_grid, model.K);
            add_mean = mean(fsim.payoff - immPayoff);

            if (strcmp(method, 'probit') || strcmp(method, 'mprobit') )
                add_mean = 2*(add_mean >= 0) - 1;
            end

            add_var = var(fsim.payoff - immPayoff)+0.00001; % avoid unstable results %
            all_X(k+1,:) = [add_grid(1,:), add_mean, add_var, model.km_batch];

            % update fit %
            if (ismember(step, update_kernel_iters) || running_budget == budget)
                fits(i) = gp_setup(model.method, all_X(1:k+1,1:d), all_X(1:k+1,d+1));
            else
                fits(i) = updateFit(fits(i), add_grid(1,:), add_mean);
            end

            % resample the candidate set
            k = k+1;
            running_budget = running_budget + model.km_batch; 
                
            if (k > model.adaptive_grid_loop || running_budget >= budget)
                add_more_sites = false;
            end 
            
            step = step + 1;
        end
    end
    fits(i).t = toc(t_one_pass);
end


fit = fits;
timeElapsed = toc;
nsims = cur_sim;
empLoss = emp_loss;
end