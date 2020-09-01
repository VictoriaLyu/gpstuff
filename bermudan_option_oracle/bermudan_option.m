clear all; close all; clc;

currentFolder = pwd;
parentFolder = fileparts(currentFolder);
addpath(parentFolder);
addpath(fullfile(parentFolder, 'others'));
addpath(fullfile(parentFolder, 'adaptive_design'));

% % 
startup;

n = 10;

method = {'gauss', 't', 'probit', 'mgauss', 'mprobit'};
design = {'MCU', 'tMSE', 'cSUR', 'ICU'};

%%%%%% 2D bermudan put %%%%%%
model.km_batch = 15;

model.budget = 1200;
model.adaptive_grid_loop = 80;

model.look_ahead = 1;
model.init_size = 10;
model.final_runs = 0;

model.cand_len = 1000;
model.K = 40;
model.x0 = [40,40];
model.sigma = [0.2,0.2];
model.r = 0.06;
model.div = 0;
model.T = 1;
model.dt = 0.04;
model.dim = 2;
model.sim_func = @sim_gbm;
model.option_payoff = @put_payoff;
model.search = false;
model.el_thresh = 0.0001;

rng(1);

lhs20 = lhsdesign(20,2);
lhs20 = 25 + 30*lhs20;
lhs20 = lhs20(lhs20(:,1) + lhs20(:,2) <= 80, :);

model.init_grid = lhs20;

model.lhs_rect = repmat([25,55], 2, 1);

rng(1);

NN = 160000;
MM = 25;

mygr = zeros(NN,2,MM+1);

mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);

for i = 2:(MM+1)
   mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
end

M = model.T / model.dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% 3D Max Call %%%%%%%
% 
% model.km_batch = 20;
% model.budget = 2000;
% model.adaptive_grid_loop = 100;
% 
% model.look_ahead = 1;
% model.init_size = 30;
% model.final_runs = 0;
% 
% model.cand_len = 1000;
% model.K = 100;
% model.x0 = [90,90,90];
% model.sigma = [0.2,0.2,0.2];
% model.r = 0.05;
% model.div = 0.1;
% model.T = 3;
% model.dt = 1/3;
% model.dim = 3;
% model.sim_func = @sim_gbm;
% model.option_payoff = @maxCall;
% model.search = false;
% model.el_thresh = 0.0001;
% 
% rng(6);
% 
% lhs30 = lhsdesign(35,3);
% lhs30 = 50 + 100*lhs30;
% lhs30 = lhs30(max([lhs30(:,1), lhs30(:,2), lhs30(:,3)], [], 2) >= 100, :);
% 
% model.init_grid = lhs30;
% 
% model.lhs_rect = repmat([50,150], 3, 1);
% 
% NN = 160000;
% MM = 9;
% 
% M = model.T / model.dt;
% 
% rng(10);
%     
% mygr = zeros(NN,3,MM);
% 
% mygr(:,:,1) = model.sim_func( repmat(model.x0, NN, 1), model, model.dt);
% 
% for i = 2:MM
%    mygr(:,:,i) = model.sim_func( mygr(:,:,i-1), model, model.dt);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i =  1:size(method, 2)
    for k = 1:size(design, 2)
        model.design = char(design(k));
        model.method = char(method(i));
        
        fits = repmat(struct('gp',[],'x',[],'y',[], 't', []), M, n);
        timeElapsed = zeros(n,1);
        nsims = zeros(n,1);
        empLoss = zeros(M, model.adaptive_grid_loop - model.init_size, n);
        payoff_gp = zeros(NN,n);

        parfor j = 1:n
            [fits(:,j), timeElapsed(j), nsims(j), empLoss(:,:,j)] = osp_seq_design(model);
            oos = forward_sim_policy( mygr, MM, fits(:,j), model);
            payoff_gp(:,j) = oos.payoff;
        end
        
        [mean(mean(payoff_gp)) std(mean(payoff_gp))]
        plot_2d_put(fits(15))
    end
end
