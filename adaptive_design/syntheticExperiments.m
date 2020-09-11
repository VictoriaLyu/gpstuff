% SYNTHETICEXPERIMENTS the main file used to generate results for synthetic
% experiments: 1d Exponential, 2d Brannin-Hoo, 2d Michal, 2d GP, 2d TP and 
% 6d Hartman.

% Description:
%   SYNTHETICEXPERIMENTS calculates the synthetic results for known true
%   functions. There are four cases for the noise structure: normal
%   distributed noise with small/large variance, t distributed noise with
%   small/large variance, mixture normal distributed noise and
%   heteroskedastic t distributed noise. GP, t-GP and Cl-GP are supported
%   likelihood for Gaussian Process simulator. The parameters are set at
%   the same value as in the paper (also support user defined parameters).
%   Here is a demo for 2d Branin-Hoo function.

%%% k is the number of runs 
%%% I is the initial design
%%% budget is the total budget
%%% m0 is the test size

clear all; close all; clc;

currentFolder = pwd;
parentFolder = fileparts(currentFolder);
addpath(parentFolder);

startup;


%%%%%% dimensions %%%%%%

% %%%%% 1d %%%%%
% 
% k = 100; % runs of experiments
% I = 10; 
% d = 1;
% budget = 80;
% m0 = 1000;
% fun = @(x) (x+0.75).*(x-0.75); 

%%%%%%%%%%%%%%%%%%%%%%% other choices %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 2d %%%%%

k = 100; % runs of experiments
I = 20;
d = 2;
budget = 150;
m0 = 50;
fun = @braninsc2;
% fun = @michal;
% fun = @gp;
% fun = @tp;

%%%% 6d %%%%%
% 
% k = 100; % runs of experiments
% I = 60;
% d = 6;
% budget = 1000;
% m0 = 0;
% fun = @hart6; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% noise cases %%%

noisestructure = {'t_constdf', 't_constdf', 'mixed', 't_heterodf'}; 
noisevar = {'small', 'large', 'mixed', 'hetero'}; 

% %%% noise cases for GP/TP true function %%%
% 
% noisestructure = {'normal', 'normal', 't_constdf', 't_constdf'}; 
% noisevar = {'small', 'large', 'small', 'large'}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% model %%%

model = {'gauss', 't', 'probit'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% design %%%

design = {'cUCB', 'tMSE', 'gSUR', 'SUR'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(noisestructure,2)
    
    if (isequal(fun, @gp) || isequal(fun, @tp))
        fun = initializeGPTP(d, 40, fun, char(noisestructure(i)), char(noisevar(i)));
    end
    
    for j = 1:size(design,2)
        for jj = 1:size(model, 2)
            [x_seq, y_seq, ee, er, bias, metric, t, Ef, Varf] = updategppar(I, k, m0, d, budget, fun, char(noisestructure(i)), char(noisevar(i)), char(model(jj)), char(design(j)));
        end
    end
end
