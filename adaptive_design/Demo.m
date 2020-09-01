clear all; close all; clc;

currentFolder = pwd;
parentFolder = fileparts(currentFolder);
addpath(parentFolder);

startup;

% 2d experiments%
I = 20;
d = 2;
budget = 150;
m0 = 50;
fun = @braninsc2;

rng(0);
Xint = lhsdesign(I, 2);

x1=repmat(linspace(0,1,m0)',1,m0);
x2=repmat(linspace(0,1,m0)',1,m0)';
xt=[x1(:) x2(:)];
ft = fun(xt);

model = 'gauss';
design = 'tMSE';

noisestructure = 't_constdf';
noisevar = 'small';

[x_seq, y_seq, metric, time, ee, er, bias, Ef, Varf, l, sigma2, sigman, t_optim, t_gen] = updategp(fun, noisestructure, noisevar, Xint, xt, model, design, budget, ft);

m = Ef(:, end);
s = abs(Varf(:, end));

figure(1), hold on, box on
contour(reshape(xt(:,1),m0,m0),reshape(xt(:,2),m0,m0),reshape(ft,m0,m0),[0 0], 'k-', 'linewidth', 2.5);
contour(reshape(xt(:,1),m0,m0),reshape(xt(:,2),m0,m0),reshape(m,m0,m0),[0 0], 'r-', 'linewidth', 2.5);
contour(reshape(xt(:,1),m0,m0),reshape(xt(:,2),m0,m0),reshape(m+2*sqrt(s),m0,m0),[0 0],'r--', 'linewidth', 2.5);
contour(reshape(xt(:,1),m0,m0),reshape(xt(:,2),m0,m0),reshape(m-2*sqrt(s),m0,m0),[0 0], 'r--', 'linewidth', 2.5);
scatter(x_seq(:,1), x_seq(:,2), 'b', 'filled');
xlabel('$x_1$', 'Interpreter','latex');
ylabel('$x_2$', 'Interpreter','latex');
set(gca, 'FontSize', 18);
