clear all;

currentFolder = pwd;
parentFolder = fileparts(currentFolder);
addpath(parentFolder);
resultsFolder = fullfile(parentFolder,'results/testBatch/am_put_3d_r30_n1000');

n = 20;

% load(strcat(resultsFolder, '/gauss_lhs_plain'));
% put_1 = put3d;
% m = size(put_1,1);
% t_1 = zeros(m-1, n);
% for i = 1:n
%     for j = 1:m-1
%         if (~isempty(put_1(j, i).t))
%             t_1(j, i) = put_1(j, i).t;
%         end
%     end
% end
% 
% m_t_1 = mean(mean(t_1));
% sd_t_1 = std(mean(t_1));
% 
% n_1 = zeros(m-1, n);
% for i = 1:n
%     for j = 1:m-1
%         if (~isempty(put_1(j, i).t))
%             n_1(j, i) = size(put_1(j, i).r,1);
%         end
%     end
% end
% 
% m_n_1 = mean(mean(n_1));
% sd_n_1 = std(mean(n_1));
% 
% m_pay_1 = mean(mean(payoff_gp));
% sd_pay_1 = std(mean(payoff_gp));
% 
% [m_pay_1 sd_pay_1]
% [m_t_1 sd_t_1 m_n_1 sd_n_1]

load(strcat(resultsFolder, '/gauss_MCU_FB'));
put_2 = put3d;
m = size(put_2,1);
t_2 = zeros(m-1, n);
for i = 1:n
    for j = 1:m-1
        if (~isempty(put_2(j, i).t))
            t_2(j, i) = put_2(j, i).t;
        end
    end
end

m_t_2 = mean(mean(t_2));
sd_t_2 = std(mean(t_2));

n_2 = zeros(m-1, n);
for i = 1:n
    for j = 1:m-1
        if (~isempty(put_2(j, i).t))
            n_2(j, i) = size(put_2(j, i).r,1);
        end
    end
end

m_n_2 = mean(mean(n_2));
sd_n_2 = std(mean(n_2));

m_pay_2 = mean(mean(payoff_gp));
sd_pay_2 = std(mean(payoff_gp));

[m_pay_2 sd_pay_2]
[m_t_2 sd_t_2 m_n_2 sd_n_2]

load(strcat(resultsFolder, '/gauss_MCU_RB'));
put_3 = put3d;
m = size(put_3,1);
t_3 = zeros(m-1, n);
for i = 1:n
    for j = 1:m-1
        if (~isempty(put_3(j, i).t))
            t_3(j, i) = put_3(j, i).t;
        end
    end
end

m_t_3 = mean(mean(t_3));
sd_t_3 = std(mean(t_3));

n_3 = zeros(m-1, n);
for i = 1:n
    for j = 1:m-1
        if (~isempty(put_3(j, i).t))
            n_3(j, i) = size(put_3(j, i).r,1);
        end
    end
end

m_n_3 = mean(mean(n_3));
sd_n_3 = std(mean(n_3));

m_pay_3 = mean(mean(payoff_gp));
sd_pay_3 = std(mean(payoff_gp));

[m_pay_3 sd_pay_3]
[m_t_3 sd_t_3 m_n_3 sd_n_3]

load(strcat(resultsFolder, '/gauss_MCU_MLB'));
put_4 = put3d;
m = size(put_4,1);
t_4 = zeros(m-1, n);
for i = 1:n
    for j = 1:m-1
        if (~isempty(put_4(j, i).t))
            t_4(j, i) = put_4(j, i).t;
        end
    end
end

m_t_4 = mean(mean(t_4));
sd_t_4 = std(mean(t_4));

n_4 = zeros(m-1, n);
for i = 1:n
    for j = 1:m-1
        if (~isempty(put_4(j, i).t))
            n_4(j, i) = size(put_4(j, i).r,1);
        end
    end
end

m_n_4 = mean(mean(n_4));
sd_n_4 = std(mean(n_4));

m_pay_4 = mean(mean(payoff_gp));
sd_pay_4 = std(mean(payoff_gp));

[m_pay_4 sd_pay_4] 
[m_t_4 sd_t_4 m_n_4 sd_n_4]

load(strcat(resultsFolder, '/gauss_ABSUR_ABSUR'));
put_5 = put3d;
m = size(put_5,1);
t_5 = zeros(m-1, n);
for i = 1:n
    for j = 1:m-1
        if (~isempty(put_5(j, i).t))
            t_5(j, i) = put_5(j, i).t;
        end
    end
end

m_t_5 = mean(mean(t_5));
sd_t_5 = std(mean(t_5));

n_5 = zeros(m-1, n);
for i = 1:n
    for j = 1:m-1
        if (~isempty(put_5(j, i).t))
            n_5(j, i) = size(put_5(j, i).r,1);
        end
    end
end

m_n_5 = mean(mean(n_5));
sd_n_5 = std(mean(n_5));

m_pay_5 = mean(mean(payoff_gp));
sd_pay_5 = std(mean(payoff_gp));

[m_pay_5 sd_pay_5] 
[m_t_5 sd_t_5 m_n_5 sd_n_5]

load(strcat(resultsFolder, '/gauss_MCU_ADSA'));
put_6 = put3d;
m = size(put_6,1);
t_6 = zeros(m-1, n);
for i = 1:n
    for j = 1:m-1
        if (~isempty(put_6(j, i).t))
            t_6(j, i) = put_6(j, i).t;
        end
    end
end

m_t_6 = mean(mean(t_6));
sd_t_6 = std(mean(t_6));

n_6 = zeros(m-1, n);
for i = 1:n
    for j = 1:m-1
        if (~isempty(put_6(j, i).t))
            n_6(j, i) = size(put_6(j, i).r,1);
        end
    end
end

m_n_6 = mean(mean(n_6));
sd_n_6 = std(mean(n_6));

m_pay_6 = mean(mean(payoff_gp));
sd_pay_6 = std(mean(payoff_gp));

[m_pay_6 sd_pay_6] 
[m_t_6 sd_t_6 m_n_6 sd_n_6]

load(strcat(resultsFolder, '/gauss_MCU_DDSA'));
put_7 = put3d;
m = size(put_7,1);
t_7 = zeros(m-1, n);
for i = 1:n
    for j = 1:m-1
        if (~isempty(put_7(j, i).t))
            t_7(j, i) = put_7(j, i).t;
        end
    end
end

m_t_7 = mean(mean(t_7));
sd_t_7 = std(mean(t_7));

n_7 = zeros(m-1, n);
for i = 1:n
    for j = 1:m-1
        if (~isempty(put_7(j, i).t))
            n_7(j, i) = size(put_7(j, i).r,1);
        end
    end
end

m_n_7 = mean(mean(n_7));
sd_n_7 = std(mean(n_7));

m_pay_7 = mean(mean(payoff_gp));
sd_pay_7 = std(mean(payoff_gp));

[m_pay_7 sd_pay_7] 
[m_t_7 sd_t_7 m_n_7 sd_n_7]

%%
figure (1);
[ax1, c1] = plot_amput(put_5(15, 15), 1, 2, 1);
[ax2, c2, h, pos] = plot_amput(put_6(15, 19), 1, 2, 2);
c = [min([c1 c2]), max([c1 c2])];
caxis(c)
set(ax2,'YTick', [], 'YLabel' ,[]);
subplot(1, 2, 1)
colorbar off
% set(ax1, 'Units', 'Normalized', 'OuterPosition', [0.10, 0.10, 0.4, 0.7]);
ps = ax1.Position;
% set(ax1,'YTick', [],'YLabel' ,[],  'Position', [ps(1), ps(2), ps(3), 0.8*ps(4)]);
set(ax2,'YTick', [],'YLabel' ,[],  'Position', [ps(1) + 0.4, ps(2), ps(3) + 0.1, ps(4)]);
% set(ax3,'XTick',[], 'YTick', [], 'XLabel' ,[], 'YLabel' ,[],  'OuterPosition', [0.61, 0.48, 0.30, 0.5]);
% set(ax4, 'YTick', [0, 0.5, 1], 'OuterPosition', [0.01, 0.01, 0.32, 0.5]);
% set(ax5,'YTick', [], 'YLabel' ,[], 'OuterPosition', [0.33, 0.01, 0.28, 0.5]);
% set(ax6,'YTick', [], 'YLabel' ,[], 'OuterPosition', [0.61, 0.01, 0.30, 0.5]);
% h.Position = [0.91 pos(2) 0.9*pos(3) pos(4)];
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% xlabel('$x_1$', 'Interpreter','latex');
% ylabel('$x_2$', 'Interpreter','latex');
saveas(gcf, fullfile(resultsFolder,'fitted_all_am_put.png'))

%%
[c1] = plot_2d_amput_batch(put_5(15, 15));
saveas(gcf, fullfile(resultsFolder,'fitted_absur_am_put.pdf'))
%%
[c3] = plot_2d_amput_batch(put_6(15, 19));
saveas(gcf, fullfile(resultsFolder,'fitted_adsa_am_put.pdf'))