function plot_2d_amcall(fit)

% Plots the fitted Ameriacal call option at one time stamp

xt1=repmat(linspace(50,150,50)',1,50);
xt2=repmat(linspace(50,150,50)',1,50)';
xt=[xt1(:) xt2(:)];

[m, s] = gp_pred(fit.gp, fit.x, fit.y, xt);

figure, hold on;
h1=pcolor(reshape(xt(:,1),50,50),reshape(xt(:,2),50,50),reshape(m,50,50));
set(h1, 'edgealpha', 0), set(h1, 'facecolor', 'interp');
contour(reshape(xt(:,1),50,50),reshape(xt(:,2),50,50),reshape(m,50,50),[0 0], 'k-', 'linewidth', 3);
contour(reshape(xt(:,1),50,50),reshape(xt(:,2),50,50),reshape(m+2*sqrt(s),50,50),[0 0],'k--', 'linewidth', 3);
contour(reshape(xt(:,1),50,50),reshape(xt(:,2),50,50),reshape(m-2*sqrt(s),50,50),[0 0], 'k--', 'linewidth', 3);
plot(fit.x(fit.y <= 0,1),fit.x(fit.y <= 0,2),'ro', 'markersize', 8, 'linewidth', 2);
plot(fit.x(fit.y > 0,1),fit.x(fit.y > 0,2),'rx', 'markersize', 8, 'linewidth', 2);
xlabel('$x_1$', 'Interpreter','latex');
ylabel('$x_2$', 'Interpreter','latex');
colorbar
set(gca, 'FontSize', 18);