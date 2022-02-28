TableDab = readtable('2PDabPlot');
% T - Empirical temperatures (ºC)
% Dab - Empirical oxygen-water diffusion coefficient (m^2/s)
% Fit Dab data:
PolyFitDab = polyfit(TableDab.T, TableDab.Dab*10^-4, 1);

figure;
scatter(TableDab.T, TableDab.Dab*10^-4, 'filled')
hold on
T = 10:1:30;
plot(T, PolyFitDab(1).*T+PolyFitDab(2), '--', 'linewidth', 2, 'Color', [.620 .620 .620])

% box off;
% ----------------------------------- Axis and axis grid configurations
axis tight;
ax = gca;
ax.FontSize = 13;
ax.XGrid = 'off';
ax.YGrid = 'off';
ax.GridColor = [0.1490    0.1490    0.1490];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
% ----------------------------------------------- Legend configurations
% legend('boxoff');
% legend('Position', [0.5, 0.54, 0.2, 0]);
% legend('Orientation','horizontal');
%     title(legend,'pO2');
annotation('textbox', [0.185, 0.54, 0.1, 0.1], 'String', ...
    ['D$_{AB}$ = ' num2str(PolyFitDab(1),3) '$\times$T' '+' num2str(PolyFitDab(2),3)], 'interpreter', 'latex', ...
    'Edgecolor','none', 'FontSize', 14, 'Position', [0.18, 0.8, 0.5, 0]);
% ---------------------------------------------------------- Axis range
% y-Axis
xrange = [0 inf];
xlim([10 inf]);
% y-Axis
yrange = [0 inf];
ylim([1.5e-9 inf]);
% --------------------------------------- Tick range, spacing and angle
% x-Axis
xtickangle(0);
% y-Axis
ytickangle(0);
title('Oxygen-water Diffusivity Fit')
% x-Axis
xlabel('Time (s)','interpreter','latex');
% y-Axis
ylabel('D$_{AB}$ (m$^2$/s)','interpreter','latex');