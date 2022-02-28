folder = {'probes_Refined', 'probes_pSlip125'};

opts = load('opts.mat');

figure;
hold on

for n = 1:2
    
    file = [char(folder(n)) '/0/kcInstF'];
    Table = readtable(file, opts.opts);
    Table(1:2,:) = [];
    Table(:,3:end) = [];
    
    if n == 1
        plot(Table.Var1, Table.Var2, 'color', [1.0000 0.4118 0.1608], 'linewidth', 3, 'DisplayName', 'Refined Mesh')
    elseif n == 2
        plot(Table.Var1, Table.Var2, '--', 'color', 'black', 'linewidth', 3, 'DisplayName', 'Coarse Mesh')
    end
end

% y-Axis
xrange = [0 7.5];
xtick = 0.5;
% y-Axis
yrange = [0 inf];
% ytick = 1e-5;
% ------------------------------------------- Tick range, spacing and angle
% x-Axis
xtl=(xrange(1):xtick:xrange(2));
xticks(xtl);
xtickangle(0);
% y-Axis
% ytl=(yrange(1):ytick:9e-5);
% yticks(ytl);
ytickangle(0);
% ------------------------------------------------------ Label descriptions
% x-Axis
xlabel('Time (s)','interpreter','latex')
% y-Axis
ylabel('k$_{L}$ (m/s)','interpreter','latex')
% ------------------------------------------------ Figure Position and Size
%[HorzPs, VertPs, HorzSz, VertSz]
set(gcf,'Position',[10 350 1900 500]);
% ------------------------------------------------------ Box configurations
box off
% ----------------------------------------------------- Axis configurations
axis tight
ax = gca;
ax.FontSize = 13;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridColor = [0.1490    0.1490    0.1490];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
% --------------------------------------------------- Legend configurations
legend('boxoff')
legend('Location','northeast')
legend('Orientation','Vertical')
% title(legend,'RPM')
title('Instantaneous Mass transfer Coefficient')
% --------------------------------------------------------------- Color map
map = brewermap(7,'Set1');
xlim(xrange);
ylim(yrange);