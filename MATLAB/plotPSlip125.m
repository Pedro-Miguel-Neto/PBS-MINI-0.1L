file = 'Partial125100';
load('constants')
% With (true) or without (false) optimization results output
doOptResults2 = false;
% Fit 4 (true) or 2 (false) Sherwood constants
do4const = false;
% Constants [T, P, A, V]
conditions = [21 1 1.635565e-3 9.86746e-5];

% Read file
Table_PSlip125 = readtable(file);

% Table_PSlip125(end,:) = [];

% Compute adimensional numbers
[Table_PSlip125.Re, Table_PSlip125.Sc, Table_PSlip125.G, Table_PSlip125.Sh] = ...
    calcAdimensionalNumbers(Table_PSlip125.T,   Table_PSlip125.P, ...
                            Table_PSlip125.rpm, Table_PSlip125.A, ...
                            Table_PSlip125.V,  Table_PSlip125.kL);

% Fitt Sherwood constants
if doOptResults2 == true
    options = optimset('PlotFcns',@optimplotfval,'TolX', 1e-8, ...
        'MaxIter', 10000, 'MaxFunEvals', 100000);
else
    options = optimset('TolX', 1e-8, 'MaxIter', 10000, ...
        'MaxFunEvals', 100000);
end

% Define objective function and fit constants
if do4const == true
    fun = @(k) sum((Table_PSlip125.Sh - (k(1).*Table_PSlip125.Re.^k(2).*Table_PSlip125.Sc.^k(3).*Table_PSlip125.G.^k(4))).^2);
    k0 = [kinit(1), kinit(2), kinit(3), kinit(4)];
    k = fminsearch(fun, k0, options);
    Table_PSlip125.Shfit = (k(1).*Table_PSlip125.Re.^k(2).*Table_PSlip125.Sc.^k(3).*Table_PSlip125.G.^k(4));
else
    fun = @(k) sum((Table_PSlip125.Sh - (k(1).*Table_PSlip125.Re.^k(2).*Table_PSlip125.Sc.^kinit(3).*Table_PSlip125.G.^kinit(4))).^2);
    k0 = [kinit(1), kinit(2)];
    k = fminsearch(fun, k0, options);
    Table_PSlip125.Shfit = (k(1).*Table_PSlip125.Re.^k(2).*Table_PSlip125.Sc.^kinit(3).*Table_PSlip125.G.^kinit(4));
end

% Oxygen-water diffusion coefficient (m^2/s)
Table_PSlip125.Dab(1:5) = calcDab(Table_PSlip125.T(1:5));

% Compute fitted mass transfer coefficients (m7s)
Table_PSlip125.kLfit(1:5) = Table_PSlip125.Shfit(1:5).*Table_PSlip125.Dab(1:5)./const.D;

% Sort Re and keep the sort index in "sortIdx"
[Table_PSlip125.Re, sortIdx] = sort(Table_PSlip125.Re, "ascend");
    
% Sort table using the sorting index except for Re (Re is already sorted)
for j = 1:width(Table_PSlip125)
    if j ~= find(string(Table_PSlip125.Properties.VariableNames) == 'Re')
        Table_PSlip125(:,j) = Table_PSlip125(sortIdx,j);
    end
end

%Plot experimental kL, fitted kL and power regressions
figure;
hold on;

scatter(Table_PSlip125.rpm, Table_PSlip125.kL, 50, 'MarkerEdgeColor', [0.420    0.420    0.420],   ...
    'DisplayName', 'Sim')

rpm = 0:1:50;
    
[Re, Sc, G, Sh] = ...
    calcAdimensionalNumbers([zeros(1,length(rpm))+conditions(1)]', conditions(2), [rpm]', ...
    conditions(3), [zeros(1,length(rpm))+conditions(4)]', [zeros(1,length(rpm))]');

if do4const == true
    plot(rpm, (k(1).*Re.^k(2).*Sc.^k(3).*G.^k(4)).*calcDab([zeros(1,length(rpm))+conditions(1)]')./const.D, '--', 'linewidth', 2, 'DisplayName', 'Fitted')
    hold on
else
%     plot(rpm, (k(1).*Re.^k(2).*Sc.^kinit(3).*G.^kinit(4)).*calcDab([zeros(1,length(rpm))+conditions(1)]')./const.D, '--', 'linewidth', 2, 'DisplayName', 'Fitted')
    hold on
%     scatter(Table_PSlip125.rpm, (kinit(1).*Table_PSlip125.Re.^kinit(2).*Table_PSlip125.Sc.^kinit(3).*Table_PSlip125.G.^kinit(4)).*calcDab(Table_PSlip125.T)./const.D, 150, 'x', 'LineWidth', 1.5, 'DisplayName', 'Protocol 1')
    scatter(Table_PSlip125.rpm, (kExp(1).*Table_PSlip125.Re.^kExp(2).*Table_PSlip125.Sc.^kinit(3).*Table_PSlip125.G.^kExp(3)).*calcDab(Table_PSlip125.T)./const.D, 75, 'x', 'MarkerEdgeColor', 'red', 'DisplayName', 'Experimental')
end

ax = gca;
ax.FontSize = 13;
title('Slip Coefficient = 12.5%')
legend('boxoff')
legend('Location','northwest')
legend('Orientation','horizontal')
ylabel('k$_L$ (m/s)','interpreter','latex')
xlabel('Agitation Rate (rpm)','interpreter','latex')
% annotation('textbox', [0.185, 0.54, 0.1, 0.1], 'String', ...
%     ['Sh = ' num2str(k(1), 3) 'Re$^{' num2str(k(2), 3) '}$Sc$^{' num2str(kinit(3), 3) '}$G$^{' num2str(kinit(4), 3) '}$'], ...
%     'interpreter', 'latex', 'Edgecolor','none', 'FontSize', 13, ...
%     'Position', [0.2, 0.8, 0.5, 0]);
xlim([0 70])
ylim([0 4e-5])

% % Compute maximum error between correlated kL and simulated kL
% [Err(1,1), Err(2,1)] = max(abs(((k(1).*Table_PSlip125.Re.^k(2).*Table_PSlip125.Sc.^kinit(3) ...
%     .*Table_PSlip125.G.^kinit(4)).*calcDab([zeros(1,length(Table_PSlip125.rpm))+conditions(1)]') ...
%     ./const.D - Table_PSlip125.kL)./Table_PSlip125.kL.*100));
% 
% % Compute maximum error between experimental kL (1st protocol) and simulated kL
% [Err(1,2), Err(2,2)] = max(abs(((kinit(1).*Table_PSlip125.Re.^kinit(2).*Table_PSlip125.Sc.^kinit(3) ...
%     .*Table_PSlip125.G.^kinit(4)).*calcDab([zeros(1,length(Table_PSlip125.rpm))+conditions(1)]') ...
%     ./const.D - Table_PSlip125.kL)./((kinit(1).*Table_PSlip125.Re.^kinit(2).*Table_PSlip125.Sc.^kinit(3) ...
%     .*Table_PSlip125.G.^kinit(4)).*calcDab([zeros(1,length(Table_PSlip125.rpm))+conditions(1)]') ...
%     ./const.D).*100));
% 
% % Compute maximum error between experimental kL (2nd protocol) and simulated kL
% [Err(1,3), Err(2,3)] = max(abs(((kExp(1).*Table_PSlip125.Re.^kExp(2).*Table_PSlip125.Sc.^kinit(3) ...
%     .*Table_PSlip125.G.^kExp(3)).*calcDab([zeros(1,length(Table_PSlip125.rpm))+conditions(1)]') ...
%     ./const.D - Table_PSlip125.kL)./((kExp(1).*Table_PSlip125.Re.^kExp(2).*Table_PSlip125.Sc.^kinit(3) ...
%     .*Table_PSlip125.G.^kExp(3)).*calcDab([zeros(1,length(Table_PSlip125.rpm))+conditions(1)]') ...
%     ./const.D).*100));
% 
% % Compute mean error between correlated kL and simulated kL
% Err(3,1) = mean(abs(((k(1).*Table_PSlip125.Re.^k(2).*Table_PSlip125.Sc.^kinit(3) ...
%     .*Table_PSlip125.G.^kinit(4)).*calcDab([zeros(1,length(Table_PSlip125.rpm))+conditions(1)]') ...
%     ./const.D - Table_PSlip125.kL)./Table_PSlip125.kL.*100));
% 
% % Compute mean error between experimental kL (1st protocol) and simulated kL
% Err(3,2) = mean(abs(((kinit(1).*Table_PSlip125.Re.^kinit(2).*Table_PSlip125.Sc.^kinit(3) ...
%     .*Table_PSlip125.G.^kinit(4)).*calcDab([zeros(1,length(Table_PSlip125.rpm))+conditions(1)]') ...
%     ./const.D - Table_PSlip125.kL)./((kinit(1).*Table_PSlip125.Re.^kinit(2).*Table_PSlip125.Sc.^kinit(3) ...
%     .*Table_PSlip125.G.^kinit(4)).*calcDab([zeros(1,length(Table_PSlip125.rpm))+conditions(1)]') ...
%     ./const.D).*100));
% 
% % Compute mean error between experimental kL (2nd protocol) and simulated kL
% Err(3,3) = mean(abs(((kExp(1).*Table_PSlip125.Re.^kExp(2).*Table_PSlip125.Sc.^kinit(3) ...
%     .*Table_PSlip125.G.^kExp(3)).*calcDab([zeros(1,length(Table_PSlip125.rpm))+conditions(1)]') ...
%     ./const.D - Table_PSlip125.kL)./((kExp(1).*Table_PSlip125.Re.^kExp(2).*Table_PSlip125.Sc.^kinit(3) ...
%     .*Table_PSlip125.G.^kExp(3)).*calcDab([zeros(1,length(Table_PSlip125.rpm))+conditions(1)]') ...
%     ./const.D).*100));
% 
% Err(2,:) = Table_PSlip125.rpm(Err(2,:))