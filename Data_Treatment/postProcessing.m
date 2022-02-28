load('constants')
%% PART 1 - SIMULATION VALUES

%% Initial inputs

folder = {'Slip100', 'Slip80', 'Slip60'};
% folder = {'Partial20100', 'Partial2080', 'Partial2060'};

% With (true) or without (false) optimization results output
doOptResults = false;
% Do (true) or do not (false) show luminescence plots
doPlots = false;
% Define optimization bounds (min,max)
kLbounds = [1e-6 1e-4];
% Time from which kL evaluation starts
t0 = 4;

% Constants [T, P, A, V]
conditions = [21 1 1.635565e-3 9.86746e-5; 21 1 1.60185e-3 7.997321e-5; 37 1 1.56514e-3 5.996167e-5];
conditionsExp = [24.3125 1 1.643433e-3 0.000103; 24.5500 1 1.610688e-3 0.000085; 24.4000 1 1.56514e-3 0.00006];

% With (true) or without (false) optimization results output
doOptResults2 = false;
% Fit 4 (true) or 2 (false) Sherwood constants
do4const = false;
% Initial estimates [k, Re, Sc, G]
kinit = [1.6784 0.5508 1/3 1.1817];
% Set range (L) (ex: 100mL +- 5mL)
setRange = 1e-5;

%% Receives a cell list with folders to analyse, the transport conditions for the cases and the t0
% Constructs a struct with all input data and computes and stores kInstF and kcF
for n = 1:length(folder)
    
    % Initiate struct/substruct
    Table_kcInstF.(char(folder(n))) = table2struct(getFileNumber(char(folder(n))));
    % Declare figure
    if n == 1
        figure;
        hold on;
    end
    % Readtable Delimiter Options
    opts = load('opts.mat');
    
    % Computes and stores kcInstF and kcF
    for k = 1:length([Table_kcInstF.(char(folder(n))).rpm])
        % Add conditions to struct
        Table_kcInstF.(char(folder(n)))(k).T = conditions(n,1);
        Table_kcInstF.(char(folder(n)))(k).P = conditions(n,2);
        Table_kcInstF.(char(folder(n)))(k).A = conditions(n,3);
        Table_kcInstF.(char(folder(n)))(k).V = conditions(n,4);
        % Assign path
        path = [char(folder(n)) '/' char(Table_kcInstF.(char(folder(n)))(k).name)];
        % Create path variable
        S = dir(path);
        % Compute number of path subfolders
        numfolders = sum([S.isdir]);
        % Initialize table variables
        t = [];
        kcInstF = [];
        % Concatenate variable values from all subfloders into one list
        for m = (3:numfolders)
            file = [path '/' S(m).name '/kcInstF'];
            trial = readtable(file, opts.opts);
            trial(1,:) = [];
            t = [t; trial.Var1];
            kcInstF = [kcInstF; trial.Var2];
        end
        % Plot kcInstF results for 100 mL
        if n == 1
            plot(t,kcInstF, 'DisplayName', num2str(Table_kcInstF.(char(folder(n)))(k).rpm, 3), 'linewidth', 1);  %For 2-D plot
        end
        % Compute time-averaged kL from starting time t0
        kcF = zeros(length(kcInstF),1);
        r = find(t == t0);
        kcF(r) = kcInstF(r);
        for m = r+1:length(kcInstF)
            kcF(m) = (kcF(m-1).*(t(m-1)-t(1))+kcInstF(m).*(t(m)-t(m-1)))./(t(m)-t(1));
        end
        % Store all data in struct
        Table_kcInstF.(char(folder(n)))(k).trial = table(t, kcInstF, kcF);
        Table_kcInstF.(char(folder(n)))(k).kL =  Table_kcInstF.(char(folder(n)))(k).trial.kcF(end);
        
    end
    
    %% Figure settings
    % -------------------------------------------------------------- Axis range
    % y-Axis
    xrange = [0 60];
    xlim(xrange);
    xtick = 10;
    % y-Axis
    yrange = [0 inf];
    ylim(yrange);
    ytick = 1e-5;
    % ------------------------------------------- Tick range, spacing and angle
    % x-Axis
    xtl=(xrange(1):xtick:xrange(2));
    xticks(xtl);
    xtickangle(0);
    % y-Axis
    ytl=(yrange(1):ytick:9e-5);
    yticks(ytl);
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
    legend('Orientation','horizontal')
    title(legend,'RPM')
    title('Instantaneous Mass transfer Coefficient')
    % --------------------------------------------------------------- Color map
    map = brewermap(7,'Set1');
    xlim([0 60])
end

%% Receives a struct and merges substructs to form a table
% Construct table fom struct
Table_kcInstFSh = struct2table(Table_kcInstF.(char(folder(1))));   
for n = 2:length(folder)
    Table_kcInstFSh.rpm(end+1) = Table_kcInstF.(char(folder(n))).rpm;
    Table_kcInstFSh.T(end) = Table_kcInstF.(char(folder(n))).T;
    Table_kcInstFSh.P(end) = Table_kcInstF.(char(folder(n))).P;
    Table_kcInstFSh.A(end) = Table_kcInstF.(char(folder(n))).A;
    Table_kcInstFSh.V(end) = Table_kcInstF.(char(folder(n))).V;
    Table_kcInstFSh.kL(end) = Table_kcInstF.(char(folder(n))).kL;
end
% Remove unwanted fields
Table_kcInstFSh = removevars(Table_kcInstFSh,{'name', 'trial'});

%% Computes Reynolds, Schmidt, Geometrical term and Sherwood
[Table_kcInstFSh.Re, Table_kcInstFSh.Sc, Table_kcInstFSh.G, Table_kcInstFSh.Sh] = ...
    calcAdimensionalNumbers(Table_kcInstFSh.T, Table_kcInstFSh.P, Table_kcInstFSh.rpm, ...
    Table_kcInstFSh.A, Table_kcInstFSh.V, Table_kcInstFSh.kL);

%% Fits Sherwood constants
% With (true) or without (false) optimization results output
if doOptResults2 == true
    options = optimset('PlotFcns',@optimplotfval,'TolX', 1e-8, ...
        'MaxIter', 10000, 'MaxFunEvals', 100000);
else
    options = optimset('TolX', 1e-8, 'MaxIter', 10000, ...
        'MaxFunEvals', 100000);
end

% Define objective function and fit constants
if do4const == true
    fun = @(k) sum((Table_kcInstFSh.Sh - (k(1).*Table_kcInstFSh.Re.^k(2).*Table_kcInstFSh.Sc.^k(3).*Table_kcInstFSh.G.^k(4))).^2);
    k0 = [kinit(1), kinit(2), kinit(3), kinit(4)];
    k = fminsearch(fun, k0, options)
    Table_kcInstFSh.Shfit = (k(1).*Table_kcInstFSh.Re.^k(2).*Table_kcInstFSh.Sc.^k(3).*Table_kcInstFSh.G.^k(4));
else
    fun = @(k) sum((Table_kcInstFSh.Sh - (k(1).*Table_kcInstFSh.Re.^k(2).*Table_kcInstFSh.Sc.^kinit(3).*Table_kcInstFSh.G.^k(3))).^2);
    k0 = [kinit(1), kinit(2), kinit(4)];
    k = fminsearch(fun, k0, options)
    Table_kcInstFSh.Shfit = (k(1).*Table_kcInstFSh.Re.^k(2).*Table_kcInstFSh.Sc.^kinit(3).*Table_kcInstFSh.G.^k(3));
end

%% Computes Dab and kLfit for the table

%Oxygen-water diffusion coefficient (m^2/s)
Table_kcInstFSh.Dab = calcDab(Table_kcInstFSh.T);
% Compute fitted mass transfer coefficients (m/s)
Table_kcInstFSh.kLfit = Table_kcInstFSh.Shfit.*Table_kcInstFSh.Dab./const.D;

%% Sorts the table based on Re to enable plots

% Sort Re and keep the sort index in "sortIdx"
[Table_kcInstFSh.Re, sortIdx] = sort(Table_kcInstFSh.Re, "ascend");

% Sort table using the sorting index except for Re (Re is already sorted)
for j = 1:width(Table_kcInstFSh)
    if j ~= find(string(Table_kcInstFSh.Properties.VariableNames) == 'Re')
        Table_kcInstFSh(:,j) = Table_kcInstFSh(sortIdx,j);
    end
end

%% Plots simulations' kL and fitted kL
figure;
hold on;

% Simulations'
scatter(Table_kcInstFSh.rpm, Table_kcInstFSh.kL, 50, 'MarkerEdgeColor', [0.420    0.420    0.420],   ...
    'DisplayName', 'Sim.')

% Compute fitted kL for desired range and conditions
rpm = 0:1:120;
for n = 1:length(folder)    
    [Re, Sc, G, Sh] = ...
        calcAdimensionalNumbers([zeros(1,length(rpm))+conditions(n,1)]', conditions(n,2), [rpm]', ...
        conditions(n,3), [zeros(1,length(rpm))+conditions(n,4)]', [zeros(1,length(rpm))]');    
    if do4const == true
        plot(rpm, (k(1).*Re.^k(2).*Sc.^k(3).*G.^k(4)).*calcDab([zeros(1,length(rpm))+conditions(n,1)]')./const.D, '--', 'linewidth', 2, 'DisplayName', 'Sim. Corr.')
        hold on
    else
        if n == 1
            plot(rpm, (k(1).*Re.^k(2).*Sc.^kinit(3).*G.^k(3)).*calcDab([zeros(1,length(rpm))+conditions(n,1)]')./const.D, '--', 'linewidth', 1, 'Color', [0    0.4471    0.7412], 'DisplayName', 'Sim. Corr.')
        hold on
        else
            plot(rpm, (k(1).*Re.^k(2).*Sc.^kinit(3).*G.^k(3)).*calcDab([zeros(1,length(rpm))+conditions(n,1)]')./const.D, '--', 'linewidth', 1, 'Color', [0    0.4471    0.7412], 'HandleVisibility', 'off')    
        end
    end
end

%% Figure settings

ax = gca;
ax.FontSize = 13;
title('Simulations')
legend('boxoff')
legend('Location','northwest')
legend('Orientation','horizontal')
ylabel('k$_L$ (m/s)','interpreter','latex')
xlabel('Agitation Rate (rpm)','interpreter','latex')
annotation('textbox', [0.185, 0.54, 0.1, 0.1], 'String', ...
    ['Sh = ' num2str(k(1), 3) 'Re$^{' num2str(k(2), 3) '}$Sc$^{' num2str(kinit(3), 3) '}$G$^{' num2str(k(3), 3) '}$'], ...
    'interpreter', 'latex', 'Edgecolor','none', 'FontSize', 13, ...
    'Position', [0.2, 0.8, 0.5, 0]);
xlim([0 120])

%% PART 2 - EXPERIMENTAL VALUES

%% Initial inputs

% With (true) or without (false) optimization results output
doOptResults = false;
% Do (true) or do not (false) show luminescence plots
doPlots = false;
% Define optimization bounds (min,max)
kLbounds = [1e-6 1e-4];

% Constants
P = 1;          % Pressure (atm)

% With (true) or without (false) optimization results output
doOptResults2 = false;
% Fit 4 (true) or 2 (false) Sherwood constants
do4const2 = false;
% Set range (L) (ex: 100mL +- 5mL)
setRange = 1e-5;

%% Reads experimental files and fits a kL for each one

% Retrieve trial conditions and analysis range
Table_I = cell2struct([struct2cell(table2struct ...
    (readtable('P2Dict.txt')));struct2cell(table2struct ...
    (readtable('P2Range.txt')))],[fieldnames(table2struct ...
    (readtable('P2Dict.txt')));fieldnames(table2struct ...
    (readtable('P2Range.txt')))]);


for n = 1:length({Table_I.date})
    % Retrieve kcInstF data from files
    Table_I(n).table = readtable(num2str(n));
    % Compute luminescence deltas
    Table_I(n).table.Iexp = Table_I(n).I0 - Table_I(n).table.Imed;
    % Define search range from range fields
    eval(['search = [' char(Table_I(n).range1) ', ' char(Table_I(n).range2) '];']);
    % Define initial time for kL regression
    if find(Table_I(n).table.Iexp >= 0, 1) == 1
        t0 = Table_I(n).table.t(1);
    else
        t0 = Table_I(n).table.t(find(Table_I(n).table.Iexp >= 0, 1)-1);
    end
    
    %% Fitt the kL
    % Define objective function
    fun = @(x)sum((Table_I(n).table.Iexp(search) - (Table_I(n).I0 - Table_I(n).If) ...
        .*(1-exp(-x*Table_I(n).A/Table_I(n).V.*(Table_I(n).table.t(search) - t0 )))).^2);
    
    % Define optimization options
    % With (true) or without (false) optimization results output    
    if doOptResults == true
        options = optimset('PlotFcns',@optimplotfval,'TolX', 1e-8);
        [kL,fval,exitflag,output] = fminbnd(fun,kLbounds(1), kLbounds(2),options);
    else
        options = optimset('TolX', 1e-10);
        [kL] = fminbnd(fun,kLbounds(1), kLbounds(2),options);
    end
    % Store kL value in struct
    Table_I(n).kLfit = kL;    
    % Compute and store the calculated luminescence (photon counts)
    Table_I(n).table.Icalc = (Table_I(n).I0 - Table_I(n).If) ...
        .*(1-exp(-kL*Table_I(n).A/Table_I(n).V.*(Table_I(n).table.t - t0)));
    
    %% Plot the results
    if doPlots == true
        figure;
        hold on;
        plot(Table_I(n).table.t, Table_I(n).table.Icalc, '--', 'Color', [.620 .620 .620], 'linewidth', 3, ...
            'DisplayName', ['Fitted curve - kL = ' num2str(Table_I(n).kLfit, 3)]);
        
%         plot(Table_I(n).table.t, Table_I(n).table.Iexp, 'DisplayName', [char(Table_I(n).date) ' - ' ...
%             num2str(Table_I(n).rpm) ' rpm - ' num2str(Table_I(n).V*1e+6) ' mL']);
        
        plot(Table_I(n).table.t, Table_I(n).table.Iexp, 'DisplayName', 'Experimental');
        
        legend('boxoff');        
        legend('Position', [0.65 0.45, 0, 0]);
        legend('Orientation','vertical');
        title(legend, 'Mass transfer coefficient (m/s)')
        ylim([0 inf])
        ax = gca;
        ax.FontSize = 13;
        title([num2str(Table_I(n).rpm) ' rpm - ' num2str(Table_I(n).V*1e+6) ' mL'])
        ylabel('Luminescence (photon counts)','interpreter','latex')
        xlabel('Time (s)','interpreter','latex')
    end
end

%% Constructs table from struct

Table_ISh = struct2table(Table_I);
% Remove unwanted variables
Table_ISh = removevars(Table_ISh,{'date', 'I0', 'If', 'range1', 'range2', ...
    'ignore', 'table'});
% Rename kLfit to kL (because there will be a new kLfit) 
Table_ISh.Properties.VariableNames{'kLfit'} = 'kL';
% Remove outliers before optimization
Table_ISh([15,16],:) = [];

%% Averages repeated trials' results

% Find unique rows
[uniqueValues, uniqueRows] = unique([Table_ISh.V, ...
    Table_ISh.rpm], 'rows');

% Find and average duplicate sets (stored in the set's first ocurrence)
for n = uniqueRows'
    setRows = find(Table_ISh.rpm == Table_ISh.rpm(n) &  ...
        Table_ISh.V == Table_ISh.V(n));
    Table_ISh{n,:} = mean(Table_ISh{setRows,:}, 1);
end

% Find and delete duplicate rows
duplicateRows = setdiff(1:height(Table_ISh), uniqueRows);
Table_ISh(duplicateRows,:) = [];

%% Computes adimensional numbers

% Reynolds, Schmidt, Geometrical term, Sherwood
[Table_ISh.Re, Table_ISh.Sc, Table_ISh.G, Table_ISh.Sh] = ...
    calcAdimensionalNumbers(Table_ISh.T, P, Table_ISh.rpm, ...
    Table_ISh.A, Table_ISh.V, Table_ISh.kL);

%% Fits Sherwood constants

% With (true) or without (false) optimization results output
if doOptResults2 == true
    options = optimset('PlotFcns',@optimplotfval,'TolX', 1e-8, ...
        'MaxIter', 10000, 'MaxFunEvals', 100000);
else
    options = optimset('TolX', 1e-8, 'MaxIter', 10000, ...
        'MaxFunEvals', 100000);
end

% Define objective function and fit constants
if do4const2 == true
    fun = @(k) sum((Table_ISh.Sh - (k(1).*Table_ISh.Re.^k(2).*Table_ISh.Sc.^k(3).*Table_ISh.G.^k(4))).^2);
    k0 = [kinit(1), kinit(2), kinit(3), kinit(4)];
    kExp = fminsearch(fun, k0, options)
    Table_ISh.Shfit = (kExp(1).*Table_ISh.Re.^kExp(2).*Table_ISh.Sc.^kExp(3).*Table_ISh.G.^kExp(4));
else
    fun = @(k) sum((Table_ISh.Sh - (k(1).*Table_ISh.Re.^k(2).*Table_ISh.Sc.^kinit(3).*Table_ISh.G.^k(3))).^2);
    k0 = [kinit(1), kinit(2), kinit(4)];
    kExp = fminsearch(fun, k0, options)
    Table_ISh.Shfit = (kExp(1).*Table_ISh.Re.^kExp(2).*Table_ISh.Sc.^kinit(3).*Table_ISh.G.^kExp(3));
end

%% Computes and store Dab and kLfit
% Oxygen-water diffusion coefficient (m^2/s)
Table_ISh.Dab = calcDab(Table_ISh.T);
% Compute fitted mass transfer coefficients (m7s)
Table_ISh.kLfit = Table_ISh.Shfit.*Table_ISh.Dab./const.D;

%% Computes exponential regressions

% Find unique rows
[uniqueValues, uniqueRows] = unique([Table_ISh.V], 'rows');

PolyFitkL = [];
% Find volume sets and compute exponential regressions
for n = uniqueRows'
    setRows = find(abs(Table_ISh.V - Table_ISh.V(n)) <= setRange);
    if length(setRows) > 1
        PolyFitkL = [PolyFitkL; polyfit(log(Table_ISh.rpm(setRows,:)), ...
            log(Table_ISh.kL(setRows,:)), 1)];
    end
end

%% Sorts the table based on Re to enable plots

% Sort Re and keep the sort index in "sortIdx"
[Table_ISh.Re, sortIdx] = sort(Table_ISh.Re, "ascend");

% Sort table using the sorting index except for Re (Re is already sorted)
for n = 1:width(Table_ISh)
    if n ~= find(string(Table_ISh.Properties.VariableNames) == 'Re')
        Table_ISh(:,n) = Table_ISh(sortIdx,n);
    end
end

%% Plot fitted kL, correlated kL and power regressions

% Initialize figure
figure;
hold on;

% Plot regressed kL values
scatter(Table_ISh.rpm, Table_ISh.kL, 75, 'o', 'MarkerEdgeColor', [0.420    0.420    0.420], ...
    'DisplayName', 'Experimental')

% Define range and plot correlated values and/or exponential trendline for each volume set
rpm = 0:1:120;
% for n = 1:length(PolyFitkL(:,1))-1
%     plot(rpm, exp(PolyFitkL(n,2)).*rpm.^PolyFitkL(n,1), '--', 'Color', 'blue', ...
%         'DisplayName', '', 'HandleVisibility', 'off')
% end

for n = 1:length(folder)
    [Re, Sc, G, Sh] = ...
        calcAdimensionalNumbers([zeros(1,length(rpm))+conditionsExp(n,1)]', conditionsExp(n,2), [rpm]', ...
        conditionsExp(n,3), [zeros(1,length(rpm))+conditionsExp(n,4)]', [zeros(1,length(rpm))]');
    if do4const == true
        plot(rpm, (kExp(1).*Re.^kExp(2).*Sc.^kExp(3).*G.^kExp(4)).*calcDab([zeros(1,length(rpm))+conditionsExp(n,1)]')./const.D, '--', 'linewidth', 2, 'DisplayName', 'Sim. Corr.')
        hold on
    else
        if n == 1
            plot(rpm, (kExp(1).*Re.^kExp(2).*Sc.^kinit(3).*G.^kExp(3)).*calcDab([zeros(1,length(rpm))+conditionsExp(n,1)]')./const.D, '--', 'linewidth', 1, 'Color', [0    0.4471    0.7412], 'DisplayName', 'Exp. Corr.')
            hold on
        else
            plot(rpm, (kExp(1).*Re.^kExp(2).*Sc.^kinit(3).*G.^kExp(3)).*calcDab([zeros(1,length(rpm))+conditionsExp(n,1)]')./const.D, '--', 'linewidth', 1, 'Color', [0    0.4471    0.7412], 'HandleVisibility', 'off')
        end
    end
end

% Plot correlated kL values
scatter(Table_ISh.rpm, Table_ISh.kLfit, 150, 'x', 'MarkerEdgeColor', 'red', ...
    'DisplayName', 'Fitted')

%% Figure settings

ax = gca;
ax.FontSize = 13;
legend('boxoff')
legend('Location','northwest')
legend('Orientation','horizontal')
title('Experimental')
ylabel('k$_L$ (m/s)','interpreter','latex')
xlabel('Agitation Rate (rpm)','interpreter','latex')
annotation('textbox', [0.185, 0.54, 0.1, 0.1], 'String', ...
        ['Sh = ' num2str(kExp(1), 3) 'Re$^{' num2str(kExp(2), 3) '}$Sc$^{' num2str(kinit(3), 3) '}$G$^{' num2str(kExp(3), 3) '}$'], ...
        'interpreter', 'latex', 'Edgecolor','none', 'FontSize', 13, ...
        'Position', [0.2, 0.8, 0.5, 0]);

%% Plots correlated kL vs fitted kL

% Initialize figure
figure;
hold on
scatter(Table_ISh.kL, Table_ISh.kLfit, 75)
% Compute R^2
sqrR = 1-(1-rsquare(Table_ISh.kL, Table_ISh.kLfit)).*(length(Table_ISh.kL)-1)./(length(Table_ISh.kL)-2);
% Plot graph diagonal
plot([100:500]*10^-7, [100:500]*10^-7, '--', 'color', 'black')

%% Figure settings

ax = gca;
ax.FontSize = 13;
ylabel('Correlated k$_L$','interpreter','latex')
xlabel('Experimental k$_L$','interpreter','latex')
title('Experimental')
annotation('textbox', [0.185, 0.54, 0.1, 0.1], 'String', ...
    ['R${^2}$ = ' num2str(sqrR, 3)], 'interpreter', 'latex', ...
    'Edgecolor','none', 'FontSize', 13, 'Position', [0.3, 0.8, 0.5, 0]);
xlim([1e-5 4.5e-5])
ylim([1e-5 4.5e-5])

%% PART 3 - COMPARE SIMULATION VS EXPERIMENTAL

W = [];

% doExpBOSim = true;
for doExpBOSim = 0:1

%% Fits Sim/Exp ratio

if doExpBOSim == true
    if do4const == true
        fun = @(W) sum((Table_ISh.kL - W.*(k(1).*Table_ISh.Re.^k(2).*Table_ISh.Sc.^k(3).*Table_ISh.G.^k(4)).*Table_ISh.Dab./const.D).^2);
        w = fminsearch(fun, 0.5, options)
        Table_ISh.kLSimFit = w.*(k(1).*Table_ISh.Re.^k(2).*Table_ISh.Sc.^k(3).*Table_ISh.G.^k(4)).*Table_ISh.Dab./const.D;
    else
        fun = @(W) sum((Table_ISh.kL - W.*(k(1).*Table_ISh.Re.^k(2).*Table_ISh.Sc.^kinit(3).*Table_ISh.G.^k(3)).*Table_ISh.Dab./const.D).^2);
        w = fminsearch(fun, 0.5, options)
        Table_ISh.kLSimFit = w.*(k(1).*Table_ISh.Re.^k(2).*Table_ISh.Sc.^kinit(3).*Table_ISh.G.^k(3)).*Table_ISh.Dab./const.D;
    end
else
    if do4const == true
        fun = @(W) sum((Table_kcInstFSh.kL - W.*(kExp(1).*Table_kcInstFSh.Re.^kExp(2).*Table_kcInstFSh.Sc.^kExp(3).*Table_kcInstFSh.G.^kExp(4)).*Table_kcInstFSh.Dab./const.D).^2);
        w = fminsearch(fun, 2, options)
        Table_kcInstFSh.kLSimFit = w.*(kExp(1).*Table_kcInstFSh.Re.^kExp(2).*Table_kcInstFSh.Sc.^kExp(3).*Table_kcInstFSh.G.^kExp(4)).*Table_kcInstFSh.Dab./const.D;
    else
        fun = @(W) sum((Table_kcInstFSh.kL - W.*(kExp(1).*Table_kcInstFSh.Re.^kExp(2).*Table_kcInstFSh.Sc.^kinit(3).*Table_kcInstFSh.G.^kExp(3)).*Table_kcInstFSh.Dab./const.D).^2);
        w = fminsearch(fun, 2, options)
        Table_kcInstFSh.kLSimFit = w.*(kExp(1).*Table_kcInstFSh.Re.^kExp(2).*Table_kcInstFSh.Sc.^kinit(3).*Table_kcInstFSh.G.^kExp(3)).*Table_kcInstFSh.Dab./const.D;
    end
end

W = [W w];

%% Plots Experimental Correlated data and Simulation Correlated data

% Initialize figure
figure;
hold on;

if doExpBOSim == true
    % Plot Experimental kL (regressed)
    scatter(Table_ISh.rpm, Table_ISh.kL, 75, 'o', 'MarkerEdgeColor', [0.420    0.420    0.420], ...
        'DisplayName', 'Experimental')
    % Plot Experimental kL (based on simulation correlation)
    scatter(Table_ISh.rpm, Table_ISh.kLSimFit, 150, 'x', 'MarkerEdgeColor', 'red', ...
        'DisplayName', 'Adj. Sim.')
else
    % Plot Experimental kL (regressed)
    scatter(Table_kcInstFSh.rpm, Table_kcInstFSh.kL, 75, 'o', 'MarkerEdgeColor', [0.420    0.420    0.420], ...
        'DisplayName', 'Simulations')
    % Plot Experimental kL (based on simulation correlation)
    scatter(Table_kcInstFSh.rpm, Table_kcInstFSh.kLSimFit, 150, 'x', 'MarkerEdgeColor', 'red', ...
        'DisplayName', 'SimBOExp')
end

% Compute fitted kL for desired range and conditions
rpm = 0:1:120;
for n = 1:length(folder)
    if doExpBOSim == false
        [Re, Sc, G, Sh] = ...
        calcAdimensionalNumbers([zeros(1,length(rpm))+conditions(n,1)]', conditions(n,2), [rpm]', ...
        conditions(n,3), [zeros(1,length(rpm))+conditions(n,4)]', [zeros(1,length(rpm))]');
        if do4const == true
            plot(rpm, (k(1).*Re.^k(2).*Sc.^k(3).*G.^k(4)).*calcDab([zeros(1,length(rpm))+conditions(n,1)]')./const.D, '--', 'linewidth', 2, 'DisplayName', 'SimCorr')
            hold on
        else
            if n == 1
                plot(rpm, (k(1).*Re.^k(2).*Sc.^kinit(3).*G.^k(3)).*calcDab([zeros(1,length(rpm))+conditions(n,1)]')./const.D, '--', 'linewidth', 1, 'Color', [0    0.4471    0.7412], 'DisplayName', 'SimCorr')
                hold on
            else
                plot(rpm, (k(1).*Re.^k(2).*Sc.^kinit(3).*G.^k(3)).*calcDab([zeros(1,length(rpm))+conditions(n,1)]')./const.D, '--', 'linewidth', 1, 'Color', [0    0.4471    0.7412], 'HandleVisibility', 'off')
            end
        end 
    else
        [Re, Sc, G, Sh] = ...
        calcAdimensionalNumbers([zeros(1,length(rpm))+conditionsExp(n,1)]', conditionsExp(n,2), [rpm]', ...
        conditionsExp(n,3), [zeros(1,length(rpm))+conditionsExp(n,4)]', [zeros(1,length(rpm))]');
        if do4const == true
            plot(rpm, (kExp(1).*Re.^kExp(2).*Sc.^kExp(3).*G.^kExp(4)).*calcDab([zeros(1,length(rpm))+conditionsExp(n,1)]')./const.D, '--', 'linewidth', 2, 'DisplayName', 'SimCorr')
            hold on
        else
            if n == 1
                plot(rpm, (kExp(1).*Re.^kExp(2).*Sc.^kinit(3).*G.^kExp(3)).*calcDab([zeros(1,length(rpm))+conditionsExp(n,1)]')./const.D, '--', 'linewidth', 1, 'Color', [0    0.4471    0.7412], 'DisplayName', 'Exp. Corr.')
                hold on
            else
                plot(rpm, (kExp(1).*Re.^kExp(2).*Sc.^kinit(3).*G.^kExp(3)).*calcDab([zeros(1,length(rpm))+conditionsExp(n,1)]')./const.D, '--', 'linewidth', 1, 'Color', [0    0.4471    0.7412], 'HandleVisibility', 'off')
            end
        end
    end
end

%% Figure Settings

ax = gca;
ax.FontSize = 13;
legend('boxoff')
legend('Location','northwest')
legend('Orientation','horizontal')
title('Adj. Sim. Correlation')
ylabel('k$_L$ (m/s)','interpreter','latex')
xlabel('Agitation Rate (rpm)','interpreter','latex')

if doExpBOSim == false
    annotation('textbox', [0.185, 0.54, 0.1, 0.1], 'String', ...
        ['Sh$^{Adj. Exp}$ = ' num2str(w, 3) '(' num2str(kExp(1), 3) 'Re$^{' num2str(kExp(2), 3) '}$Sc$^{' num2str(kinit(3), 3) '}$G$^{' num2str(kExp(3), 3) '}$)'], ...
        'interpreter', 'latex', 'Edgecolor','none', 'FontSize', 13, ...
        'Position', [0.2, 0.8, 0.9, 0]);
    annotation('textbox', [0.185, 0.54, 0.1, 0.1], 'String', ...
        ['Sh$^{Sim}$ = ' num2str(k(1), 3) 'Re$^{' num2str(k(2), 3) '}$Sc$^{' num2str(kinit(3), 3) '}$G$^{' num2str(k(3), 3) '}$'], ...
        'interpreter', 'latex', 'Edgecolor','none', 'FontSize', 13, ...
        'Position', [0.2, 0.72, 0.5, 0]);
else
    annotation('textbox', [0.185, 0.54, 0.1, 0.1], 'String', ...
        ['Sh$^{Adj. Sim}$ = ' num2str(w, 3) '(' num2str(k(1), 3) 'Re$^{' num2str(k(2), 3) '}$Sc$^{' num2str(kinit(3), 3) '}$G$^{' num2str(k(3), 3) '}$)'], ...
        'interpreter', 'latex', 'Edgecolor','none', 'FontSize', 13, ...
        'Position', [0.2, 0.8, 0.9, 0]);
    annotation('textbox', [0.185, 0.54, 0.1, 0.1], 'String', ...
        ['Sh$^{Exp}$ = ' num2str(kExp(1), 3) 'Re$^{' num2str(kExp(2), 3) '}$Sc$^{' num2str(kinit(3), 3) '}$G$^{' num2str(kExp(3), 3) '}$'], ...
        'interpreter', 'latex', 'Edgecolor','none', 'FontSize', 13, ...
        'Position', [0.2, 0.72, 0.5, 0]);
end
xlim([0 120])
if doExpBOSim == true
    ylim([0 1e-4])
else
    ylim([0 2.5e-4])
end
end

%% PART 4 - ERRORS AND OTHER VALUES

%% Computes Experimental Errors

% ExpDiffList = Table_ISh.kLfit - Table_ISh.kL;
% meanDiffErr = mean(abs(ExpDiffList));
% maxDiffErr = max(abs(ExpDiffList));

ExpErrList = (Table_ISh.kLfit - Table_ISh.kL)./Table_ISh.kL*100;
meanExpErr = mean(abs(ExpErrList));
maxExpErr = max(abs(ExpErrList));

%% Computes Simulations Differences and Errors

SimDiffList = Table_kcInstFSh.kLfit - Table_kcInstFSh.kL;
meanSimDiff = mean(abs(SimDiffList));
maxSimDiff = max(abs(SimDiffList));

% SimErrList = (Table_kcInstFSh.kLfit - Table_kcInstFSh.kL)./Table_kcInstFSh.kL*100;
% meanSimErr = mean(abs(SimErrList));
% maxSimErr = max(abs(SimErrList));

AlphaErr = (k(2) - kExp(2))/kExp(2)*100;

%% Computes ExpBOSim or SimBOExp Errors

% doExpBOSim = true;
for doExpBOSim = 0:1
    if doExpBOSim == true
        ExpBOSimkL = W(2).*k(1).*Table_ISh.Re.^k(2).*Table_ISh.Sc.^kinit(3).*Table_ISh.G.^k(3).*calcDab(Table_ISh.T)./const.D;
        ExpBOSimErrList = (Table_ISh.kL - ExpBOSimkL)./ExpBOSimkL*100;
        meanExpBOSimErr = mean(abs(ExpBOSimErrList));
        maxExpBOSimErr = max(abs(ExpBOSimErrList));
    else
        SimBOExpkL = W(1).*kExp(1).*Table_kcInstFSh.Re.^kExp(2).*Table_kcInstFSh.Sc.^kinit(3).*Table_kcInstFSh.G.^kExp(3).*calcDab(Table_kcInstFSh.T)./const.D;
        SimBOExpErrList = (Table_kcInstFSh.kL - SimBOExpkL)./SimBOExpkL*100;
        meanSimBOExpErr = mean(abs(SimBOExpErrList));
        maxSimBOExpErr = max(abs(SimBOExpErrList));
    end
end

%% PART 5 - PRODUCE HISTOGRAMS

%% Initial inputs

% Select which graphs to plot [bar, area, plot, smoothPlot]
graphType = [0, 0, 1, 0];

%% Retrieves histogram values from files 

%Select histogram data: 'Kolmogorov' = 1 or 'EDR' = 2
file = {'Hist_kmgS_100mL', 'Hist_TKEDR_100mL'};

DataAnalysis.rpm = [unique(Table_kcInstFSh.rpm)];

%% Plots Kolmogorov Scale and EDR histograms and calculates ALP values

for n = 1:2
    
    %% Retrieve and compute vars
    
    opts = detectImportOptions(char(file(n)));
    histg = readtable(char(file(n)), opts);
    
    %Axis range
    if n == 1
        histg.(1) = histg.(1)*10^6;
        xrange = [0 2000];
        yrange = [0 inf];
        xtick = 200;
        ytick = 1;
    else
        xrange = [0 4e-3];
        yrange = [0 5];
        xtick = 5e-4;
        ytick = 1;
    end
    
    for m = 2:width(histg)
        if n == 1
            DataAnalysis.kmgS_ALP(m-1) = 1 - sum(histg.(m)(histg.(1) <= 350));
            DataAnalysis.meankmgS(m-1) = sum(histg.(1).*histg.(m));
            DataAnalysis.peakskmgS(m-1,:) = [histg.(1)(smooth(histg.(m)) == max(smooth(histg.(m)))) max(smooth(histg.(m)))];
        else
            DataAnalysis.TKEDR_ALP(m-1) = 1 - sum(histg.(m)(histg.(1) <= 1.5e-3));
            DataAnalysis.meanTKEDR(m-1) = sum(histg.(1).*histg.(m));
            DataAnalysis.calckmgS(m-1) = ((calcNuH2O(21, 1))^3/sum(histg.(1).*histg.(m)))^(1/4)*10^6;
        end
    end
    
    %Label tick range and spacing
    xtl=(xrange(1):xtick:xrange(2));
    ytl=(yrange(1):ytick:100);
    
    map = brewermap(7,'Set1');
    
    %% Bar graph with shade
    
    if graphType(1) == true
        figure;
        hold on
        for m = 2:width(histg)
            bar(histg.(1),histg.(m)*100,'facecolor',map(m-1,:),'facealpha',.7-(m-2)/10,'edgecolor','none', 'DisplayName', num2str(DataAnalysis.rpm(m-1), 3))
        end
        
        xticks(xtl)
        xtickangle(45)
        yticks(ytl)
        
        xlabel('Length ($ \mu m $)','interpreter','latex')
        ylabel('Distribution (\%)','interpreter','latex')
        
        box off
        axis tight
        legend('location','northeast')
        legend boxoff
        
        xlim(xrange)
        ylim(yrange)
    end
    
    %% Area graph with shade
    
    if graphType(2) == true
        figure;
        hold on
        for m = 2:width(histg)
            area(histg.(1),histg.(m)*100,'facecolor',map(m-1,:),'facealpha',.7-(m-2)/10,'edgecolor','none', 'DisplayName', num2str(DataAnalysis.rpm(m-1), 3))
        end
        
        xticks(xtl)
        xtickangle(45)
        yticks(ytl)
        
        xlabel('Length ($ \mu m $)','interpreter','latex')
        ylabel('Distribution (\%)','interpreter','latex')
        
        box off
        axis tight
        legend('location','northeast')
        legend boxoff
        
        xlim(xrange)
        ylim(yrange)
    end
    
    %% Line graph
    
    if graphType(3) == true
        figure;
        hold on
        for m = 2:width(histg)
            plot(histg.(1),histg.(m)*100, 'DisplayName', num2str(DataAnalysis.rpm(m-1), 3))
        end
        
        xticks(xtl)
        if n == 1
            xtickangle(45)
        end    
        yticks(ytl)
        
        xlabel('Length ($ \mu m $)','interpreter','latex')
        ylabel('Distribution (\%)','interpreter','latex')
        
        box off
        axis tight
%         legend('10 rpm','20 rpm','30 rpm','40 rpm','60 rpm','80 rpm','100 rpm','location','northeast')
        legend boxoff
        
        xlim(xrange)
        ylim(yrange)
        
        ax = gca;
        ax.FontSize = 13;
        
        if n == 1
            title('Kolmogorov Length Scale Distribution')
        elseif n ==2
            title('Energy Dissipation Rate Distribution')    
        end
    end
    
    %% Smooth line graph
    
    if graphType(4) == true
        figure;
        hold on
        for m = 2:width(histg)
            plot(histg.(1),smooth(histg.(m)*100), 'DisplayName', num2str(DataAnalysis.rpm(m-1), 3))
        end
        
        %set(gca, 'YScale', 'log')
        
        xticks(xtl)
        xtickangle(45)
        yticks(ytl)
        
        xlabel('Length ($ \mu m $)','interpreter','latex')
        ylabel('Distribution (\%)','interpreter','latex')
        
        box off
        axis tight
        legend('location','northeast')
        legend boxoff
        
        xlim(xrange)
        ylim(yrange)
    end
    
end

%% Computes total ALP %  

DataAnalysis.population_ALP = DataAnalysis.kmgS_ALP + DataAnalysis.TKEDR_ALP;

%% Plots ALP %

figure;
plot(DataAnalysis.rpm, DataAnalysis.kmgS_ALP)
hold on
plot(DataAnalysis.rpm, DataAnalysis.TKEDR_ALP)
plot(DataAnalysis.rpm, DataAnalysis.population_ALP)

%% Figure settings

xlabel('Agitation Rate (rpm)','interpreter','latex')
ylabel('Distribution (\%)','interpreter','latex')

ax = gca;
ax.FontSize = 13;
box off
axis tight
title('Above Limit Precentages')
legend('Kolmogorov','EDR','location','northeast')
legend boxoff

%% PART 6 - PRODUCE TEMPERATURE DEPENDENCY GRAPH

%% Initial inputs

t0 = 8;

%% Retrieve plot data

load('Dab.mat')
load('Dab2.mat')
load('Nut.mat')

%% Sort list data

[Dab(:,1), Idx] = sort(Dab(:,1));
Dab(:,2) = Dab(Idx,2);
[Dab2(:,1), Idx] = sort(Dab2(:,1));
Dab2(:,2) = Dab2(Idx,2);
Dab = [Dab; Dab2];
[Nut(:,1), Idx] = sort(Nut(:,1));
Nut(:,2) = Nut(Idx,2);

%% Plots Temperature Dependency

figure;
hold on;
plot(Dab(:,1), Dab(:,2), 'DisplayName', 'Diffusion Coefficient');
plot(Nut(:,1), Nut(:,2), 'DisplayName', 'Kinematic Viscosity');

% yline(0, 'HandleVisibility', 'off');

yline(1-1/3, '--', ['$1 - \beta$ = ' num2str(1 - 1/3, 3)], ...
    'interpreter', 'latex', 'LabelVerticalAlignment', 'top', ...
    'FontSize', 14, 'DisplayName', 'Average', 'HandleVisibility', 'on');
yline(1/3 - 0.551, '--', ['$\beta - \alpha$ = ' num2str(1/3 - 0.551, 3)], ...
    'interpreter', 'latex', 'LabelVerticalAlignment', 'bottom', ...
    'FontSize', 14, 'HandleVisibility', 'off');

% [c idx] = min(abs(Dab(:,1) - t0));
% yline(mean(Dab([idx:end],2)), '--', ['$1 - \beta$ = ' num2str(mean(Dab([idx:end],2)), 3)], ...
%     'interpreter', 'latex', 'LabelVerticalAlignment', 'top', ...
%     'FontSize', 14, 'DisplayName', 'Average', 'HandleVisibility', 'on');
% 
% [c idx] = min(abs(Nut(:,1) - t0));
% yline(mean(Nut([idx:end],2)), '--', ['$\beta - \alpha$ = ' num2str(mean(Nut([idx:end],2)), 3)], ...
%     'interpreter', 'latex', 'LabelVerticalAlignment', 'bottom', ...
%     'FontSize', 14, 'HandleVisibility', 'off');

%% Figure settings

set(gcf,'Position',[10 350 1900 500]);
ax = gca;
ax.FontSize = 13;
title('Temperature Variables'' Dependency')
legend('boxoff')
legend('Location','Northeast')
legend('Orientation','Horizontal')
ylabel('Exponent','interpreter','latex')
xlabel('Time (s)','interpreter','latex')
axis tight
ax.XGrid = 'off';
ax.YGrid = 'off';
ax.GridColor = [0.1490    0.1490    0.1490];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.2;
ax.Layer = 'top';
xlim([0 inf])
ylim([-0.4 1])

%% PART 7 - Other functions
% 
plotDab
plotCO2
plotMRS_kcInstF
plotPSlip125