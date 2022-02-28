%% Program Options
% List fraction
frac = 0.05;
% Define optimization bounds (min,max)
bounds = [1e-6 1e-4];
% With (true) or without (false) optimization results output
showOptResults = false;
% Plot (true) or not (false) fitted CO2 in separate plots
doSepPlots = false;

% Constants [T, P, A, V]
conditions = [21 1 1.635565e-3 9.86746e-5; 21 1 1.60185e-3 7.997321e-5; 37 1 1.56514e-3 5.996167e-5];

% With (true) or without (false) optimization results output
doOptResults2 = false;
% Fit 4 (true) or 2 (false) Sherwood constants
do4const = false;
% Initial estimates [k, Re, Sc, G]
kinit = [1.6497    0.5541  1/3  1.1658];
% Set range (L) (ex: 100mL +- 5mL)
setRange = 1e-5;

folder = {'Slip100', 'Slip80', 'Slip60'};

 for m = 1:length(folder)

    
    Table_CO2.(char(folder(m))) = table2struct(getFileNumber(char(folder(m))));
    
    % %% Trial Conditions
    % Temperature (ºC)
    % Pressure (atm)
    % Contact area (m^2)
    % Liquid volume (m^3)
    conditions = [21 1 0.001635587 0.0001; 21 1 0.001599264 0.00008; 37 1 0.00156294 0.00006];
    
    
    % Constants
    load('constants.mat')
    
    % Load Delimeted Text Import Options
    load('opts.mat');
    
    if doSepPlots == false
        figure;
    end
    
    % --------------------------------------------------------- Main Function
    for n = 1:length([Table_CO2.(char(folder(m))).rpm])
        
        %% Add conditions
        
        Table_CO2.(char(folder(m)))(n).T = conditions(m,1);
        Table_CO2.(char(folder(m)))(n).P = conditions(m,2);
        Table_CO2.(char(folder(m)))(n).A = conditions(m,3);
        Table_CO2.(char(folder(m)))(n).V = conditions(m,4);
        
        %% Compute
        % Oxygen saturation concentration (mol/m^3)
        Table_CO2.(char(folder(m)))(n).Csat = calcCsat(Table_CO2.(char(folder(m)))(n).T, Table_CO2.(char(folder(m)))(n).P);
        
        %% Get CO2sim and plot it
        path = dir([char(folder(m)) '/' Table_CO2.(char(folder(m)))(n).name]);
        numfolders = sum([path.isdir]);
        t = [];
        CO2sim = [];
        
        for k = (3:numfolders)
            trial = readtable([char(folder(m)) '/' Table_CO2.(char(folder(m)))(n).name '/' path(k).name '/CO2bF'], opts);
            trial([1],:) = [];
            t = [t; trial.Var1];
            CO2sim = [CO2sim; trial.Var2];
        end
        
        Table_CO2.(char(folder(m)))(n).trial = table(t, CO2sim);
        
        %% Fitt the kL
        % Define objective function (In case of error, check if you are in the
        % zPostProcessing folder)
        fun = @(kL)sum((Table_CO2.(char(folder(m)))(n).trial.CO2sim ...
            (round(length(Table_CO2.(char(folder(m)))(n).trial.CO2sim)*(1-frac)):end,:) ...
            - calcCO2(Table_CO2.(char(folder(m)))(n).Csat, kL, Table_CO2.(char(folder(m)))(n).A, Table_CO2.(char(folder(m)))(n).V, ...
            Table_CO2.(char(folder(m)))(n).trial.t(round(length(Table_CO2.(char(folder(m)))(n).trial.t) ...
            *(1-frac)):end,:))).^2);
        
        
        %% Define optimization options
        if showOptResults == true
            options = optimset('PlotFcns',@optimplotfval,'TolX', 1e-7);
            [kL,fval,exitflag,output] = fminbnd(fun,bounds(1), bounds(2), options);
        else
            options = optimset('TolX', 1e-7);
            [kL] = fminbnd(fun,bounds(1), bounds(2),options);
        end
        
        % Store kL value in kLlist
        Table_CO2.(char(folder(m)))(n).kL = kL;
        
        %% Calculate the CO2calc
        % Calculated oxygen concentration (mol/m^3)
        Table_CO2.(char(folder(m)))(n).trial.CO2calc = calcCO2(Table_CO2.(char(folder(m)))(n).Csat, ...
            Table_CO2.(char(folder(m)))(n).kL, Table_CO2.(char(folder(m)))(n).A, Table_CO2.(char(folder(m)))(n).V, Table_CO2.(char(folder(m)))(n).trial.t);
        
        %% Plot pO2exp and pO2calc vs t
        % ----------------------------------------------------------- Color map
        map = brewermap(7,'Set1');
        % -------------------------------------------------------- Figure input
        if doSepPlots == true
            figure;
        end
        
        % -------------------------------------------------------- Merge graphs
        hold on;
        
        
        plot(t, CO2sim, '--', 'Color', [.620 .620 .620], 'linewidth', 3, ...
            'DisplayName', '', 'HandleVisibility', 'off');
        % Smooth plot pO2calc
        plot(t, smooth(Table_CO2.(char(folder(m)))(n).trial.CO2calc), 'linewidth', 1, ...
            'DisplayName', [num2str(Table_CO2.(char(folder(m)))(n).rpm) ' rpm - kL = ' num2str(kL, 3)]);
        
        
        % -------------------------------------------------- Label descriptions
        % x-Axis
        xlabel('Time (s)','interpreter','latex');
        % y-Axis
        ylabel('C$_{O_2}$ (mol/m$^3$)','interpreter','latex');
        % -------------------------------------------- Figure Position and Size
        %[HorzPs, VertPs, HorzSz, VertSz]
        set(gcf,'Position',[325 150 1200 800]);
        % -------------------------------------------------- Box configurations
        box off;
        % ----------------------------------- Axis and axis grid configurations
        axis tight;
        ax = gca;
        ax.XGrid = 'off';
        ax.YGrid = 'off';
        ax.GridColor = [0.1490    0.1490    0.1490];
        ax.GridLineStyle = '--';
        ax.GridAlpha = 0.2;
        ax.Layer = 'top';
        % ----------------------------------------------- Legend configurations
        legend('boxoff');
        legend('Position', [0.35 0.75, 0, 0]);
        legend('Orientation','vertical');
        title(legend,'Mass transfer coefficient (m/s)');
        % ---------------------------------------------------------- Axis range
        % y-Axis
        xrange = [0 inf];
        xlim([0 inf]);
        % y-Axis
        yrange = [0 inf];
        ylim([0 inf]);
        % --------------------------------------- Tick range, spacing and angle
        % x-Axis
        xtickangle(0);
        ytickangle(0);
    end    
 end
 
%% Receives a struct and merges substructs to form a table
% Construct table fom struct
Table_CO2Sh = struct2table(Table_CO2.(char(folder(1))));   
for n = 2:length(folder)
    Table_CO2Sh.rpm(end+1) = Table_CO2.(char(folder(n))).rpm;
    Table_CO2Sh.T(end) = Table_CO2.(char(folder(n))).T;
    Table_CO2Sh.P(end) = Table_CO2.(char(folder(n))).P;
    Table_CO2Sh.A(end) = Table_CO2.(char(folder(n))).A;
    Table_CO2Sh.V(end) = Table_CO2.(char(folder(n))).V;
    Table_CO2Sh.kL(end) = Table_CO2.(char(folder(n))).kL;
end
% Remove unwanted fields
Table_CO2Sh = removevars(Table_CO2Sh,{'name', 'trial'});

%% Computes Reynolds, Schmidt, Geometrical term and Sherwood
[Table_CO2Sh.Re, Table_CO2Sh.Sc, Table_CO2Sh.G, Table_CO2Sh.Sh] = ...
    calcAdimensionalNumbers(Table_CO2Sh.T, Table_CO2Sh.P, Table_CO2Sh.rpm, ...
    Table_CO2Sh.A, Table_CO2Sh.V, Table_CO2Sh.kL);

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
    fun = @(k) sum((Table_CO2Sh.Sh - (k(1).*Table_CO2Sh.Re.^k(2).*Table_CO2Sh.Sc.^k(3).*Table_CO2Sh.G.^k(4))).^2);
    k0 = [kinit(1), kinit(2), kinit(3), kinit(4)];
    k = fminsearch(fun, k0, options)
    Table_CO2Sh.Shfit = (k(1).*Table_CO2Sh.Re.^k(2).*Table_CO2Sh.Sc.^k(3).*Table_CO2Sh.G.^k(4));
else
    fun = @(k) sum((Table_CO2Sh.Sh - (k(1).*Table_CO2Sh.Re.^k(2).*Table_CO2Sh.Sc.^kinit(3).*Table_CO2Sh.G.^k(3))).^2);
    k0 = [kinit(1), kinit(2), kinit(4)];
    k = fminsearch(fun, k0, options)
    Table_CO2Sh.Shfit = (k(1).*Table_CO2Sh.Re.^k(2).*Table_CO2Sh.Sc.^kinit(3).*Table_CO2Sh.G.^k(3));
end

%% Computes Dab and kLfit for the table

%Oxygen-water diffusion coefficient (m^2/s)
Table_CO2Sh.Dab = calcDab(Table_CO2Sh.T);
% Compute fitted mass transfer coefficients (m/s)
Table_CO2Sh.kLfit = Table_CO2Sh.Shfit.*Table_CO2Sh.Dab./const.D;
