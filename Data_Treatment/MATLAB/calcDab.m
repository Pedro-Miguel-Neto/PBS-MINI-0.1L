function [Dab] = calcDab(T)
% Read text file
TableDab = readtable('2PDabPlot');
% T - Empirical temperatures (ºC)
% Dab - Empirical oxygen-water diffusion coefficient (m^2/s)
% Fit Dab data:
PolyFitDab = polyfit(TableDab.T, TableDab.Dab*10^-4, 1);
%Compute oxygen-water diffusion coefficient (m^2/s)
% Oxygen-water diffusion coefficient (m/s)
Dab = PolyFitDab(1)*T + PolyFitDab(2);