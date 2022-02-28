function [pO2sat] = calcPO2sat(T, P)
% Calculate oxygen saturation concentration (mol/m^3)
Csat = calcCsat(T, P);
% Compute Henry's constant (mol/L/bar)
H = calcHenryH2O(T, P);
% Oxygen saturation partial pressure (bar)
pO2sat = Csat./H.*10^-3;