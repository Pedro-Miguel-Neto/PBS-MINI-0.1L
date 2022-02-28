function [CO2calc] = calcCO2(CO2sat, kL, A, V, t)
% Compute the calculated oxygen concentration (mol/m^3)
CO2calc = CO2sat.*(1-exp(-kL*A/V.*t));