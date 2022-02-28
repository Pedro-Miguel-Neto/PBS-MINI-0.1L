function [H] = calcHenryH2O(T, P)
load('constants')
%Compute water density (kg/m^3)
rhoH2O = calcRhoH2O(T, P);
% Henry's constant (mol/L/bar)
H = const.H25/const.Pbar.*rhoH2O./1000 ...
    .*exp(const.DsolHR.*(1./(T + 273.15) - 1/(25 + 273.15)));