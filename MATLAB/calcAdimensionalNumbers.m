function [Re, Sc, G, Sh] = calcAdimensionalNumbers(T, P, rpm, A, V, kL)
load('constants')
% Water kinematic viscosity (m^2/s)
nu = calcNuH2O(T, P);
% Oxygen-water diffusion coefficient (m^2/s)
Dab = calcDab(T);
% Reynolds number
Re = const.D^2*rpm/60./nu;
% Schmidt number
Sc = nu./Dab;
% Geometrical term
G = const.D.*A./V;
%Sherwood number
Sh = kL*const.D./Dab;