function [nu] = calcNuH2O(T, P)
load('constants')
% Compute water density (kg/m^3)
rhoH2O = calcRhoH2O(T, P);
% Water kinematic viscosity (m^2/s)
nu = const.mu20*10.^((20-T)./(96+T).*(const.nu(1)+const.nu(2).*(20-T) ...
    +const.nu(3).*(20-T).^2+const.nu(4).*(20-T).^3))./rhoH2O.*10^-3;