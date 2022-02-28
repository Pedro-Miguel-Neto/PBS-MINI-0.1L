function [rhoH2O] = calcRhoH2O(T, P)

load('constants')
% Uncorrected water density (kg/m^3)
rhoH2O = const.a(1) + const.a(2).*T + const.a(3).*T.^2 ...
    + const.a(4).*T.^3 + const.a(5).*T.^4;

kt = const.b(1) + const.b(2).*T + const.b(3).*T.^2 ...
    + const.b(4).*T.^3 + const.b(5).*T.^4;

% Water density (kg/m^3)
rhoH2O = rhoH2O.*(1+kt.*(P*100-101.325));