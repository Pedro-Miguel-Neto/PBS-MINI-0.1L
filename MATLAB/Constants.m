%%% Constants
%Impeller diameter (m)
const.D = 0.0415;
% Ref Henry's constant (mol/L/atm)
const.H25 = 0.0013;
% -DeltaSolH/R (K)
const.DsolHR = 1500;
% P(atm) -> P(bar)
const.Pbar = 1.01325;
% Molecular oxygen molar mass (g/mol)
const.MMO2 = 2*15.999;
% Water density coefficients
const.a = [999.84847 6.337563e-2 -8.523829e-3 6.943248e-5 -3.821216e-7];
% Thermal expansion coefficients
const.b = [50.83101e-8 -3.68293e-9 7.263725e-11 -6.597702e-13 2.87767e-15];
% Reference water dynamic viscosity (mPa.s)
const.mu20 = 1.002;
% Water viscosity coefficients
const.nu = [1.2364 -1.37e-3 5.7e-6 2.55e-8];