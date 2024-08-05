function Phi = GHKflux_fast(Vm, Cin, Cout, P, charge, T)
%GHKFLUX
% charge ... charge [K]
% T ... temperature [K]
%
% Phi is the current density (flux) across the membrane carried by ion S, measured in amperes per square meter (A·m−2)
% arguments
%     Vm double % transmembrane potential in volts
%     Cin double % intracellular concentration of ion S, measured in mol·m−3 or mmol·l−1
%     Cout double % extracellular concentration of ion S, measured in mol·m−3 or mmol·l−1
%     P double = 1; % permeability of the membrane for ion S measured in m·s−1
% end

F = 96485.3328959;   % [C/mol]
R = 8.314459848;     % [J/mol/K]

% Phi = P * charge^2 * (Vm*F^2 / (R*T)) * (Cin - Cout*exp(-charge*Vm*F/(R*T))) / (1 - exp(-charge*Vm*F/(R*T)));

U = charge * Vm * F / R / T;
emU = exp(-U);
% epU = exp(U);

Phi = P .* charge * F * U .* (Cin - Cout.*emU) ./ (1 - emU);

% Phi_in = P .* charge * F * U .* Cin ./ (1 - emU);
% Phi_out = P .* charge * F * U .* Cout ./ (1 - epU);

end

