function outputArg1 = Pitching_moment_coefficients_diagram(y, cm_Symm, cm_Unsymm, Point)
% Unsymm_load_diagram(y, Torsion_load, Reg) 
% Torsion load distribution along the span of the main wing diagram. 
%
%  INPUT 
%  y            --> Distribution of station along the spanwise direction.
%  Torsion_load --> Torsion load distribution along the spanwise direction.
%  Point        --> Point of the flight envelope.
%
%  OUTPUT
%  Torsion load spanwise distribution diagram.

outputArg1 = figure;
hold on
grid on 
grid minor
plot(y, cm_Unsymm, '-r', 'LineWidth', 1.5)
plot(y, cm_Symm, '-.k', 'LineWidth', 1.0)
xlabel("Wing semispan - $y$ $(m)$", "Interpreter", "latex")
ylabel("Pitching mom. coeff. - $C_m$", "Interpreter", "latex")
title("Pitching moment coefficient comparison at point", Point, "Interpreter", "latex")
legend({'$C_m$ Unsymmetrical','$C_m$ Symmetrical'}, 'Interpreter', 'latex', 'Location', 'northwest')

end

