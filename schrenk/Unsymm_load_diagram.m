function outputArg1 = Unsymm_load_diagram(y, Torsion_load_unsymm, Torsion_load, Point)
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
plot(y, Torsion_load_unsymm, '-r', 'LineWidth', 1.5)
plot(y, Torsion_load, '-.k', 'LineWidth', 1.0)
xlabel("Wing semispan - $y$ $(m)$", "Interpreter", "latex")
ylabel("Unsymmetrical Torsion load $(daN \cdot m)$", "Interpreter", "latex")
title("Unsymmetrical Torsion load due to aileron deflection at", Point, "Interpreter", "latex")
legend({'Unsymmetrical Torsion','Symmetrical Torsion'}, 'Interpreter', 'latex', 'Location', 'southeast')

end

