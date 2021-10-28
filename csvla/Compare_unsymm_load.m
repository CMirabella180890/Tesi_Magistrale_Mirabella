function outputArg1 = Compare_unsymm_load(y, ...
    Torsion_load_unsymm_PointA, ...
    Torsion_load_unsymm_PointC, ...
    Torsion_load_unsymm_PointD, ... 
    PointA, PointC, PointD, Load)
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
plot(y, Torsion_load_unsymm_PointA, '-g', 'LineWidth', 1.5)
plot(y, Torsion_load_unsymm_PointC, '-b', 'LineWidth', 1.5)
plot(y, Torsion_load_unsymm_PointD, '-r', 'LineWidth', 1.5)
xlabel("Wing semispan - $y$ $(m)$", "Interpreter", "latex")
ylabel("Unsymmetrical Torsion load $(daN \cdot m)$", "Interpreter", "latex")
title("Unsymmetrical Torsion loads comparison - ", Load, "Interpreter", "latex")
legend({PointA,PointC,PointD}, 'Interpreter', 'latex', 'Location', 'southeast')

end