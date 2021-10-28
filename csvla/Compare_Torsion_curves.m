function outputArg1 = Compare_Torsion_curves(y, ...
    Torsion_S, Torsion_A, Torsion_C, Torsion_D, Torsion_F, Torsion_G, Torsion_E, ...
    PointS, PointA, PointC, PointD, PointF, PointG, PointE)
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
plot(y, Torsion_S, 'LineWidth', 1.5)
plot(y, Torsion_A, 'LineWidth', 1.5)
plot(y, Torsion_C, 'LineWidth', 1.5)
plot(y, Torsion_D, 'LineWidth', 1.5)
plot(y, Torsion_F, 'LineWidth', 1.5)
plot(y, Torsion_G, 'LineWidth', 1.5)
plot(y, Torsion_E, 'LineWidth', 1.5)

xlabel("Wing semispan - $y$ $(m)$", "Interpreter", "latex")
ylabel("Torsion load $(daN \cdot m)$", "Interpreter", "latex")
title("Torsion loads comparison", "Interpreter", "latex")
legend({PointS,PointA,PointC,PointD,PointF,PointG,PointE}, 'Interpreter', 'latex', 'Location', 'southeast')

end
