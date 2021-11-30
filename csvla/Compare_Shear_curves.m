function outputArg1 = Compare_Shear_curves(y, ...
    Shear_S, Shear_A, Shear_C, Shear_D, Shear_F, Shear_G, Shear_E, ...
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
plot(y, Shear_S, 'LineWidth', 1.5)
plot(y, Shear_A, 'LineWidth', 1.5)
plot(y, Shear_C, 'LineWidth', 1.5)
plot(y, Shear_D, 'LineWidth', 1.5)
plot(y, Shear_F, 'LineWidth', 1.5)
plot(y, Shear_G, 'LineWidth', 1.5)
plot(y, Shear_E, 'LineWidth', 1.5)

xlabel("Wing semispan - $y$ $(m)$", "Interpreter", "latex")
ylabel("Shear load $(daN)$", "Interpreter", "latex")
title("Shear loads comparison", "Interpreter", "latex")
legend({PointS,PointA,PointC,PointD,PointF,PointG,PointE}, 'Interpreter', 'latex', 'Location', 'southeast')

end
