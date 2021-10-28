function outputArg1 = Compare_Bending_curves(y, ...
    Bending_S, Bending_A, Bending_C, Bending_D, Bending_F, Bending_G, Bending_E, ...
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
plot(y, Bending_S, 'LineWidth', 1.5)
plot(y, Bending_A, 'LineWidth', 1.5)
plot(y, Bending_C, 'LineWidth', 1.5)
plot(y, Bending_D, 'LineWidth', 1.5)
plot(y, Bending_F, 'LineWidth', 1.5)
plot(y, Bending_G, 'LineWidth', 1.5)
plot(y, Bending_E, 'LineWidth', 1.5)

xlabel("Wing semispan - $y$ $(m)$", "Interpreter", "latex")
ylabel("Bending load $(daN \cdot m)$", "Interpreter", "latex")
title("Bending loads comparison", "Interpreter", "latex")
legend({PointS,PointA,PointC,PointD,PointF,PointG,PointE}, 'Interpreter', 'latex', 'Location', 'southeast')

end
