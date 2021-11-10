function outputArg1 = Compare_cm_curves(y, ...
    cm_S, cm_A, cm_C, cm_D, cm_F, cm_G, cm_E, ...
    PointS, PointA, PointC, PointD, PointF, PointG, PointE)
% Unsymm_load_diagram(y, cd_load, Reg) 
% cl, cd, cm load distribution along the span of the main wing diagram. 
%    cd_S, cd_A, cd_C, cd_D, cd_F, cd_G, cd_E, ...
%    cm_S, cm_A, cm_C, cm_D, cm_F, cm_G, cm_E, ...
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
plot(y, cm_S, 'LineWidth', 1)
plot(y, cm_A, 'LineWidth', 1)
plot(y, cm_C, 'LineWidth', 1)
plot(y, cm_D, 'LineWidth', 1)
plot(y, cm_F, 'LineWidth', 1)
plot(y, cm_G, 'LineWidth', 1)
plot(y, cm_E, 'LineWidth', 1)
xlabel("Wing semispan - $y$ $(m)$", "Interpreter", "latex")
ylabel("$C_m = C_{m}(y)$", "Interpreter", "latex")
title("Pitch mom. coefficient spanwise ditribution", "Interpreter", "latex")
legend({PointS,PointA,PointC,PointD,PointF,PointG,PointE}, 'Interpreter', 'latex', 'Location', 'southeast')