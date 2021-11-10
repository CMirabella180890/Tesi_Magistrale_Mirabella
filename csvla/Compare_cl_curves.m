function outputArg1 = Compare_cl_curves(y, ...
    cl_S, cl_A, cl_C, cl_D, cl_F, cl_G, cl_E, ...
    PointS, PointA, PointC, PointD, PointF, PointG, PointE)
% Unsymm_load_diagram(y, cl_load, Reg) 
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
plot(y, cl_S, 'LineWidth', 1)
plot(y, cl_A, 'LineWidth', 1)
plot(y, cl_C, 'LineWidth', 1)
plot(y, cl_D, 'LineWidth', 1)
plot(y, cl_F, 'LineWidth', 1)
plot(y, cl_G, 'LineWidth', 1)
plot(y, cl_E, 'LineWidth', 1)
xlabel("Wing semispan - $y$ $(m)$", "Interpreter", "latex")
ylabel("$C_l = C_{l}(y)$", "Interpreter", "latex")
title("Lift coefficient spanwise ditribution", "Interpreter", "latex")
legend({PointS,PointA,PointC,PointD,PointF,PointG,PointE}, 'Interpreter', 'latex', 'Location', 'southeast')

end