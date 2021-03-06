% function outputArg1 = Lift_coefficients_Points(alpha_S, alpha_A, alpha_C, alpha_D, alpha_F, alpha_G, alpha_E, ...
%                                                CL_S,    CL_A,    CL_C,    CL_D,    CL_F,    CL_G,    CL_E, ...
%                                                Point_S, Point_A, Point_C, Point_D, Point_F, Point_G, Point_E, ...
%                                                alpha_dot, CL_dot, alpha_full, CL_full)
                                           
function outputArg1 = Lift_coefficients_Points(alpha_S, alpha_A, alpha_C, alpha_D, ...
                                               CL_S,    CL_A,    CL_C,    CL_D, ...
                                               Point_S, Point_A, Point_C, Point_D, ...
                                               alpha_dot, CL_dot, alpha_full, CL_full)
% Lift_comparison(AoA, CL_NonLin, CL_Lin, CL_Data)
%  Wing-body lift curve diagram. 
%
%  INPUT 
%  AoA          --> Distribution of angle of attack.
%  CL_fullmodel --> Complete lift coefficient curve.
%  CL_inverted  --> Complete lift coefficient curve for inverted flight.
%  CL_Data      --> Data points about lift coefficient available.
%
%  OUTPUT
%  Lift curve diagram.

outputArg1 = figure;
hold on
grid on 
grid minor
% plot(AoA, CL_fullmodel, '-r', 'LineWidth', 1.5)
% plot(AoA, CL_inverted, '-b', 'LineWidth', 1.5)
plot(alpha_S, CL_S, '.', ...
    'LineWidth', 2, ...
    'MarkerSize', 20, ...
    'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', [0.5,0.5,0.5])
plot(alpha_A, CL_A, '.', ...
    'LineWidth', 2, ...
    'MarkerSize', 20, ...
    'MarkerEdgeColor', 'g', ...
    'MarkerFaceColor', [0.5,0.5,0.5])
plot(alpha_C, CL_C, '.', ...
    'LineWidth', 2, ...
    'MarkerSize', 20, ...
    'MarkerEdgeColor', 'b', ...
    'MarkerFaceColor', [0.5,0.5,0.5])
plot(alpha_D, CL_D, '.', ...
    'LineWidth', 2, ...
    'MarkerSize', 20, ...
    'MarkerEdgeColor', 'c', ...
    'MarkerFaceColor', [0.5,0.5,0.5])
% plot(alpha_F, CL_F, '.', ...
%     'LineWidth', 2, ...
%     'MarkerSize', 20, ...
%     'MarkerEdgeColor', 'm', ...
%     'MarkerFaceColor', [0.5,0.5,0.5])
% plot(alpha_G, CL_G, '.', ...
%     'LineWidth', 2, ...
%     'MarkerSize', 20, ...
%     'MarkerEdgeColor', 'y', ...
%     'MarkerFaceColor', [0.5,0.5,0.5])
% plot(alpha_E, CL_E, '.', ...
%     'LineWidth', 2, ...
%     'MarkerSize', 20, ...
%     'MarkerEdgeColor', 'k', ...
%     'MarkerFaceColor', [0.5,0.5,0.5])
plot(alpha_dot, CL_dot, '.', ...
    'LineWidth', 2, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0.5,0.5,0.5])
plot(alpha_full, CL_full, '-k', ...
    'LineWidth', 0.25, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0.5,0.5,0.5])
% xlim([]);
% ylim([0.0 1.80]);
xlabel("Angle of attack - $\alpha$ $(deg)$", "Interpreter", "latex")
ylabel("Wing Body Lift coefficient - $(C_L)_wb$", "Interpreter", "latex")
title("Lift curve model", "Interpreter", "latex")
% legend({'Full model','Inverted','Data'}, 'Interpreter', 'latex', 'Location', 'northeast')
legend({Point_S, Point_A, Point_C, Point_D}, 'Interpreter', 'latex', 'Location', 'northwestoutside')
% legend({Point_S, Point_A, Point_C, Point_D, Point_F, Point_G, Point_E}, 'Interpreter', 'latex', 'Location', 'northwestoutside')

end