function outputArg1 = Lift_fullmodel_curve(AoA, CL_fullmodel, CL_Data, alpha)
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
plot(AoA,   CL_fullmodel, '-r', 'LineWidth', 1.5)
plot(alpha, CL_Data,      'k.', 'MarkerSize', 10)
    xlim 'padded' ;
    ylim 'padded' ;
xlabel("Angle of attack - $\alpha$ $(deg)$", "Interpreter", "latex")
ylabel("Lift coefficient - $C_{L_{wb}}$", "Interpreter", "latex")
title("Lift curve model", "Interpreter", "latex")
legend({'Model','Data points'}, 'Interpreter', 'latex', 'Location', 'southeast')

end