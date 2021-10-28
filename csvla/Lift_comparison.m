function outputArg1 = Lift_comparison(AoA, CL_NonLin, CL_Lin, CL_Data, alpha)
% Lift_comparison(AoA, CL_NonLin, CL_Lin, CL_Data)
%  Graphical comparison between lift models and data. 
%
%  INPUT 
%  AoA       --> Distribution of angle of attack.
%  CL_NonLin --> Non linear lift coefficient curve.
%  CL_Lin    --> Linear lift coefficient curve.
%  CL_Data   --> Data points about lift coefficient available
%
%  OUTPUT
%  Lift comparison diagram.

outputArg1 = figure;
hold on
grid on 
grid minor
plot(AoA, CL_NonLin, '-r', 'LineWidth', 1.5)
plot(AoA, CL_Lin, '-b', 'LineWidth', 1.0)
plot(alpha, CL_Data, 'k.', 'MarkerSize', 10)
% xlim([]);
ylim([-2.0 2.0]);
xlabel("Angle of attack - $\alpha$ $(deg)$", "Interpreter", "latex")
ylabel("Lift coefficient - $C_L$", "Interpreter", "latex")
title("Lift curve model comparison", "Interpreter", "latex")
legend({'Non linear','Linear','Data'}, 'Interpreter', 'latex', 'Location', 'southeast')

end

