function OutputArg1 = cm_interpolation_graph(x, y, Cmy, Cm_interp, XI, YI)

% Output = cm_interpolation_graph(x, y, Cmy, Cm_interp, XI, YI)
%   A function which is able to plot a 3-dimensional interpolation of
%   values of pitching moment coefficient distributions along the span.
%
%   To illustrate the results of the 3-dimensional interpolation a carpet
%   plot of interpolated values is implemented in this function. 
%
%   INPUT
%   x         --> Point along x-axis of the meshgrid; it must be remembered
%                 that the total number of elements along x inside the
%                 vector is related to the number of subdivisions adopted
%                 with OPEN VSP.
%   y         --> Point along y-axis of the meshgrid; it must be remembered
%                 that y values corresponds to angle of attack at whic the
%                 coefficient is evaluated.
%   Cmy       --> Cm = Cm(y) obtained from OPEN VSP; remember that these
%                 values are referred to the wing semi-span.
%   Cm_interp --> Interpolated pitching moment coefficients, from interp2
%                 MATLAB function.
%   XI        --> Point at which the interpolations are performed along the
%                 x direction; in this particular case, for convenience, XI
%                 will have the same number of elements.
%   YI        --> Point at which the interpolations are performed along the
%                 y direction; in this case YI will have a greater number
%                 of elements to cover more angle of attack.
%
%   OUTPUT
%   Carpet plot of pitching moment coefficients.

OutputArg1 = figure;
plot3(x, y(1)*ones(1, length(x)), Cmy(1,:), '-b.', 'MarkerSize', 10, 'LineWidth', 2)
hold on

    for i = 2:length(y)
        plot3(x, y(i)*ones(1, length(x)), Cmy(i,:), '-b.', 'MarkerSize', 10, 'LineWidth', 2)
    end
    
grid on 
mesh(XI, YI, Cm_interp, 'FaceAlpha', '0.15')

% Axis and title
xlabel("X", "Interpreter", "latex")
ylabel("Y", "Interpreter", "latex")
zlabel("Cm = Cm(y)", "Interpreter", "latex")
title("Interpolation of the Cm values", "Interpreter", "latex") % Applied regulation from 'Aircraft' struct

end
