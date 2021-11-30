function out = AeroPlot(alpha, CL, CD, CM)

% out = AeroPlot(alpha, CL, CD, CM)
% A convenient function that creates a separate figure with subplots
% representing aerodynamic polars from the aerod. data angle of attack
% (alpha), lift coefficient (CL), drag coefficient (CD), pitching moment
% coefficient (CM) relative to the aircraft.
% INPUT
%   alpha --> An [n X 1] vector (column vector) which contains angle of
%             attack (expressed in degrees or radians)
%   CL    --> An [n X 1] vector (column vector) which contains aircraft
%             lift coefficient evaluated at the prescribed alpha (Non dim.)
%   CD    --> An [n X 1] vector (column vector) which contains aircraft 
%             drag coefficient evaluated at the prescribed alpha (Non dim.)
%   CM    --> An [n X 1] vector (column vector) which contains aircraft 
%             pitching moment coefficient evaluated at the prescribed alpha
%             (Non dim.)
% OUTPUT 
%   A figure containing subplots of aerodynamic aircraft polars. 
%
% To plot user-defined aero data change axis limit and other parameters
% inside the function related to the subplots.

out = figure;
% -------------------------------------------------------------------------
% Lift coefficient vs Angle of attack
subplot(2,2,1)
plot(alpha, CL)
ylim([0.0 1.8])    % y-axis limit
xlim([-5 15])      % x-axis limit
grid on
grid minor
xlabel("$\alpha$", "Interpreter", "latex")
ylabel("$C_L$", "Interpreter", "latex")
title('Subplot 1: $C_L$ vs $\alpha$', "Interpreter", "latex")
% -------------------------------------------------------------------------
% Drag coefficient vs Angle of attack 
subplot(2,2,2)
plot(alpha, CD)
ylim([0.0 0.3])    % y-axis limit
xlim([-5 15])      % x-axis limit
grid on
grid minor
xlabel("$\alpha$", "Interpreter", "latex")
ylabel("$C_D$", "Interpreter", "latex")
title('Subplot 2: $C_D$ vs $\alpha$', "Interpreter", "latex"')
% -------------------------------------------------------------------------
% Lift coefficient vs Drag coefficient 
subplot(2,2,3)
plot(CL, CD)
ylim([0.0 0.3])     % y-axis limit
xlim([0.0 1.8])     % x-axis limit
grid on
grid minor
xlabel("$C_L$", "Interpreter", "latex")
ylabel("$C_D$", "Interpreter", "latex")
title('Subplot 2: $C_D$ vs $C_L$', "Interpreter", "latex"')
% -------------------------------------------------------------------------
% Pitching moment coefficient vs Angle of attack 
subplot(2,2,4)
plot(alpha, CM)
ylim([-0.3 0.0])    % y-axis limit
xlim([-5 15])       % x-axis limit
grid on
grid minor
xlabel("$\alpha$", "Interpreter", "latex")
ylabel("$C_M$", "Interpreter", "latex")
title('Subplot 3: $C_M$ vs $\alpha$', "Interpreter", "latex"')
exportgraphics(out,'Polars.pdf','ContentType','vector')
exportgraphics(out,'Polars.png','ContentType','vector')
end

