function fig = Schrenk_Cl_diagram(eta, Cl)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   TEST CASE: 
%   S      = 66.5 squared feet |  S      = 6.178052 squared meters
%   b      = 19.0 feet         |  b      = 5.7912 meters
%   c_root = 5.0 feet          |  c_root = 1.524 meters
%   c_tip  = 2.0 feet          |  c_tip  = 0.6096 meters

   fig = figure;
   hold on 
   grid on 
   grid minor 
   xlim([0 1.0+0.05]);
   plot(eta, Cl, "-k", "LineWidth", 1.0);
   xlabel("Non-dimensional semi-span - $\eta$", "Interpreter", "latex")
   ylabel("Span-wise lift coefficient - $C_l $", "Interpreter", "latex")
   title("Schrenk lift coefficient across the span", "Interpreter", "latex")
   
end
