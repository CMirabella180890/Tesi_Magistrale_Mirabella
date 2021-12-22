function fig = Schrenk_plots(eta, chord, elliptical_load, schrenk_load)
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
   plot(eta, chord, "-.k", "LineWidth", 1.5);
   plot(eta, elliptical_load, "-r", "LineWidth", 1); 
   plot(eta, schrenk_load, "-b", "LineWidth", 1);
   xlabel("Non-dimensional semi-span - $\eta$", "Interpreter", "latex")
   ylabel("Span-wise load and chord - $(m)$", "Interpreter", "latex")
   title("Schrenk load distribution across the span", "Interpreter", "latex")
   legend(["Chord", "Elliptical load", "Schrenk load"], 'Interpreter', 'latex', 'Location', 'northeast')
   
end