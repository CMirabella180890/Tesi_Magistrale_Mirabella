function out1 = Schrenk_load_distr(b, S, c_root, c_tip, n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   TEST CASE: 
%   S      = 66.5 squared feet |  S      = 6.178052 squared meters
%   b      = 19.0 feet         |  b      = 5.7912 meters
%   c_root = 5.0 feet          |  c_root = 1.524 meters
%   c_tip  = 2.0 feet          |  c_tip  = 0.6096 meters

   y   = linspace(0, b/2, n)';
   y   = flip(y);
   % Ellipse 
   eta = y/(b/2);
   % Chord distr. for tapered wing 
   taper_ratio = c_tip/c_root;
   c_y         = ((2*S)/((1 + taper_ratio)*b))*(1 - (1 - taper_ratio)*abs(eta));
   ell_height  = ((4*S)/(pi*b));
   x           = sqrt(1 - (eta.^2));
   Ellipse     = ell_height*x;
   Schrenk_cCl = (c_y + Ellipse)*0.5;
   Unit_Cl     = (Schrenk_cCl)./(c_y);
   disp(" ++++ CHECK ++++ ")
   disp(" Cl = Cl(y) ")
   fprintf(" %f", trapz(flip(eta), Unit_Cl))
   out1 = [y eta c_y Ellipse Schrenk_cCl Unit_Cl];
   
end

