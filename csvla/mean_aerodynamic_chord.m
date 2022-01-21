function MAC = mean_aerodynamic_chord(c_root, taper_ratio)
% MAC = mean_aerodynamic_chord(c_root, taper_ratio)
%   A mean aerodynamic chord calculator.
%   INPUT 
%   taper_ratio = Taper ratio of the wing.
%   c_root      = Root chord of the wing.
%   OUTPUT
%   MAC         = Mean aerodynamic chord.

x   = (2 / 3) * c_root;
y   = 1 + taper_ratio + taper_ratio^2;
z   = 1 + taper_ratio;
MAC = x * (y / z);
end

