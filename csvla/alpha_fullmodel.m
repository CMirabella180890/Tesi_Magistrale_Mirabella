function alfa_interp = alpha_fullmodel(CL, a2, b2, c2, CL_max, ...
    CL_star, CL0, CLalfa)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

        alpha_lin   = @(CL, CL0, CLalfa) (CL - CL0) / CLalfa;
        alpha_plus  = @(CL, a, b, c) (-b + sqrt(b^2 - 4*a*(c - CL)))/(2*a);

            if CL < CL_star 
                alfa_interp = alpha_lin(CL, CL0, CLalfa);
            elseif CL > CL_star
                alfa_interp = alpha_plus(CL, a2, b2, c2);
                if ~isreal(alfa_interp) == 1
%                     CL = CL_max + 0.033735;
                    CL = CL_max;
                    alfa_interp = alpha_plus(CL, a2, b2, c2);
                elseif ~isreal(alfa_interp) == 0
                    alfa_interp = alpha_plus(CL, a2, b2, c2);
                end
            end

end

