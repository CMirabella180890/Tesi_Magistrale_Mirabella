function alfa_interp = alpha_fullmodel(CL, a2, b2, c2, CL_max, ...
    CL_star, CL0, CLalfa)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

        alpha_lin   = @(CL, CL0, CLalfa) (CL - CL0) / CLalfa;
        alpha_plus  = @(CL, a, b, c) (-b + sqrt(b^2 - 4*a*(c - CL)))/(2*a);

        if CL > CL_max
            disp("ERROR: The code will stop now!")
            error("WARNING: There is no solution for alpha with this lift coefficient! ")
%            return
        end
        
            if CL < CL_star 
                alfa_interp = alpha_lin(CL, CL0, CLalfa);
            elseif CL > CL_star
                alfa_interp = alpha_plus(CL, a2, b2, c2);
                if ~isreal(alfa_interp) == 1
                    disp(" WARNING: negative number under square root ")
                    disp("          CL is greater than CLmax wingbody ")
                    return
%                     CL = CL_max;
%                     alfa_interp = alpha_plus(CL, a2, b2, c2);
                elseif ~isreal(alfa_interp) == 0
                    alfa_interp = alpha_plus(CL, a2, b2, c2);
                end
            end

end

