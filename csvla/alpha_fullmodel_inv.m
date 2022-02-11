function alfa_interp_inv = alpha_fullmodel_inv(CL, CL_max_inv, CL0, CLalfa, alfa0l)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

        alpha_lin   = @(CL, CL0, CLalfa) (CL - CL0) / CLalfa;

            if CL < CL_max_inv 
                alfa_interp_inv = alpha_lin(CL, CL0, CLalfa) - alfa0l;
            elseif CL > CL_max_inv
                CL = CL_max_inv;
                alfa_interp_inv = alpha_lin(CL, CL0, CLalfa) - alfa0l;
            end

end