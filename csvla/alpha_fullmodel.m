function alfa_interp = alpha_fullmodel(CL, a2, b2, c2, CL_max, ...
    CL_star, alfa_aux, p_cl_lin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

        alpha_plus  = @(CL, a, b, c) (-b + sqrt(b^2 - 4*a*(c - CL)))/(2*a);

        alfa_interp = zeros(length(alfa_aux), 1);

        for i = 1:length(alfa_aux)
            if CL < CL_star 
                alfa_interp(i) = polyval(p_cl_lin, CL);
            elseif CL > CL_star
                alfa_interp(i) = alpha_plus(CL, a2, b2, c2);
                if ~isreal(alfa_interp(i)) == 1
%                     CL = CL_max + 0.033735;
                    CL = CL_max;
                    alfa_interp(i) = alpha_plus(CL, a2, b2, c2);
                elseif ~isreal(alfa_interp(i)) == 0
                    alfa_interp(i) = alpha_plus(CL, a2, b2, c2);
                end
            end
        end

end

