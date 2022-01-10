function CL_interp = CL_fullmodel(alfa)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

        a           = -0.021837866;
        b           = 0.436064773;
        c           = -0.56312855; 
        alpha_plus  = @(CL) (-b + sqrt(b^2 - 4*a*(c - CL)))/(2*a);

        numb        = 1e3;
        
        alfa_data   = [-4; 0; 4; 8; 10; 12; 13];
        CL_data     = [0.43483; 0.84054; 1.21541; 1.52777; 1.61373; 1.52500; 1.47190];

        alfa_dot1   = [-4;       8];
        CL_dot1     = [0.43483;  1.52777];

        p1          = polyfit(alfa_dot1, CL_dot1, 1);

%         alfa_i      = - 4.0; 
%         alfa_f      = 13.0; 
%         alfa_interp = linspace(alfa_i, alfa_f, numb)';

        CL_star     = 1.47;
        alfa_star   = alpha_plus(CL_star);

        alfa_dot2   = [alfa_star; 8;       10;      12;      13];
        CL_dot2     = [CL_star;   1.52777; 1.61373; 1.52500; 1.47190];
        p2          = polyfit(alfa_dot2, CL_dot2, 2);

        % CL_interp   = ones(length(numb), 1);

        for i = 1:length(CL_interp)
            if alfa < alfa_star 
                CL_interp = polyval(p1, alfa);
            elseif alfa(i) > alfa_star
                CL_interp = polyval(p2, alfa);
            end
        end

end

