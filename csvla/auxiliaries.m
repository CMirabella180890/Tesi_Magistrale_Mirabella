classdef auxiliaries
    % auxiliaries Class with simple functions 
    %   A generic collection of useful function to evaluate Shear and
    %   Bending moment distributions along the main wing span.
    
    properties
        chord
    end
    
    methods
        function obj = calc_chord(obj, Swing, taper_ratio, span, y)
            % c(y) = calc_chord(Swing, taper_ratio, span, y) 
            %    A simple chord distribution calculator. 
            %
            % INPUT
            %    Swing       --> Planform wing surface
            %    taper_ratio --> Wing taper ratio, defined as the ratio
            %                    Ctip/Croot
            %    span        --> Wing span 
            %    y           --> A vector along the span defined as follow:
            %                     * start: at the root; 
            %                     * end: at the tip; 
            %                     * n: station along the main wing span.
            % OUTPUT 
            %    c = c(y)    --> A chord distribution along the span
            
            Ref_surf = 2*Swing;
            u        = 1 - taper_ratio;
            v        = 1 + taper_ratio; 
            obj      = (Ref_surf/(v*span))*(1 - (u/span)*abs(2*y));
        end
    end
        
    properties
        forces
    end
    
    methods
       function obj = calc_normal_force(obj, AoA_Tot, cCl, cCd)
            % N = calc_normal_force(AoA_Tot, cCl, cCd)
            %   Normal force calculator function.
            %   Lift and drag non dimensional force coefficient will be
            %   projected along wing axes; with the normal force is
            %   possible to evaluate shear and bending moment distribution
            %   along the span.
            %   INPUT 
            %     AoA_Tot --> Total angle of attack defined as the sum of
            %                 the local angle of attack and the wing twist
            %                 angle
            %     cCl     --> Product of the local chord and the local wing
            %                 lift coefficient
            %     cCd     --> Product of the local chord and the local wing
            %                 drag coefficient
            %   OUTPUT
            %     N(y)    --> This is the aerodynamic force component
            %                 normal to the plane of the wing; it is
            %                 measured in [N * m^-1], so it is a force per
            %                 unit length
            
            obj = cCl*cos(AoA_Tot) + cCd*sin(AoA_Tot); 
       end
       function obj = calc_axial_force(obj, AoA_Tot, cCl, cCd)
            % A = calc_normal_force(AoA_Tot, cCl, cCd)
            %   Axial force calculator function.
            %   Lift and drag non dimensional force coefficient will be
            %   projected along wing axes; with the axial force is
            %   possible to evaluate shear and bending moment distribution
            %   along the span.
            %   INPUT 
            %     AoA_Tot --> Total angle of attack defined as the sum of
            %                 the local angle of attack and the wing twist
            %                 angle
            %     cCl     --> Product of the local chord and the local wing
            %                 lift coefficient
            %     cCd     --> Product of the local chord and the local wing
            %                 drag coefficient
            %   OUTPUT
            %     A(y)    --> This is the aerodynamic force component
            %                 which in the plane of the wing; it is
            %                 measured in [N * m^-1], so it is a force per
            %                 unit length
            
            obj = cCd*cos(AoA_Tot) - cCl*sin(AoA_Tot); 
       end
       function obj = calc_shear_force(obj, y, FZ)
            % A = calc_normal_force(AoA_Tot, cCl, cCd)
            %   Axial force calculator function.
            %   Lift and drag non dimensional force coefficient will be
            %   projected along wing axes; with the axial force is
            %   possible to evaluate shear and bending moment distribution
            %   along the span.
            %   INPUT 
            %     y    --> A vector with distances from the body
            %              longitudinal axis of symmetry
            %     FZ   --> A vector which contains normal forces at the 
            %              local station y
            %   OUTPUT
            %     S(y) --> Shear force distribution along the span 
            
            array_dim = size(y);
            if array_dim == [1, length(y)]
                disp("Dimension must be correct");
                y  = y'; 
                y  = flip(y); 
            elseif array_dim == [length(y), 1]
                disp("Vector already in correct dimension");
                y = flip(y);
            end
            
            FZ = flip(FZ);
            S  = zeros(length(y),1); 
            for i = 2:length(y) 
                a = 0.5*(FZ(i-1)+FZ(i))*(y(i-1) - y(i));
                S(i) = S(i-1) + a; 
            end
            obj = S;
       end   
       function obj = calc_bend_mom(obj, y, S)
            % A = calc_normal_force(AoA_Tot, cCl, cCd)
            %   Axial force calculator function.
            %   Lift and drag non dimensional force coefficient will be
            %   projected along wing axes; with the axial force is
            %   possible to evaluate shear and bending moment distribution
            %   along the span.
            %   INPUT 
            %     y     --> A vector with distances from the body
            %              longitudinal axis of symmetry
            %     S     --> A vector which contains local shear forces
            %               along the main lifting wing
            %   OUTPUT
            %     BM(y) --> Bending moment distribution along the span 
            
            array_dim = size(y);
            if array_dim == [1, length(y)]
                disp("Dimension must be correct");
                y  = y'; 
                y  = flip(y); 
            elseif array_dim == [length(y), 1]
                disp("Vector already in correct dimension");
                y = flip(y);
            end
            
            % S = flip(S);
            Bend_mom  = zeros(length(y),1); 
            for i = 2:length(y) 
                a = 0.5*(S(i-1)+S(i))*(y(i-1) - y(i));
                Bend_mom(i) = Bend_mom(i-1) + a; 
            end
            obj = Bend_mom;
       end         
    end
end

