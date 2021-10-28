classdef ShearBendingTorsion
% SHEAR - BENDING - TORSION Class with simple functions.
%   A generic collection of useful function to evaluate Shear and
%   Bending moment distributions along the main wing span. 
%   In this class, one can find: 
%
%    1. A function to calculate a chord distribution, knowing the wing
%       surface Sw, the taper ratio lambda, the span b and a vector
%       which contains stations along the span, called y inside this
%       file.
%    2. Two functions to evaluate axial and normal forces acting on the
%       wing per unit length along the wing span. To be used the
%       function requires a knowledge of the product c(y)*Cl(y) and 
%       c(y)*Cd(y) plus the actual total angle of attack of the wing
%       with respect to the flow direction, defined as 
%       
%                  Alpha_TOT = Alpha(y) + iw(y)
%       
%       with iw the twist angle of the wing.
%    3. Two functions to evaluate shear and bending moment distribution
%       along the span of the wing. Shear force distribution is based
%       on the normal force distribution along the wing, while bending
%       moment is calculated from the shear distribution and a proper
%       arm from the wing root/aircraft longitudinal simmetry axis. 
%    4. A function that automatically plot the shear and bending moment
%       diagram along the wing span. 
%
% It is important to remember the following flight airworthiness rules:
%
% FLIGHT LOADS
%  CS - VLA 321 
%  (a) Flight load factors represent the ratio of the aerodynamic force
%      component (acting normal to the assumed longitudinal axis of the
%      aeroplane) to the weight of the aeroplane. A positive flight load
%      factor is one in which the aerodynamic force acts upward, with
%      respect to the aeroplane. 
%  (b) Compliance with the flight load requirements of this subpart must
%      shown - 
%      (1) At each critical altitude within the range in which the
%          aeroplane may be expected to operate;
%      (2) At each practicable combination of weight and disposable load
%          within the operating limitations specified in the Flight Manual.
%  
%  CS - VLA 331
%  (a) The appropriate balancing horizontal tail load must be accounted for
%      in a rational or conservative manner when determining the wing loads
%      and linear inertia loads corresponding to any of the symmetrical
%      flight conditions specified in CS - VLA 331 to 345.
%  (b) The incremental horizontal tail loads due to manoeuvring and gusts
%      must be reacted by the angular inertia of the aeroplane in a
%      rational or conservative manner.

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
            
            obj = times(cCl, cosd(AoA_Tot)) + times(cCd, sind(AoA_Tot)); 
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
            
            obj = times(cCd, cosd(AoA_Tot)) - times(cCl, sind(AoA_Tot)); 
       end
       function obj = calc_shear_force(obj, y, FZ)
            % S = calc_shear_force(y, FZ)
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
            % BM = calc_bend_mom(y, S)
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
        %% BALANCING LOADS - METHOD
        function obj = Shear_Bending_Torsion_diag(obj, y, Shear, Bend_mom, Torsion, Point)
        % FUNCTION DETAILED DESCRIPTION
        %  fig1 = Bending_Shear_diag(y, Shear, Bend_mom, Torsion, Point)
        %
        %  This function plot the shear and bending moment distributions
        %  along the main wing span. 
        %
        %      INPUTS 
        %       y             --> A vector which contains stations along
        %                         the wing span
        %       Shear         --> A vector which contains shear forces at
        %                         local stations along the wing span  
        %       Bend_mom      --> A vector which contains bending moment at
        %                         local stations along the wing span 
        %       Torsion       --> A vector which contains torsion at local
        %                         stations along the wing span
        %       Point         --> Point of the final envelope diagram
        %
        %      OUTPUT
        %       fig5          --> Shear, bending moment and torsion diagram
        %            
        %      FURTER REFERENCES
        %       The V-N diagram is based on: 
        %         
        %              CS-VLA 333 FLIGHT ENVELOPE, pag. 39/190
        %                 
        %       EASA Airworthiness rules can be find at 
        %        url: https://www.easa.europa.eu/sites/default/files/dfu/Easy Access Rules CS-VLA (Amendment 1).pdf

            obj = figure;
            hold on
            grid on 
            grid minor
            % ylim([min(HT_Lift_posstall)-0.5 max(HT_Lift_fromDtoE)+0.5])
            % xlim([0 max(V_fromDtoE)+10])
            subplot(3,1,1);
            plot(y, Shear, '-r', 'LineWidth', 1)
            grid on 
            grid minor
            % text(y(1)-0.2, Shear(1)+20.0, Reg) % Regulation applied
            % xlabel("Position along the span - $y$ (m)", "Interpreter", "latex")
            ylabel("Shear - $S$ $(daN/m)$", "Interpreter", "latex")
            title('Subplot 1: Shear distribution')

            subplot(3,1,2); 
            plot(y, Bend_mom, '-r', 'LineWidth', 1)
            grid on 
            grid minor
            % text(y(1)-0.2, Bend_mom(1)+20.0, Aircraft_name)                                 % Aircraft name inside the plot
            % xlabel("Position along the span - $y$ (m)", "Interpreter", "latex")
            ylabel("Bending - $M$ $(daN*m)$", "Interpreter", "latex")
            title('Subplot 2: Bending moment distribution')
            
            subplot(3,1,3); 
            plot(y, Torsion, '-r', 'LineWidth', 1)
            grid on 
            grid minor
            % text(y(1)-0.2, Torsion(1)-20.0, Point)                                           % Regulation applied
            xlabel("Position along the semispan - $y$ (m)", "Interpreter", "latex")
            ylabel("Torsion - $T$ $(daN*m)$", "Interpreter", "latex")
            title('Subplot 3: Torsion distribution')
            
            sgtitle(Point) 
        end
        %% BALANCING LOADS - METHOD      
        function obj = calc_cl_distr(obj, CLcentral, CL1, CL2, ...
                                     cl1, cl2)
            
        % FUNCTION DETAILED DESCRIPTION
        %  cl(y) = Bending_Shear_diag(y, Shear, Bend_mom, ...
        %                                    Reg, Aircraft_name)
        %
        %  This function evaluates a lift coefficient along the span
        %  distribution with the following technique:
        %  
        %                   cl_centr = k1*cl1 + k2*cl2 
        %
        %  where k1 and k2 are properly defined weights, in the following
        %  simple manner:
        %                   k1 = CLcentral - CL1
        %                   k2 = CL2 - CLcentral
        % 
        %      INPUTS 
        %       CLcentral --> Value of the lift coefficient at which we
        %                     want to calculate cl = cl(y)
        %       CL1       --> An already calculated lift coefficient value;
        %                     notice that CLcentral > CL1
        %       CL2       --> An already calculated lift coefficient value;
        %                     notice that CLcentral < CL2  
        %       cl1       --> Lift coefficient distribution along the span
        %                     at a global CL = CL1
        %       cl2       --> Lift coefficient distribution along the span
        %                     at a global CL = CL2
        %
        %      OUTPUT
        %       cl        --> Lift coefficient along the span distribution
        %                     at the desired global lift coeff. CLcentral
        %            
        %      FURTER REFERENCES
        %       The V-N diagram is based on: 
        %         
        %              CS-VLA 333 FLIGHT ENVELOPE, pag. 39/190
        %                 
        %       EASA Airworthiness rules can be find at 
        %        url: https://www.easa.europa.eu/sites/default/files/dfu/Easy Access Rules CS-VLA (Amendment 1).pdf
        
        k1 = CLcentral - CL1;
        k2 = CL2 - CLcentral; 
        
        obj = k1*cl1 + k2*cl2;
        
        end
        %% BALANCING LOADS - METHOD      
        function obj = calc_cd_distr(obj, CLcentral, CL1, CL2, ...
                                     cd1, cd2)
            
        % FUNCTION DETAILED DESCRIPTION
        %  cl(y) = Bending_Shear_diag(y, Shear, Bend_mom, ...
        %                                    Reg, Aircraft_name)
        %
        %  This function evaluates a lift coefficient along the span
        %  distribution with the following technique:
        %  
        %                   cl_centr = k1*cl1 + k2*cl2 
        %
        %  where k1 and k2 are properly defined weights, in the following
        %  simple manner:
        %                   k1 = CLcentral - CL1
        %                   k2 = CL2 - CLcentral
        % 
        %      INPUTS 
        %       CLcentral --> Value of the lift coefficient at which we
        %                     want to calculate cl = cl(y)
        %       CL1       --> An already calculated lift coefficient value;
        %                     notice that CLcentral > CL1
        %       CL2       --> An already calculated lift coefficient value;
        %                     notice that CLcentral < CL2  
        %       cd1       --> Drag coefficient distribution along the span
        %                     at a global CL = CL1
        %       cd2       --> Drag coefficient distribution along the span
        %                     at a global CL = CL2
        %
        %      OUTPUT
        %       cd        --> Drag coefficient along the span distribution
        %                     at the desired global lift coeff. CLcentral
        %            
        %      FURTER REFERENCES
        %       The V-N diagram is based on: 
        %         
        %              CS-VLA 333 FLIGHT ENVELOPE, pag. 39/190
        %                 
        %       EASA Airworthiness rules can be find at 
        %        url: https://www.easa.europa.eu/sites/default/files/dfu/Easy Access Rules CS-VLA (Amendment 1).pdf
        
        k1 = CLcentral - CL1;
        k2 = CL2 - CLcentral; 
        
        obj = k1*cd1 + k2*cd2;
        
        end
        %% NORMAL FORCE CALCULATOR 
        function N_distr_along_wing(obj, cCz, qX)
            
        % FUNCTION DETAILED DESCRIPTION
        % Simple normal force calculator. N is a force per unit length.
        % 
        % INPUT 
        % cCz(y) --> A vector which contains product of c(y)*CZ(y) where 
        %            c(y)   --> Chord distribution along the span 
        %            cCz(y) --> Normal force coefficient distribution along
        %                       the span
        %            which multiplied by qX gives a force per unit length
        % qX     --> A vector which contains dynamic pressure per unit of
        %            aerodynamic surface 
        % OUTPUT 
        % N(y)   --> Normal force (per unit length) distribution along the
        %            span 
        
        obj = zeros(length(cCz), 1);
        
            for i = 1:length(cCz)
                obj(i) = cCz(i)*qX;
            end
        
        end
        %% AXIAL FORCE CALCULATOR 
        function A_distr_along_wing(obj, cCa, qX)
         
        % FUNCTION DETAILED DESCRIPTION
        % Simple axial force calculator. A is a force per unit length.
        % 
        % INPUT 
        % cCa(y) --> A vector which contains product of c(y)*CA(y) where 
        %            c(y)   --> Chord distribution along the span 
        %            cCa(y) --> Axial force coefficient distribution along
        %                       the span
        %            which multiplied by qX gives a force per unit length
        % qX     --> A vector which contains dynamic pressure per unit of
        %            aerodynamic surface 
        % OUTPUT 
        % A(y)   --> Axial force (per unit length) distribution along the
        %            span
        
        obj = zeros(length(cCa), 1);
        
            for i = 1:length(cCz)
                obj(i) = cCa(i)*qX;
            end
        
        end
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function obj = calc_tors_mom(obj, y, m)
            
            % TORSION MOMENT DISTRIBUTION
            % T = calc_tors_mom(obj, y, m)
            %   Torsion moment calculator function.
            %   Lift and drag non dimensional force coefficient will be
            %   projected along wing axes; with the axial force is
            %   possible to evaluate shear and bending moment distribution
            %   along the span.
            %   INPUT 
            %     y     --> A vector with distances from the body
            %              longitudinal axis of symmetry
            %     m     --> A vector which contains local pitching mom. 
            %               along the main lifting wing
            %   OUTPUT
            %     T(y)  --> Torsion moment distribution along the span 
            
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
            Tors_mom  = zeros(length(y),1); 
            for i = 2:length(y) 
                a = 0.5*(m(i-1)+m(i))*(y(i-1) - y(i));
                Tors_mom(i) = Tors_mom(i-1) + a; 
            end
            obj = Tors_mom;
        end 
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function obj = Tor_mom_per_unit_length(obj, CM, Dyn_Press, c_of_y)
            
            % TORSION MOMENT DISTRIBUTION
            % m(y) = Tor_mom_per_unit_length(obj, CM, Dyn_Press, c_of_y)
            %   Pitching moment per unit length along the wing span.
            %   Lift and drag non dimensional force coefficient will be
            %   projected along wing axes; with the axial force is
            %   possible to evaluate shear and bending moment distribution
            %   along the span.
            %   INPUT 
            %     CM        --> Pitching moment coefficient distribution 
            %                   per unit length along the wing span 
            %     Dyn_Press --> Dynamic pressure evaluated at the selected 
            %                   flight envelope point
            %     c_of_y    --> Chord distribution per unit length along
            %                   the wing span
            %   OUTPUT
            %     m(y) --> Pitching moment distribution along the wing span
            
            M_of_y = zeros(length(CM), 1);
            for i = 1:length(CM)
                x         = c_of_y(i);
                M_of_y(i) = CM(i)*Dyn_Press*(x^2);
            end
            obj = M_of_y;
        end 
    end
end

