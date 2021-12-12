%% Script to evaluate shear and bending moment distr. along the main wing span

%   DESCRIPTION
%    In this script, the lift and drag distribution along the span will be 
%    used to evaluate shear and bending moment along the main wing span.
%    First, it is necessary to store the critical point, which can be
%    exctracted from the V-N diagram. 

%% STORE CRITICAL POINTS OF THE FINAL ENVELOPE
% It is crucial to store inside the struct variable 'Aircraft' all the
% critical points coming from the Final envelope diagram; this points are
% associated with the most critical, in-flight structural loads for the 
% main wing and horizontal tail structures. These loads are used to size
% structural elements (i.e. main wing spar, main wing ribs...) of the
% aircraft. The points are (V, n) couples.
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Interpolation to obtain Cd at the desired global CL
x_cl   = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value(7,:));
y_cl   = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value(:,1));
xi_cl  = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value(7,:));
yi_cl  = 1:0.1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value(:,1));
[XX_cl,YY_cl] = meshgrid(x_cl, y_cl);
[XI_cl,YI_cl] = meshgrid(xi_cl, yi_cl);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value = interp2(XX_cl, YY_cl, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value, ...
                                                                        XI_cl, YI_cl, 'spline');
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.Attributes.unit = "Non dimensional";
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% FIGURE 8 - INTERPOLATION OF SPANWISE LIFT DISTRIBUCTION
disp(" ")
disp(" ++++ FIGURE 8 - 3D INTERPOLATION OF SPANWISE LIFT DISTR. ++++ ");
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.Attributes.graph = cl_interpolation_graph(x_cl, y_cl, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value, ... 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value, XI_cl, YI_cl); 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.Attributes.graph, 'ClInterpolation3dplot.pdf', 'ContentType', 'vector')
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.Attributes.graph, 'ClInterpolation3dplot.png', 'ContentType', 'vector')
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Saving figures inside correct folder
fprintf('Saving ClInterpolation3dplot.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile ClInterpolation3dplot.pdf Output
movefile ClInterpolation3dplot.png Output
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INTERPOLATION OF THE DRAG AND PITCHING MOMENT COEFFICIENTS
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Interpolation to obtain Cd at the desired global CL
x   = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:));
y   = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(:,1));
xi  = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:));
yi  = 1:0.1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(:,1));
[XX,YY] = meshgrid(x,y);
[XI,YI] = meshgrid(xi,yi);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Interpolated_Cd.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value, ...
                                                                        XI, YI, 'spline');
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Interpolated_Cd.Attributes.unit = 'Non dimensional';
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% FIGURE 9 - INTERPOLATION OF SPANWISE DRAG DISTRIBUCTION
disp(" ")
disp(" ++++ FIGURE 9 - 3D INTERPOLATION OF SPANWISE LIFT DISTR. ++++ ");
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.CD_Interpolation_Graph.value = cd_interpolation_graph(x, y, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value, ... 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Interpolated_Cd.value, XI, YI);            
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.CD_Interpolation_Graph.value, 'CdInterpolation3dplot.pdf', 'ContentType', 'vector')
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.CD_Interpolation_Graph.value, 'CdInterpolation3dplot.png', 'ContentType', 'vector')
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Saving figures inside correct folder
fprintf('Saving CdInterpolation3dplot.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile CdInterpolation3dplot.pdf Output
movefile CdInterpolation3dplot.png Output            
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
% Pitching moment interpolation along the span
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Interpolated_Cm.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value, ...
                                                                        XI, YI, 'spline'); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Interpolated_Cm.Attributes.unit = 'Non dimensional';
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% FIGURE 10 - INTERPOLATION OF SPANWISE DRAG DISTRIBUCTION
disp(" ")
disp(" ++++ FIGURE 10 - 3D INTERPOLATION OF SPANWISE LIFT DISTR. ++++ ");
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.CM_Interpolation_Graph.value = cm_interpolation_graph(x, y, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value, ... 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Interpolated_Cm.value, XI, YI);
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.CM_Interpolation_Graph.value, 'CmInterpolation3dplot.pdf', 'ContentType', 'vector')
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.CM_Interpolation_Graph.value, 'CmInterpolation3dplot.png', 'ContentType', 'vector')
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Saving figures inside correct folder
fprintf('Saving CmInterpolation3dplot.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile CmInterpolation3dplot.pdf Output
movefile CmInterpolation3dplot.png Output
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Half-span stations vector 
b      = Aircraft.Geometry.Wing.b.value;
b_half = b*0.5;
N      = length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value(1,:));
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Half-span vector
half_span = linspace(0, b_half, N);
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% STARTING CALCULATIONS OF SHEAR, BENDING, TORSION

% ---------------- **** **** STRAIGHT FLIGHT **** **** ----------------
switch (Straight_flight_Case)
    % CASE 1: VA greater than the intercept
    case 'Case 1'
        % =================================================================
        if max(n_gust_cruise_plus) > nmax
        % ================================================================= 
        
            % STORE USEFUL VALUES INSIDE THE AIRCRAFT STRUCT VARIABLE 
            % Point S
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.yS.value = half_span;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.yS.Attributes.unit = "m";
            % Point A1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.yA1.value = half_span;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.yA1.Attributes.unit = "m";  
            % Point C1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.yC1.value = half_span;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.yC1.Attributes.unit = "m";
            % Point C
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.yC.value = half_span;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.yC.Attributes.unit = "m";
            % Point C2
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.yC2.value = half_span;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.yC2.Attributes.unit = "m";
            % Point D
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.yD.value = half_span;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.yD.Attributes.unit = "m";

            % LIFT CURVE AND LIFT COEFFICIENT EVALUATED AT THE FINAL ENVELOPE POINTS
            alpha_dot = Aircraft.Certification.Aerodynamic_data.alpha.value;
            CL_dot    = Aircraft.Certification.Aerodynamic_data.CL.value;
            disp(" ")
            disp(" ++++ FIGURE 11 - LIFT MODELS AND FLIGHT ENVELOPE POINTS ++++ ");
            % ---------------------------------------------------------------------------------------------
            CL_S  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value;
            CL_A1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.value;
            CL_C  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value;
            CL_D  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value;
            % ---------------------------------------------------------------------------------------------
            alfaS  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.value;
            alfaA1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.value;
            alfaC  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value;
            alfaD  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value;
            % ---------------------------------------------------------------------------------------------
            PointS  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value;
            PointA1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.point_name.value;
            PointC  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value;
            PointD  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value;
            % ---------------------------------------------------------------------------------------------
            AOA_aux_fullmodel = Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.value;
            CL_fullmodel      = Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value;
            % ---------------------------------------------------------------------------------------------
            Diagram_lift_coefficient_comparison = Lift_coefficients_Points(alfaS, ...
                                       alfaA1, alfaC, alfaD, ...
                                       CL_S, CL_A1, CL_C, CL_D, ...
                                       PointS, PointA1, PointC, PointD, ...
                                       str2num(alpha_dot), str2num(CL_dot), ...
                                       AOA_aux_fullmodel, CL_fullmodel);

            exportgraphics(Diagram_lift_coefficient_comparison, 'LiftComparisonWithPoints.pdf', 'ContentType', 'vector')
            exportgraphics(Diagram_lift_coefficient_comparison, 'LiftComparisonWithPoints.png', 'ContentType', 'vector')

            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Diagram_lift_coefficient_comparison.value = Diagram_lift_coefficient_comparison;

            % Saving figures inside correct folder
            fprintf('Saving LiftComparisonWithPoints.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile LiftComparisonWithPoints.pdf Output
            movefile LiftComparisonWithPoints.png Output    
            
            % LOAD A CLASS OF FUNCTIONS USEFUL TO EVALUATE ALL THE REQUIRED DATA
            obj2 = ShearBendingTorsion; 

            % Lift coefficient distribution at a global CL equal to one
            % A simple function to evaluate the lift coefficient distribution along the
            % span cl = cl(y) when the associated global lift coefficient of the whole
            % wing is equal to 1.0; the function use a method similar to that suggested
            % by Abbott in Theory of Wing Section. See the complete documentation
            % inside the cl_unit_lift.m file
            CL_equal_to_one = 0.0;
            global_CL = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value(1,:)), 1);
            for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value(1,:))
                global_CL(i) = trapz(half_span, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value(i,:))/S;
                if (global_CL(i) >= 1.0-1e-2) && (global_CL(i) <= 1.0+1e-2)
                        CL_equal_to_one = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value(i,:)';
                end
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value = CL_equal_to_one;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.Attributes.unit = "Non dimensional"; 

            % Main wing taper ratio
            ctip        = Aircraft.Geometry.Wing.ctip.value;
            croot       = Aircraft.Geometry.Wing.croot.value;
            taper_ratio = ctip/croot;
            % STORE INSIDE THE AIRCRAFT STRUCT VARIABLE
            Aircraft.Geometry.Wing.taper.value = taper_ratio;
            Aircraft.Geometry.Wing.taper.Attributes.unit = "Non dimensional";

            % Calculation of a chord distribution with a convenient, simple function.
            % 
            % c(y) = calc_chord(Swing, taper_ratio, span, y)
            % A complete documentation of this function is included inside the class
            % ShearBendingTorsion.m
            S           = Aircraft.Geometry.Wing.S.value;
            b           = Aircraft.Geometry.Wing.b.value;
            chord_distr = calc_chord(obj2, S, taper_ratio, b, half_span);
            % STORE INSIDE THE AIRCRAFT STRUCT VARIABLE
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value = chord_distr';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.Attributes.unit = "m";            
        
        % =================================================================
        elseif max(n_gust_cruise_plus) < nmax
        % ================================================================= 
        
            % STORE USEFUL VALUES INSIDE THE AIRCRAFT STRUCT VARIABLE 
            % Point S
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.yS.value = half_span;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.yS.Attributes.unit = "m";
            % Point A1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.yA1.value = half_span;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.yA1.Attributes.unit = "m";
            % Point C
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.yC.value = half_span;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.yC.Attributes.unit = "m";
            % Point D
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.yD.value = half_span;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.yD.Attributes.unit = "m";

            % LIFT CURVE AND LIFT COEFFICIENT EVALUATED AT THE FINAL ENVELOPE POINTS
            alpha_dot = Aircraft.Certification.Aerodynamic_data.alpha.value;
            CL_dot    = Aircraft.Certification.Aerodynamic_data.CL.value;
            disp(" ")
            disp(" ++++ FIGURE 11 - LIFT MODELS AND FLIGHT ENVELOPE POINTS ++++ ");
            % ---------------------------------------------------------------------------------------------
            CL_S  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value;
            CL_A1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.value;
            CL_C  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value;
            CL_D  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value;
            % ---------------------------------------------------------------------------------------------
            alfaS  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.value;
            alfaA1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.value;
            alfaC  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value;
            alfaD  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value;
            % ---------------------------------------------------------------------------------------------
            PointS = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value;
            PointA1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.point_name.value;
            PointC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value;
            PointD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value;
            % ---------------------------------------------------------------------------------------------
            AOA_aux_fullmodel = Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.value;
            CL_fullmodel      = Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value;
            % ---------------------------------------------------------------------------------------------
            Diagram_lift_coefficient_comparison = Lift_coefficients_Points(alfaS, ...
                                       alfaA1, alfaC, alfaD, ...
                                       CL_S, CL_A1, CL_C, CL_D, ...
                                       PointS, PointA1, PointC, PointD, ...
                                       str2num(alpha_dot), str2num(CL_dot), ...
                                       AOA_aux_fullmodel, CL_fullmodel);

            exportgraphics(Diagram_lift_coefficient_comparison, 'LiftComparisonWithPoints.pdf', 'ContentType', 'vector')
            exportgraphics(Diagram_lift_coefficient_comparison, 'LiftComparisonWithPoints.png', 'ContentType', 'vector')

            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Diagram_lift_coefficient_comparison.value = Diagram_lift_coefficient_comparison;

            % Saving figures inside correct folder
            fprintf('Saving LiftComparisonWithPoints.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile LiftComparisonWithPoints.pdf Output
            movefile LiftComparisonWithPoints.png Output        
        
            % LOAD A CLASS OF FUNCTIONS USEFUL TO EVALUATE ALL THE REQUIRED DATA
            obj2 = ShearBendingTorsion; 

            % Lift coefficient distribution at a global CL equal to one
            % A simple function to evaluate the lift coefficient distribution along the
            % span cl = cl(y) when the associated global lift coefficient of the whole
            % wing is equal to 1.0; the function use a method similar to that suggested
            % by Abbott in Theory of Wing Section. See the complete documentation
            % inside the cl_unit_lift.m file
            CL_equal_to_one = 0.0;
            global_CL = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value(1,:)), 1);
            for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value(1,:))
                global_CL(i) = trapz(half_span, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value(i,:))/S;
                if (global_CL(i) >= 1.0-1e-2) && (global_CL(i) <= 1.0+1e-2)
                        CL_equal_to_one = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value(i,:)';
                end
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value = CL_equal_to_one;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.Attributes.unit = "Non dimensional"; 

            % Main wing taper ratio
            ctip        = Aircraft.Geometry.Wing.ctip.value;
            croot       = Aircraft.Geometry.Wing.croot.value;
            taper_ratio = ctip/croot;
            % STORE INSIDE THE AIRCRAFT STRUCT VARIABLE
            Aircraft.Geometry.Wing.taper.value = taper_ratio;
            Aircraft.Geometry.Wing.taper.Attributes.unit = "Non dimensional";

            % Calculation of a chord distribution with a convenient, simple function.
            % 
            % c(y) = calc_chord(Swing, taper_ratio, span, y)
            % A complete documentation of this function is included inside the class
            % ShearBendingTorsion.m
            S           = Aircraft.Geometry.Wing.S.value;
            b           = Aircraft.Geometry.Wing.b.value;
            chord_distr = calc_chord(obj2, S, taper_ratio, b, half_span);
            % STORE INSIDE THE AIRCRAFT STRUCT VARIABLE
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value = chord_distr';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.Attributes.unit = "m";
        
        end
    % CASE 1: VA lower than the intercept
    case 'Case 2'  
        % STORE USEFUL VALUES INSIDE THE AIRCRAFT STRUCT VARIABLE 
        % Point S
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.yS.value = half_span;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.yS.Attributes.unit = "m";
        % Point A
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.yA.value = half_span;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.yA.Attributes.unit = "m";
        % Point A1
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.yA1.value = half_span;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.yA1.Attributes.unit = "m";
        % Point C
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.yC.value = half_span;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.yC.Attributes.unit = "m";
        % Point A2
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.yA2.value = half_span;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.yA2.Attributes.unit = "m";
        % Point D
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.yD.value = half_span;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.yD.Attributes.unit = "m";
        
        % LIFT CURVE AND LIFT COEFFICIENT EVALUATED AT THE FINAL ENVELOPE POINTS
        alpha_dot = Aircraft.Certification.Aerodynamic_data.alpha.value;
        CL_dot    = Aircraft.Certification.Aerodynamic_data.CL.value;
        disp(" ")
        disp(" ++++ FIGURE 11 - LIFT MODELS AND FLIGHT ENVELOPE POINTS ++++ ");
        % ---------------------------------------------------------------------------------------------
        CL_S = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value;
        CL_A = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value;
        CL_C = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value;
        CL_D = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value;
        % ---------------------------------------------------------------------------------------------
        alfaS  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.value;
        alfaA  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.value;
        alfaA1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.value;
        alfaC  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value;
        alfaD  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value;
        % ---------------------------------------------------------------------------------------------
        PointS  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value;
        PointA  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.point_name.value;
        PointA1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.point_name.value;
        PointC  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value;
        PointD  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value;
        % ---------------------------------------------------------------------------------------------
        AOA_aux_fullmodel = Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.value;
        CL_fullmodel      = Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value;
        % ---------------------------------------------------------------------------------------------
        Diagram_lift_coefficient_comparison = Lift_coefficients_Points(alfaS, ...
                                   alfaA, alfaC, alfaD, ...
                                   CL_S, CL_A, CL_C, CL_D, ...
                                   PointS, PointA, PointC, PointD, ...
                                   str2num(alpha_dot), str2num(CL_dot), ...
                                   AOA_aux_fullmodel, CL_fullmodel);

        exportgraphics(Diagram_lift_coefficient_comparison, 'LiftComparisonWithPoints.pdf', 'ContentType', 'vector')
        exportgraphics(Diagram_lift_coefficient_comparison, 'LiftComparisonWithPoints.png', 'ContentType', 'vector')

        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Diagram_lift_coefficient_comparison.value = Diagram_lift_coefficient_comparison;
        
        % Saving figures inside correct folder
        fprintf('Saving LiftComparisonWithPoints.pdf in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile LiftComparisonWithPoints.pdf Output
        movefile LiftComparisonWithPoints.png Output
        
        % LOAD A CLASS OF FUNCTIONS USEFUL TO EVALUATE ALL THE REQUIRED DATA
        obj2 = ShearBendingTorsion; 

        % Lift coefficient distribution at a global CL equal to one
        % A simple function to evaluate the lift coefficient distribution along the
        % span cl = cl(y) when the associated global lift coefficient of the whole
        % wing is equal to 1.0; the function use a method similar to that suggested
        % by Abbott in Theory of Wing Section. See the complete documentation
        % inside the cl_unit_lift.m file
        CL_equal_to_one = 0.0;
        global_CL = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value(1,:)), 1);
        for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value(1,:))
            global_CL(i) = trapz(half_span, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value(i,:))/S;
            if (global_CL(i) >= 1.0-1e-2) && (global_CL(i) <= 1.0+1e-2)
                    CL_equal_to_one = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value(i,:)';
            end
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value = CL_equal_to_one;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.Attributes.unit = "Non dimensional"; 
        
        % Main wing taper ratio
        ctip        = Aircraft.Geometry.Wing.ctip.value;
        croot       = Aircraft.Geometry.Wing.croot.value;
        taper_ratio = ctip/croot;
        % STORE INSIDE THE AIRCRAFT STRUCT VARIABLE
        Aircraft.Geometry.Wing.taper.value = taper_ratio;
        Aircraft.Geometry.Wing.taper.Attributes.unit = "Non dimensional";
        
        % Calculation of a chord distribution with a convenient, simple function.
        % 
        % c(y) = calc_chord(Swing, taper_ratio, span, y)
        % A complete documentation of this function is included inside the class
        % ShearBendingTorsion.m
        S           = Aircraft.Geometry.Wing.S.value;
        b           = Aircraft.Geometry.Wing.b.value;
        chord_distr = calc_chord(obj2, S, taper_ratio, b, half_span)';
        % STORE INSIDE THE AIRCRAFT STRUCT VARIABLE
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value = chord_distr;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.Attributes.unit = "m";        
        % =================================================================
        %% POINT S CALCULATIONS                 
        % Lift coefficient distribution along the span at the Point S
        cl_S = CL_S*CL_equal_to_one; % LIFT COEFFICIENT TIMES LIFT DISTRIBUTION ALONG THE SEMI-SPAN
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cl_S.value = cl_S;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cl_S.Attributes.unit = "Non dimensional";
        
        % Drag coefficient ditribution along the span at the Point S (close to
        % stall)
        cd_S = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:)';
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cd_S.value = cd_S;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cd_S.Attributes.unit = "Non dimensional";
        
        % Pitching moment coefficient distribution along the span at Point S
        cm_S = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value(7,:)';
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cm_S.value = cm_S;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cm_S.Attributes.unit = "Non dimensional";
        
        % PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
        % In this section of the code two vectors are defined to store the product 
        % c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
        cCl_distr_S = times(chord_distr, cl_S);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCl_distr.value = cCl_distr_S;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCl_distr.Attributes.unit = "m";
        cCd_distr_S = times(chord_distr, cd_S);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCd_distr.value = cCd_distr_S;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCd_distr.Attributes.unit = "m";
                
        % AoA_Tot = AoA + Twist_angle of the main wing
        twist_angle = Aircraft.Geometry.Wing.twist_angle.value;
        AoA_Tot_S     = alfaS + twist_angle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_deg.value = AoA_Tot_S;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_deg.Attributes.unit = "Degrees"; 
        
        % Convert in radians 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_rad.value = deg2rad(AoA_Tot_S);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_rad.Attributes.unit = "Radians";   
        
        % Calculation of the normal force coefficient
        % N = calc_normal_force(AoA_Tot, cCl, cCd)
        % This function will be used to evaluate the normal force coefficients
        % distribution along the span; it is possible to fin a complete
        % documentation inside the class file ShearBendingTorsion.m 
        cCz_S = calc_normal_force(obj2, AoA_Tot_S, cCl_distr_S, cCd_distr_S);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCz.value = cCz_S;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCz.Attributes.unit = "m"; 
        
        % Calculation of the axial force coefficient 
        % A = calc_normal_force(AoA_Tot, cCl, cCd)
        % This function will be used to evaluate the axial force coefficients
        % distribution along the span; it is possible to fin a complete
        % documentation inside the class file ShearBendingTorsion.m 
        cCa_S = calc_axial_force(obj2, AoA_Tot_S, cCl_distr_S, cCd_distr_S);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCa.value = cCa_S;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCa.Attributes.unit = "m"; 
        
        % Normal force 
        % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
        qS             = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.value;
        Normal_force_S = cCz_S * qS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Normal_force.value = Normal_force_S;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Normal_force.Attributes.unit = "N/m";
        
        % Axial force 
        % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
        Axial_force_S = cCa_S * qS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Axial_force.value = Axial_force_S;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Axial_force.Attributes.unit = "N/m";        
        
        % SHEAR FORCE CALCULATION 
        % A = calc_shear_force(AoA_Tot, y, cCZ)
        % A complete description of this function is available inside the class
        % file ShearBendingTorsion.m 
        Shear_distr_S = calc_shear_force(obj2, half_span, Normal_force_S)*(1e-1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_distr.value = Shear_distr_S;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_distr.Attributes.unit = "daN";
        
        % BENDING MOMENT CALCULATION 
        % BM = calc_bend_mom(y, S)
        % A complete description of this function is included inside the class file
        % ShearBendingTorsion.m
        Bend_mom_distr_S = calc_bend_mom(obj2, half_span, Shear_distr_S);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Bend_mom_distr.value = Bend_mom_distr_S;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Bend_mom_distr.Attributes.unit = "daN*m";
        
        % PLANS FOR THE STRUCTURAL DIMENSIONING
        % To correctly size the aerostructures of the main lifting surface 
        % it is necessary to apply the procedure just developed to the 
        % critical points coming from the V-N diagram. Those point represents 
        % the most demanding flight conditions that our aircraft could survive. 
        % Those points are all stored inside: 
        % --> Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope
        % Retrieve the values and apply formulas to them. 
        
        % Pitching moment per unit length
        m_distr_S = zeros(length(cm_S), 1);
        for i = 1:length(cm_S)
            m_distr_S(i) = cm_S(i) * qS *((chord_distr(i))^2);
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.m_distr.value = m_distr_S;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.m_distr.Attributes.unit = "N";   

        % Torque applied
        % T = calc_tors_mom(obj, y, m)
        % A complete distribution of this function is included inside the class
        % file ShearBendingTorsion.m
        Tors_mom_distr_S = calc_tors_mom(obj2, half_span, m_distr_S)*(1e-1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Tors_mom_distr.value = Tors_mom_distr_S;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Tors_mom_distr.Attributes.unit = "daN*m";   
        
        disp(" ")
        disp(" ++++ FIGURE 12 - POINT S SHEAR, BENDING, TORSION ++++ ");
        Shear_BendMom_diagram_S = Shear_Bending_Torsion_diag(obj2, flip(half_span), Shear_distr_S, Bend_mom_distr_S, ...
                                                           Tors_mom_distr_S, PointS);
        
        exportgraphics(Shear_BendMom_diagram_S, 'ShearBendingTorsionDiagramPointS.pdf', 'ContentType', 'vector')
        exportgraphics(Shear_BendMom_diagram_S, 'ShearBendingTorsionDiagramPointS.png', 'ContentType', 'vector')
        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_BendMom_diagram.value = Shear_BendMom_diagram_S;
        % Saving figures inside correct folder
        fprintf('Saving ShearBendingTorsionDiagramPointS.pdf in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile ShearBendingTorsionDiagramPointS.pdf Output
        movefile ShearBendingTorsionDiagramPointS.png Output        
        % =================================================================
        %% POINT A CALCULATIONS                 
        % Lift coefficient distribution along the span at the Point A
        cl_A = CL_A*CL_equal_to_one; % LIFT COEFFICIENT TIMES LIFT DISTRIBUTION ALONG THE SEMI-SPAN
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cl_A.value = cl_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cl_A.Attributes.unit = "Non dimensional";
        
        % Drag coefficient ditribution along the span at the Point A (close to
        % stall)
        cd_A = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:)';
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cd_A.value = cd_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cd_A.Attributes.unit = "Non dimensional";
        
        % Pitching moment coefficient distribution along the span at Point A
        cm_A = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value(7,:)';
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cm_A.value = cm_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cm_A.Attributes.unit = "Non dimensional";
        
        % PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
        % In this section of the code two vectors are defined to store the product 
        % c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
        cCl_distr_A = times(chord_distr, cl_A);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCl_distr.value = cCl_distr_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCl_distr.Attributes.unit = "m";
        cCd_distr_A = times(chord_distr, cd_A);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCd_distr.value = cCd_distr_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCd_distr.Attributes.unit = "m";
                
        % AoA_Tot = AoA + Twist_angle of the main wing
        AoA_Tot_A     = alfaA + twist_angle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_deg.value = AoA_Tot_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_deg.Attributes.unit = "Degrees"; 
        
        % Convert in radians 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_rad.value = deg2rad(AoA_Tot_A);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_rad.Attributes.unit = "Radians";   
        
        % Calculation of the normal force coefficient
        % N = calc_normal_force(AoA_Tot, cCl, cCd)
        % This function will be used to evaluate the normal force coefficients
        % distribution along the span; it is possible to fin a complete
        % documentation inside the class file ShearBendingTorsion.m 
        cCz_A = calc_normal_force(obj2, AoA_Tot_A, cCl_distr_A, cCd_distr_A);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCz.value = cCz_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCz.Attributes.unit = "m"; 
        
        % Calculation of the axial force coefficient 
        % A = calc_normal_force(AoA_Tot, cCl, cCd)
        % This function will be used to evaluate the axial force coefficients
        % distribution along the span; it is possible to fin a complete
        % documentation inside the class file ShearBendingTorsion.m 
        cCa_A = calc_axial_force(obj2, AoA_Tot_A, cCl_distr_A, cCd_distr_A);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCa.value = cCa_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCa.Attributes.unit = "m"; 
        
        % Normal force 
        % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
        qA             = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        Normal_force_A = cCz_A * qA;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Normal_force.value = Normal_force_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Normal_force.Attributes.unit = "N/m";
        
        % Axial force 
        % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
        Axial_force_A = cCa_A * qA;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Axial_force.value = Axial_force_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Axial_force.Attributes.unit = "N/m";        
        
        % SHEAR FORCE CALCULATION 
        % A = calc_shear_force(AoA_Tot, y, cCZ)
        % A complete description of this function is available inside the class
        % file ShearBendingTorsion.m 
        Shear_distr_A = calc_shear_force(obj2, half_span, Normal_force_A)*(1e-1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_distr.value = Shear_distr_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_distr.Attributes.unit = "daN";
        
        % BENDING MOMENT CALCULATION 
        % BM = calc_bend_mom(y, S)
        % A complete description of this function is included inside the class file
        % ShearBendingTorsion.m
        Bend_mom_distr_A = calc_bend_mom(obj2, half_span, Shear_distr_A);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Bend_mom_distr.value = Bend_mom_distr_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Bend_mom_distr.Attributes.unit = "daN*m";
        
        % PLANS FOR THE STRUCTURAL DIMENSIONING
        % To correctly size the aerostructures of the main lifting surface 
        % it is necessary to apply the procedure just developed to the 
        % critical points coming from the V-N diagram. Those point represents 
        % the most demanding flight conditions that our aircraft could survive. 
        % Those points are all stored inside: 
        % --> Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope
        % Retrieve the values and apply formulas to them. 
        
        % Pitching moment per unit length
        m_distr_A = zeros(length(cm_A), 1);
        for i = 1:length(cm_A)
            m_distr_A(i) = cm_A(i) * qA *((chord_distr(i))^2);
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.m_distr.value = m_distr_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.m_distr.Attributes.unit = "N";   

        % Torque applied
        % T = calc_tors_mom(obj, y, m)
        % A complete distribution of this function is included inside the class
        % file ShearBendingTorsion.m
        Tors_mom_distr_A = calc_tors_mom(obj2, half_span, m_distr_A)*(1e-1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Tors_mom_distr.value = Tors_mom_distr_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Tors_mom_distr.Attributes.unit = "daN*m";   
        
        disp(" ")
        disp(" ++++ FIGURE 13 - POINT A SHEAR, BENDING, TORSION ++++ ");
        Shear_BendMom_diagram_A = Shear_Bending_Torsion_diag(obj2, flip(half_span), Shear_distr_A, Bend_mom_distr_A, ...
                                                           Tors_mom_distr_A, PointA);
        
        exportgraphics(Shear_BendMom_diagram_A, 'ShearBendingTorsionDiagramPointA.pdf', 'ContentType', 'vector')
        exportgraphics(Shear_BendMom_diagram_A, 'ShearBendingTorsionDiagramPointA.png', 'ContentType', 'vector')
        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_BendMom_diagram.value = Shear_BendMom_diagram_A;
        % Saving figures inside correct folder
        fprintf('Saving ShearBendingTorsionDiagramPointA.pdf in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile ShearBendingTorsionDiagramPointA.pdf Output
        movefile ShearBendingTorsionDiagramPointA.png Output 
        
        %% POINT A1 CALCULATIONS                 
        % Lift coefficient distribution along the span at the Point A1
        cl_A1 = CL_A1 * CL_equal_to_one; % LIFT COEFFICIENT TIMES LIFT DISTRIBUTION ALONG THE SEMI-SPAN
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.cl_A1.value = cl_A1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.cl_A1.Attributes.unit = "Non dimensional";
        
        % Drag coefficient ditribution along the span at the Point A1 (close to
        % stall)
        cd_A1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:)';
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.cd_A1.value = cd_A1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.cd_A1.Attributes.unit = "Non dimensional";
        
        % Pitching moment coefficient distribution along the span at Point A1
        cm_A1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value(7,:)';
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.cm_A1.value = cm_A1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.cm_A1.Attributes.unit = "Non dimensional";
        
        % PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
        % In this section of the code two vectors are defined to store the product 
        % c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
        cCl_distr_A1 = times(chord_distr, cl_A1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.cCl_distr.value = cCl_distr_A1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.cCl_distr.Attributes.unit = "m";
        cCd_distr_A1 = times(chord_distr, cd_A1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.cCd_distr.value = cCd_distr_A1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.cCd_distr.Attributes.unit = "m";
                
        % AoA_Tot = AoA + Twist_angle of the main wing
        AoA_Tot_A1     = alfaA1 + twist_angle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.AoA_Tot_deg.value = AoA_Tot_A1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.AoA_Tot_deg.Attributes.unit = "Degrees"; 
        
        % Convert in radians 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.AoA_Tot_rad.value = deg2rad(AoA_Tot_A1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.AoA_Tot_rad.Attributes.unit = "Radians";   
        
        % Calculation of the normal force coefficient
        % N = calc_normal_force(AoA_Tot, cCl, cCd)
        % This function will be used to evaluate the normal force coefficients
        % distribution along the span; it is possible to fin a complete
        % documentation inside the class file ShearBendingTorsion.m 
        cCz_A1 = calc_normal_force(obj2, AoA_Tot_A1, cCl_distr_A1, cCd_distr_A1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.cCz.value = cCz_A1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.cCz.Attributes.unit = "m"; 
        
        % Calculation of the axial force coefficient 
        % A = calc_normal_force(AoA_Tot, cCl, cCd)
        % This function will be used to evaluate the axial force coefficients
        % distribution along the span; it is possible to fin a complete
        % documentation inside the class file ShearBendingTorsion.m 
        cCa_A1 = calc_axial_force(obj2, AoA_Tot_A1, cCl_distr_A1, cCd_distr_A1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.cCa.value = cCa_A1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.cCa.Attributes.unit = "m"; 
        
        % Normal force 
        % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
        qA1             = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.value;
        Normal_force_A1 = cCz_A1 * qA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Normal_force.value = Normal_force_A1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Normal_force.Attributes.unit = "N/m";
        
        % Axial force 
        % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
        Axial_force_A1 = cCa_A1 * qA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Axial_force.value = Axial_force_A1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Axial_force.Attributes.unit = "N/m";        
        
        % SHEAR FORCE CALCULATION 
        % A = calc_shear_force(AoA_Tot, y, cCZ)
        % A complete description of this function is available inside the class
        % file ShearBendingTorsion.m 
        Shear_distr_A1 = calc_shear_force(obj2, half_span, Normal_force_A1)*(1e-1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Shear_distr.value = Shear_distr_A1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Shear_distr.Attributes.unit = "daN";
        
        % BENDING MOMENT CALCULATION 
        % BM = calc_bend_mom(y, S)
        % A complete description of this function is included inside the class file
        % ShearBendingTorsion.m
        Bend_mom_distr_A1 = calc_bend_mom(obj2, half_span, Shear_distr_A1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Bend_mom_distr.value = Bend_mom_distr_A1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Bend_mom_distr.Attributes.unit = "daN*m";
        
        % PLANS FOR THE STRUCTURAL DIMENSIONING
        % To correctly size the aerostructures of the main lifting surface 
        % it is necessary to apply the procedure just developed to the 
        % critical points coming from the V-N diagram. Those point represents 
        % the most demanding flight conditions that our aircraft could survive. 
        % Those points are all stored inside: 
        % --> Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope
        % Retrieve the values and apply formulas to them. 
        
        % Pitching moment per unit length
        m_distr_A1 = zeros(length(cm_A1), 1);
        for i = 1:length(cm_A1)
            m_distr_A1(i) = cm_A1(i) * qA1 *((chord_distr(i))^2);
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.m_distr.value = m_distr_A1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.m_distr.Attributes.unit = "N";   

        % Torque applied
        % T = calc_tors_mom(obj, y, m)
        % A complete distribution of this function is included inside the class
        % file ShearBendingTorsion.m
        Tors_mom_distr_A1 = calc_tors_mom(obj2, half_span, m_distr_A1)*(1e-1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Tors_mom_distr.value = Tors_mom_distr_A1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Tors_mom_distr.Attributes.unit = "daN*m";   
        
        disp(" ")
        disp(" ++++ FIGURE 14 - POINT A1 SHEAR, BENDING, TORSION ++++ ");
        Shear_BendMom_diagram_A1 = Shear_Bending_Torsion_diag(obj2, flip(half_span), Shear_distr_A1, Bend_mom_distr_A1, ...
                                                           Tors_mom_distr_A1, PointA1);
        
        exportgraphics(Shear_BendMom_diagram_A1, 'ShearBendingTorsionDiagramPointA1.pdf', 'ContentType', 'vector')
        exportgraphics(Shear_BendMom_diagram_A1, 'ShearBendingTorsionDiagramPointA1.png', 'ContentType', 'vector')
        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Shear_BendMom_diagram.value = Shear_BendMom_diagram_A1;
        % Saving figures inside correct folder
        fprintf('Saving ShearBendingTorsionDiagramPointA1.pdf in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile ShearBendingTorsionDiagramPointA1.pdf Output
        movefile ShearBendingTorsionDiagramPointA1.png Output         
        
        %% POINT C CALCULATIONS                 
        % Lift coefficient distribution along the span at the Point C
        cl_C = CL_C * CL_equal_to_one; % LIFT COEFFICIENT TIMES LIFT DISTRIBUTION ALONG THE SEMI-SPAN
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cl_C.value = cl_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cl_C.Attributes.unit = "Non dimensional";
        
        % Drag coefficient ditribution along the span at the Point C (close to
        % stall)
        cd_C = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:)';
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cd_C.value = cd_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cd_C.Attributes.unit = "Non dimensional";
        
        % Pitching moment coefficient distribution along the span at Point C
        cm_C = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value(7,:)';
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.value = cm_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.Attributes.unit = "Non dimensional";
        
        % PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
        % In this section of the code two vectors are defined to store the product 
        % c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
        cCl_distr_C = times(chord_distr, cl_C);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCl_distr.value = cCl_distr_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCl_distr.Attributes.unit = "m";
        cCd_distr_C = times(chord_distr, cd_C);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCd_distr.value = cCd_distr_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCd_distr.Attributes.unit = "m";
                
        % AoA_Tot = AoA + Twist_angle of the main wing
        AoA_Tot_C     = alfaC + twist_angle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_deg.value = AoA_Tot_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_deg.Attributes.unit = "Degrees"; 
        
        % Convert in radians 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_rad.value = deg2rad(AoA_Tot_C);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_rad.Attributes.unit = "Radians";   
        
        % Calculation of the normal force coefficient
        % N = calc_normal_force(AoA_Tot, cCl, cCd)
        % This function will be used to evaluate the normal force coefficients
        % distribution along the span; it is possible to fin a complete
        % documentation inside the class file ShearBendingTorsion.m 
        cCz_C = calc_normal_force(obj2, AoA_Tot_C, cCl_distr_C, cCd_distr_C);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCz.value = cCz_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCz.Attributes.unit = "m"; 
        
        % Calculation of the axial force coefficient 
        % A = calc_normal_force(AoA_Tot, cCl, cCd)
        % This function will be used to evaluate the axial force coefficients
        % distribution along the span; it is possible to fin a complete
        % documentation inside the class file ShearBendingTorsion.m 
        cCa_C = calc_axial_force(obj2, AoA_Tot_C, cCl_distr_C, cCd_distr_C);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCa.value = cCa_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCa.Attributes.unit = "m"; 
        
        % Normal force 
        % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
        qC             = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value;
        Normal_force_C = cCz_C * qC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Normal_force.value = Normal_force_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Normal_force.Attributes.unit = "N/m";
        
        % Axial force 
        % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
        Axial_force_C = cCa_C * qC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Axial_force.value = Axial_force_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Axial_force.Attributes.unit = "N/m";        
        
        % SHEAR FORCE CALCULATION 
        % A = calc_shear_force(AoA_Tot, y, cCZ)
        % A complete description of this function is available inside the class
        % file ShearBendingTorsion.m 
        Shear_distr_C = calc_shear_force(obj2, half_span, Normal_force_C)*(1e-1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_distr.value = Shear_distr_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_distr.Attributes.unit = "daN";
        
        % BENDING MOMENT CALCULATION 
        % BM = calc_bend_mom(y, S)
        % A complete description of this function is included inside the class file
        % ShearBendingTorsion.m
        Bend_mom_distr_C = calc_bend_mom(obj2, half_span, Shear_distr_C);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Bend_mom_distr.value = Bend_mom_distr_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Bend_mom_distr.Attributes.unit = "daN*m";
        
        % PLANS FOR THE STRUCTURAL DIMENSIONING
        % To correctly size the aerostructures of the main lifting surface 
        % it is necessary to apply the procedure just developed to the 
        % critical points coming from the V-N diagram. Those point represents 
        % the most demanding flight conditions that our aircraft could survive. 
        % Those points are all stored inside: 
        % --> Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope
        % Retrieve the values and apply formulas to them. 
        
        % Pitching moment per unit length
        m_distr_C = zeros(length(cm_C), 1);
        for i = 1:length(cm_C)
            m_distr_C(i) = cm_C(i) * qC *((chord_distr(i))^2);
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.m_distr.value = m_distr_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.m_distr.Attributes.unit = "N";   

        % Torque applied
        % T = calc_tors_mom(obj, y, m)
        % A complete distribution of this function is included inside the class
        % file ShearBendingTorsion.m
        Tors_mom_distr_C = calc_tors_mom(obj2, half_span, m_distr_C)*(1e-1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Tors_mom_distr.value = Tors_mom_distr_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Tors_mom_distr.Attributes.unit = "daN*m";   
        
        disp(" ")
        disp(" ++++ FIGURE 14 - POINT C SHEAR, BENDING, TORSION ++++ ");
        Shear_BendMom_diagram_C = Shear_Bending_Torsion_diag(obj2, flip(half_span), Shear_distr_C, Bend_mom_distr_C, ...
                                                           Tors_mom_distr_C, PointC);
        
        exportgraphics(Shear_BendMom_diagram_C, 'ShearBendingTorsionDiagramPointC.pdf', 'ContentType', 'vector')
        exportgraphics(Shear_BendMom_diagram_C, 'ShearBendingTorsionDiagramPointC.png', 'ContentType', 'vector')
        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_BendMom_diagram.value = Shear_BendMom_diagram_C;
        % Saving figures inside correct folder
        fprintf('Saving ShearBendingTorsionDiagramPointC.pdf in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile ShearBendingTorsionDiagramPointC.pdf Output
        movefile ShearBendingTorsionDiagramPointC.png Output  
        
        %% POINT D CALCULATIONS                 
        % Lift coefficient distribution along the span at the Point D
        cl_D = CL_D * CL_equal_to_one; % LIFT COEFFICIENT TIMES LIFT DISTRIBUTION ALONG THE SEMI-SPAN
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cl_D.value = cl_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cl_D.Attributes.unit = "Non dimensional";
        
        % Support variables to interpolate
        x   = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:));
        y   = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(:,1));
        xi  = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:));
        yi  = 1:0.1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(:,1));        
        
        % Selection of the interpolated distribution of CD and CM
        CD_D                     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_D.value;
        Interpolated_Global_CD_D = zeros(length(yi), 1);
        for i = 1:length(yi)
            Interpolated_Global_CD_D(i) = trapz(half_span, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Interpolated_Cd.value(i,:));
            if abs(Interpolated_Global_CD_D(i) - CD_D) < 1e-2
               cd_D = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Interpolated_Cd.value(i,:)';
               cm_D = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Interpolated_Cm.value(i,:)';
            end
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Global_CD.value = Interpolated_Global_CD_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Global_CD.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cd_D.value = cd_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.value = cm_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cd_D.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.Attributes.unit = "Non dimensional";
        
        % PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
        % In this section of the code two vectors are defined to store the product 
        % c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
        cCl_distr_D = times(chord_distr, cl_D);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCl_distr.value = cCl_distr_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCl_distr.Attributes.unit = "m";
        cCd_distr_D = times(chord_distr, cd_D);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCd_distr.value = cCd_distr_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCd_distr.Attributes.unit = "m";
                
        % AoA_Tot = AoA + Twist_angle of the main wing
        AoA_Tot_D     = alfaD + twist_angle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_deg.value = AoA_Tot_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_deg.Attributes.unit = "Degrees"; 
        
        % Convert in radians 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_rad.value = deg2rad(AoA_Tot_D);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_rad.Attributes.unit = "Radians";   
        
        % Calculation of the normal force coefficient
        % N = calc_normal_force(AoA_Tot, cCl, cCd)
        % This function will be used to evaluate the normal force coefficients
        % distribution along the span; it is possible to fin a complete
        % documentation inside the class file ShearBendingTorsion.m 
        cCz_D = calc_normal_force(obj2, AoA_Tot_D, cCl_distr_D, cCd_distr_D);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCz.value = cCz_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCz.Attributes.unit = "m"; 
        
        % Calculation of the axial force coefficient 
        % A = calc_normal_force(AoA_Tot, cCl, cCd)
        % This function will be used to evaluate the axial force coefficients
        % distribution along the span; it is possible to fin a complete
        % documentation inside the class file ShearBendingTorsion.m 
        cCa_D = calc_axial_force(obj2, AoA_Tot_D, cCl_distr_D, cCd_distr_D);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCa.value = cCa_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCa.Attributes.unit = "m"; 
        
        % Normal force 
        % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
        qD             = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value;
        Normal_force_D = cCz_D * qD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Normal_force.value = Normal_force_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Normal_force.Attributes.unit = "N/m";
        
        % Axial force 
        % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
        Axial_force_D = cCa_D * qD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Axial_force.value = Axial_force_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Axial_force.Attributes.unit = "N/m";        
        
        % SHEAR FORCE CALCULATION 
        % A = calc_shear_force(AoA_Tot, y, cCZ)
        % A complete description of this function is available inside the class
        % file ShearBendingTorsion.m 
        Shear_distr_D = calc_shear_force(obj2, half_span, Normal_force_D)*(1e-1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_distr.value = Shear_distr_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_distr.Attributes.unit = "daN";
        
        % BENDING MOMENT CALCULATION 
        % BM = calc_bend_mom(y, S)
        % A complete description of this function is included inside the class file
        % ShearBendingTorsion.m
        Bend_mom_distr_D = calc_bend_mom(obj2, half_span, Shear_distr_D);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Bend_mom_distr.value = Bend_mom_distr_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Bend_mom_distr.Attributes.unit = "daN*m";
        
        % PLANS FOR THE STRUCTURAL DIMENSIONING
        % To correctly size the aerostructures of the main lifting surface 
        % it is necessary to apply the procedure just developed to the 
        % critical points coming from the V-N diagram. Those point represents 
        % the most demanding flight conditions that our aircraft could survive. 
        % Those points are all stored inside: 
        % --> Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope
        % Retrieve the values and apply formulas to them. 
        
        % Pitching moment per unit length
        m_distr_D = zeros(length(cm_D), 1);
        for i = 1:length(cm_C)
            m_distr_D(i) = cm_D(i) * qD *((chord_distr(i))^2);
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.m_distr.value = m_distr_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.m_distr.Attributes.unit = "N";   

        % Torque applied
        % T = calc_tors_mom(obj, y, m)
        % A complete distribution of this function is included inside the class
        % file ShearBendingTorsion.m
        Tors_mom_distr_D = calc_tors_mom(obj2, half_span, m_distr_D)*(1e-1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Tors_mom_distr.value = Tors_mom_distr_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Tors_mom_distr.Attributes.unit = "daN*m";   
        
        disp(" ")
        disp(" ++++ FIGURE 15 - POINT D SHEAR, BENDING, TORSION ++++ ");
        Shear_BendMom_diagram_D = Shear_Bending_Torsion_diag(obj2, flip(half_span), Shear_distr_D, Bend_mom_distr_D, ...
                                                           Tors_mom_distr_D, PointD);
        
        exportgraphics(Shear_BendMom_diagram_D, 'ShearBendingTorsionDiagramPointD.pdf', 'ContentType', 'vector')
        exportgraphics(Shear_BendMom_diagram_D, 'ShearBendingTorsionDiagramPointD.png', 'ContentType', 'vector')
        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_BendMom_diagram.value = Shear_BendMom_diagram_D;
        % Saving figures inside correct folder
        fprintf('Saving ShearBendingTorsionDiagramPointD.pdf in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile ShearBendingTorsionDiagramPointD.pdf Output
        movefile ShearBendingTorsionDiagramPointD.png Output         
        
end

% ---------------- **** **** INVERTED FLIGHT **** **** ----------------
switch (Inverted_flight_Case)
    % CASE 1: Complex solutions of the intercept
    case 'Case 1'
    % CASE 2: Real solutions of the intercept
    case 'Case 2'  
        % STORE USEFUL VALUES INSIDE THE AIRCRAFT STRUCT VARIABLE 
        % Point S_inv
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.yS_inv.value = half_span;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.yS_inv.Attributes.unit = "m";
        % Point G
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.yG.value = half_span;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.yG.Attributes.unit = "m";
        % Point G1
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.yG1.value = half_span;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.yG1.Attributes.unit = "m";
        % Point F
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.yF.value = half_span;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.yF.Attributes.unit = "m";
        % Point G2
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.yG2.value = half_span;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.yG2.Attributes.unit = "m";
        % Point E
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.yE.value = half_span;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.yE.Attributes.unit = "m";
        
end

% % For a more in depth description, look ahead in this file
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value = cl_unit_lift(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CL.value(3), ...
% %                                                                                          Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CL.value(4), ...
% %                                                                                          Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value(3,:)', ...
% %                                                                                          Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value(4,:)'); 

% % Point S of the final envelope
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.y_half_span.value = half_span;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.y_half_span.Attributes.unit = "m";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CD_S.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value(1);  
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CD_S.Attributes.unit = "Non dimensional"; 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alpha_S.value = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alpha_S.Attributes.unit = "Degrees";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTail_S.value = (0.5)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.value^2)* ...
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_positivestall.value(1))*(1E-1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTail_S.Attributes.unit = "daN";
% 
% % Point A of the final envelope
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.y_half_span.value = half_span;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.y_half_span.Attributes.unit = "m";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA_unit_load_factor.value = 0.0;
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value)
%     x = abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value - Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value(i));
%     if x < 1e-2
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA_unit_load_factor.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value(i);
%     end
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA_unit_load_factor.Attributes.unit = "m/s";  
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CD_A.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value(1);  
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CD_A.Attributes.unit = "Non dimensional"; 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alpha_A.value = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alpha_A.Attributes.unit = "Degrees";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTail_A.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_positivestall.value(end);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTail_A.Attributes.unit = "daN";
% 
% % Point C of the final envelope
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.y_half_span.value = half_span;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.y_half_span.Attributes.unit = "m";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CD_C.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value(end);  
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CD_C.Attributes.unit = "Non dimensional"; 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alpha_C.value = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alpha_C.Attributes.unit = "Degrees";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTail_C.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_positivestall.value(end);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTail_C.Attributes.unit = "daN";
% 
% % Point D of the final envelope 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.y_half_span.value = half_span;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.y_half_span.Attributes.unit = "m";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_D.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDtoE.value(1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_D.Attributes.unit = "Non dimensional";  
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alpha_D.value = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alpha_D.Attributes.unit = "Degrees";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTail_D.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromDtoE.value(1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTail_D.Attributes.unit = "daN";
% 
% % Point F of the final envelope
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.y_half_span.value = half_span;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.y_half_span.Attributes.unit = "m";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value(1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.Attributes.unit = "Non dimensional";  
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alpha_F.value = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p); 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alpha_F.Attributes.unit = "Radians";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTail_F.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromFtoE.value(1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTail_F.Attributes.unit = "daN";
% 
% % Point G of the final envelope
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.y_half_span.value = half_span;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.y_half_span.Attributes.unit = "m";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.value = Aircraft.Certification.Aerodynamic_data.CD0.value + (Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value^2)/(pi*Aircraft.Certification.Aerodynamic_data.e.value*Aircraft.Geometry.Wing.AR.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.Attributes.unit = "Non dimensional";  
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alpha_G.value = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alpha_G.Attributes.unit = "degrees"; 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTail_G.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_negativestall.value(end)*Aircraft.Geometry.Wing.S.value*(1e-1); % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromFtoE.value(end);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTail_G.Attributes.unit = "daN"; 
% 
% % Point E of the final envelope
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.y_half_span.value = half_span;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.y_half_span.Attributes.unit = "m";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CD_E.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value(end);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CD_E.Attributes.unit = "Non dimensional";  
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alpha_E.value = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alpha_E.Attributes.unit = "Degrees"; 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTail_E.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromFtoE.value(end);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LTHail_E.Attributes.unit = "daN"; 
% 
% %% PRINT POINT VALUES 
% 
% disp(" ")
% disp(" ++++ POINT OF THE FINAL ENVELOPE ++++");
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.nS.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CD_S.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alpha_S.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTail_S.value];
% disp(" ++++++++++ POINT S OF THE FINAL ENVELOPE ++++++++++ ")
% format = ' %6.6f  %6.6f  %6.6f  %6.6f  %6.6f  %6.6f  %6.6f  %6.6f\n';
% label  = ' V [m/s]    n [g]     q [Pa]      CL        CD        alpha [°]  L [daN]     LH [daN]\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CD_A.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alpha_A.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTail_A.value];
% disp(" ++++++++++ POINT A OF THE FINAL ENVELOPE ++++++++++ ")
% % format = ' %f          %f          %f          %f          %f          %f          %f\n';
% % label  = ' V [m/s]           q [Pa]             CL                CD                alpha [deg°]        L [daN]             LH [daN]\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CD_C.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alpha_C.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTail_C.value];
% disp(" ++++++++++ POINT C OF THE FINAL ENVELOPE ++++++++++ ")
% format = ' %6.6f  %6.6f  %6.6f  %6.6f  %6.6f  %6.6f  %6.6f  %6.6f\n';
% label  = ' V [m/s]    n [g]     q [Pa]       CL        CD        alpha [°]  L [daN]     LH [daN]\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_D.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alpha_D.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTail_D.value];
% disp(" ++++++++++ POINT D OF THE FINAL ENVELOPE ++++++++++ ")
% % format = ' %f          %f          %f         %f          %f          %f            %f\n';
% % label  = ' V [m/s]           q [Pa]             CL                CD                alpha [deg°]        L [daN]             LH [daN]\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alpha_F.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTail_F.value];
% disp(" ++++++++++ POINT F OF THE FINAL ENVELOPE ++++++++++ ")
% format = ' %6.6f  %6.6f  %6.6f  %6.6f  %6.6f  %6.6f  %6.6f  %6.6f\n';
% label  = ' V [m/s]    n [g]      q [Pa]       CL        CD        alpha [°]  L [daN]     LH [daN]\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.nG.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alpha_G.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTail_G.value];
% disp(" ++++++++++ POINT G OF THE FINAL ENVELOPE ++++++++++ ")
% % format = ' %f          %f          %f          %f          %f          %f          %f\n';
% % label  = ' V [m/s]           q [Pa]             CL                CD                alpha [deg°]        L [daN]             LH [daN]\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.VE.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.nE.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CD_E.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alpha_E.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTail_E.value];
% disp(" ++++++++++ POINT E OF THE FINAL ENVELOPE ++++++++++ ")
% % format = ' %f          %f          %f          %f          %f          %f          %f\n';
% % label  = ' V [m/s]           q [Pa]             CL                CD                alpha [deg°]        L [daN]             LH [daN]\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

%% LIFT CURVE AND LIFT COEFFICIENT EVALUATED AT THE FINAL ENVELOPE POINTS
% alpha_dot = Aircraft.Certification.Aerodynamic_data.alpha.value;
% CL_dot    = Aircraft.Certification.Aerodynamic_data.CL.value;
% disp(" ")
% disp(" ++++ FIGURE 11 - LIFT MODELS AND FLIGHT ENVELOPE POINTS ++++ ");
% % ---------------------------------------------------------------------------------------------
% CL_S = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value;
% CL_A = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value;
% CL_C = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value;
% CL_D = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value;
% % ---------------------------------------------------------------------------------------------
% alfaS = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.value;
% alfaA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.value;
% alfaC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value;
% alfaD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value;
% % ---------------------------------------------------------------------------------------------
% PointS = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value;
% PointA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.point_name.value;
% PointC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value;
% PointD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value;
% % ---------------------------------------------------------------------------------------------
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Diagram_lift_coefficient_comparison.value = Lift_coefficients_Points(alfaS, ...
%                            alfaA, alfaC, alfaD, ...
%                            CL_S, CL_A, CL_C, CL_D, ...
%                            PointS, PointA, PointC, PointD, ...
%                            str2num(alpha_dot), str2num(CL_dot), ...
%                            Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.value, ...
%                            Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value);
% 
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Diagram_lift_coefficient_comparison.value, 'LiftComparisonWithPoints.pdf', 'ContentType', 'vector')
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Diagram_lift_coefficient_comparison.value, 'LiftComparisonWithPoints.png', 'ContentType', 'vector')
% 
% % Saving figures inside correct folder
% fprintf('Saving LiftComparisonWithPoints.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile LiftComparisonWithPoints.pdf Output
% movefile LiftComparisonWithPoints.png Output

%% LOAD A CLASS OF FUNCTIONS USEFUL TO EVALUATE ALL THE REQUIRED DATA
% obj2 = ShearBendingTorsion; 
% % Test speed
% % To check the work, a test case issue is provided. OpenVSP is called with
% % the following value of the airspeed. Also, the chosen flight condition is
% % relative to a lift coefficient equal to one (CL = 1). All the geometrical
% % parameters are inserted inside the OpenVSP input file. 
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Test_speed.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VA.value; 
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Test_speed.Attributes.unit = "m/s";
% 
% %% DEFINE A CHORD DISTRIBUTION 
% % Lift coefficient distribution at a global CL equal to one
% %
% % cl_at_CL1 = cl_unit_lift(CL1, CL2, cl1, cl2)
% % A simple function to evaluate the lift coefficient distribution along the
% % span cl = cl(y) when the associated global lift coefficient of the whole
% % wing is equal to 1.0; the function use a method similar to that suggested
% % by Abbott in Theory of Wing Section. See the complete documentation
% % inside the cl_unit_lift.m file
% % CL1 = trapz(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.y_half_span.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value(3,:)')/Aircraft.Geometry.Wing.S.value;
% % CL2 = trapz(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.y_half_span.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value(4,:)')/Aircraft.Geometry.Wing.S.value;
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value = cl_unit_lift(CL1, ...
% %                                                                                          CL2, ...
% %                                                                                          Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value(3,:)', ...
% %                                                                                          Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value(4,:)');
% 
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value = cl_unit_lift(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CL.value(3), ...
% %                                                                                          Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CL.value(4), ...
% %                                                                                          Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value(3,:)', ...
% %                                                                                          Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value(4,:)');
% CL_equal_to_one = 0.0;
% global_CL = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value(1,:)), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value(1,:))
%     global_CL(i) = trapz(half_span, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value(i,:))/Aircraft.Geometry.Wing.S.value;
%     if (global_CL(i) >= 1.0-1e-2) && (global_CL(i) <= 1.0+1e-2)
%             CL_equal_to_one = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value(i,:);
%     end
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value = CL_equal_to_one;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.Attributes.unit = "Non dimensional";  
% 
% % A vector which contains all the stations along the main wing span.
% Aircraft.Geometry.Wing.y.value = linspace(0, ...
%                                           Aircraft.Geometry.Wing.b.value*0.5, ...
%                                           length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value));
% Aircraft.Geometry.Wing.y.Attributes.unit = "m";
%                                     
% % Main wing taper ratio
% Aircraft.Geometry.Wing.taper.value = Aircraft.Geometry.Wing.ctip.value/Aircraft.Geometry.Wing.croot.value;
% Aircraft.Geometry.Wing.taper.Attributes.unit = "Non dimensional";
% 
% % Calculation of a chord distribution with a convenient, simple function.
% % 
% % c(y) = calc_chord(Swing, taper_ratio, span, y)
% % A complete documentation of this function is included inside the class
% % ShearBendingTorsion.m
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value = calc_chord(obj2, Aircraft.Geometry.Wing.S.value, ...
%                                                                                          Aircraft.Geometry.Wing.taper.value, ...
%                                                                                          Aircraft.Geometry.Wing.b.value, ...
%                                                                                          Aircraft.Geometry.Wing.y.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.Attributes.unit = "m";
% 
% %% POINT S CALCULATIONS                 
% % Lift coefficient distribution along the span at the Point S
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cl_S.value = abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cl_S.Attributes.unit = "Non dimensional";
% 
% % Drag coefficient ditribution along the span at the Point S (close to
% % stall)
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cd_S.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:)';
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cd_S.Attributes.unit = "Non dimensional";
% 
% % Pitching moment coefficient distribution along the span at Point S
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cm_S.value =  Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value(7,:)';
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cm_S.Attributes.unit = "Non dimensional";
% 
% %% PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
% % In this section of the code two vectors are defined to store the product 
% % c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCl_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cl_S.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCl_distr.Attributes.unit = "m";
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCd_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cd_S.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCd_distr.Attributes.unit = "m";
% 
% %% ALPHA VECTOR ALONG THE MAIN WING SPAN 
% cd .. 
% cd utilities
% % The 'dir' variable contains working directory path saved as a
% % char value
% dir = pwd;
% % Store working directory inside the log file
% fprintf('-----------------');
% fprintf('\n');
% fprintf('### Current directory ###');
% fprintf('\n');
% fprintf('%s\n', dir);
% 
% % Calling the function alpha_calc
% % Angle of attack = ALPHA(obj, CL, CL0, CL_star, CLalpha, a, b, c)
% % This function is included in the class file aero_model.m; a complete
% % description of this function is included inside that file.
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_alongspan_deg.value = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cl_S.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% %                          alpha_calc_lin(obj1, ...
% %                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cl_S.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_alongspan_deg.Attributes.unit = "degrees"; 
% 
% % Convert to radians 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_alongspan_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_alongspan_deg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_alongspan_rad.Attributes.unit = "radians";
% 
% % Changing folder 
% cd .. 
% cd csvla
% % The 'dir' variable contains working directory path saved as a
% % char value
% dir = pwd;
% % Store working directory inside the log file
% fprintf('-----------------');
% fprintf('\n');
% fprintf('### Current directory ###');
% fprintf('\n');
% fprintf('%s\n', dir);
%% TOTAL ANGLE OF ATTACK DATA STORAGE
% 
% % AoA_Tot = AoA + Twist_angle of the main wing
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_deg.value = Aircraft.Geometry.Wing.twist_angle.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alpha_S.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_deg.Attributes.unit = "Degrees"; 
% 
% % Convert in radians 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_deg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_rad.Attributes.unit = "Radians";
% 
% %% PROJECTION OF FORCES ALONG WING AXES - NORMAL AND AXIAL FORCES 
% 
%% Calculation of the normal force coefficient
% % N = calc_normal_force(AoA_Tot, cCl, cCd)
% % This function will be used to evaluate the normal force coefficients
% % distribution along the span; it is possible to fin a complete
% % documentation inside the class file ShearBendingTorsion.m 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCz.value = calc_normal_force(obj2, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_deg.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCl_distr.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCd_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCz.Attributes.unit = "m"; 
% 
% % Calculation of the axial force coefficient 
% % A = calc_normal_force(AoA_Tot, cCl, cCd)
% % This function will be used to evaluate the axial force coefficients
% % distribution along the span; it is possible to fin a complete
% % documentation inside the class file ShearBendingTorsion.m 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCa.value = calc_axial_force(obj2, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_deg.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCl_distr.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCd_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCa.Attributes.unit = "m"; 
% 
% % Normal force 
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Normal_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCz.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Normal_force.Attributes.unit = "N/m";
% 
% % Axial force 
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Axial_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCa.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Axial_force.Attributes.unit = "N/m";
% 
% %% SHEAR FORCE CALCULATION 
% % A = calc_shear_force(AoA_Tot, y, cCZ)
% % A complete description of this function is available inside the class
% % file ShearBendingTorsion.m 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_distr.value = calc_shear_force(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Normal_force.value)*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_distr.Attributes.unit = "daN";
% 
% %% BENDING MOMENT CALCULATION 
% % BM = calc_bend_mom(y, S)
% % A complete description of this function is included inside the class file
% % ShearBendingTorsion.m
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Bend_mom_distr.value = calc_bend_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Bend_mom_distr.Attributes.unit = "daN*m";
% 
% %% PLANS FOR THE STRUCTURAL DIMENSIONING
% % To correctly size the aerostructures of the main lifting surface 
% % it is necessary to apply the procedure just developed to the 
% % critical points coming from the V-N diagram. Those point represents 
% % the most demanding flight conditions that our aircraft could survive. 
% % Those points are all stored inside: 
% % --> Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope
% % Retrieve the values and apply formulas to them. 
% 
% % Pitching moment per unit length
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.m_distr.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cm_S.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cm_S.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.m_distr.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cm_S.value(i)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.value*((Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value(i))^2);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.m_distr.Attributes.unit = "N";
% 
% % Torque applied
% % T = calc_tors_mom(obj, y, m)
% % A complete distribution of this function is included inside the class
% % file ShearBendingTorsion.m
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Tors_mom_distr.value = calc_tors_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.m_distr.value)*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Tors_mom_distr.Attributes.unit = "daN*m";
% 
% disp(" ");
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_distr.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Bend_mom_distr.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Tors_mom_distr.value];
% disp(" ++++++++++ POINT S OF THE FINAL ENVELOPE - SHEAR, BENDING AND TORSION ++++++++++ ")
% format = ' %6.6f          %6.6f          %6.6f\n';
% label  = ' Shear [daN]           Bending [daN * m]             Torsion [daN * m]\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % Subplots with Shear, Bending moment and Torsion
% %  fig1 = Bending_Shear_diag(y, Shear, Bend_mom, Torsion, Point)
% % A complete description of this function is included inside the class file
% % ShearBendingTorsion.m
% 
% disp(" ")
% disp(" ++++ FIGURE 12 - POINT S SHEAR, BENDING, TORSION ++++ ");
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_BendMom_diagram.value = Shear_Bending_Torsion_diag(obj2, flip(Aircraft.Geometry.Wing.y.value), ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Bend_mom_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Tors_mom_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value);
% 
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointS.pdf', 'ContentType', 'vector')
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointS.png', 'ContentType', 'vector')
% 
% % Saving figures inside correct folder
% fprintf('Saving ShearBendingTorsionDiagramPointS.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile ShearBendingTorsionDiagramPointS.pdf Output
% movefile ShearBendingTorsionDiagramPointS.png Output
% 
% %% POINT A CALCULATIONS                 
% % Lift coefficient distribution along the span at the Point A
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cl_A.value = abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cl_A.Attributes.unit = "Non dimensional";
% 
% % Drag coefficient ditribution along the span at the Point A (close to
% % stall)
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cd_A.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:)';
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cd_A.Attributes.unit = "Non dimensional";
% 
% % Pitching moment coefficient distribution along the span at Point A
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cm_A.value =  Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value(7,:)';
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cm_A.Attributes.unit = "Non dimensional";
% 
% %% PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
% % In this section of the code two vectors are defined to store the product 
% % c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCl_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cl_A.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCl_distr.Attributes.unit = "m";
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCd_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cd_A.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCd_distr.Attributes.unit = "m";
% 
% %% ALPHA VECTOR ALONG THE MAIN WING SPAN 
% cd .. 
% cd utilities
% % The 'dir' variable contains working directory path saved as a
% % char value
% dir = pwd;
% % Store working directory inside the log file
% fprintf('-----------------');
% fprintf('\n');
% fprintf('### Current directory ###');
% fprintf('\n');
% fprintf('%s\n', dir);
% 
% % Calling the function alpha_calc
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_alongspan_deg.value = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cl_A.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% %                          alpha_calc_lin(obj1, ...
% %                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cl_A.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_alongspan_deg.Attributes.unit = "degrees"; 
% 
% % Convert to radians 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_alongspan_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_alongspan_deg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_alongspan_rad.Attributes.unit = "radians";
% 
% % Changing folder 
% cd .. 
% cd csvla
% % The 'dir' variable contains working directory path saved as a
% % char value
% dir = pwd;
% % Store working directory inside the log file
% fprintf('-----------------');
% fprintf('\n');
% fprintf('### Current directory ###');
% fprintf('\n');
% fprintf('%s\n', dir);
% 
% %% TOTAL ANGLE OF ATTACK DATA STORAGE
% 
% % AoA_Tot = AoA + Twist_angle of the main wing
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_deg.value = Aircraft.Geometry.Wing.twist_angle.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alpha_A.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_deg.Attributes.unit = "Degrees"; 
% 
% % Convert in radians 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_deg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_rad.Attributes.unit = "Radians";
% 
% %% PROJECTION OF FORCES ALONG WING AXES - NORMAL AND AXIAL FORCES 
% 
% % Calculation of the normal force coefficient
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCz.value = calc_normal_force(obj2, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_deg.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCl_distr.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCd_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCz.Attributes.unit = "m"; 
% 
% % Calculation of the axial force coefficient 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCa.value = calc_axial_force(obj2, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_deg.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCl_distr.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCd_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCa.Attributes.unit = "m"; 
% 
% % Normal force 
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Normal_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCz.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Normal_force.Attributes.unit = "N/m";
% 
% % Axial force 
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Axial_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCa.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Axial_force.Attributes.unit = "N/m";
% 
% %% SHEAR FORCE CALCULATION 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_distr.value = calc_shear_force(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Normal_force.value)*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_distr.Attributes.unit = "daN";
% 
% %% BENDING MOMENT CALCULATION 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Bend_mom_distr.value = calc_bend_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Bend_mom_distr.Attributes.unit = "daN*m";
% 
% %% SHEAR, BENDING MOMENT AND TORSION DIAGRAM
% 
% % Pitching moment per unit length
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.m_distr.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cm_A.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cm_A.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.m_distr.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cm_A.value(i)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value*((Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value(i))^2);
% end
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.m_distr.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cm_A.value.*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value.^2);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.m_distr.Attributes.unit = "N";
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Tors_mom_distr.value = calc_tors_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.m_distr.value)*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Tors_mom_distr.Attributes.unit = "daN*m";
% 
% disp(" ")
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_distr.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Bend_mom_distr.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Tors_mom_distr.value];
% disp(" ++++++++++ POINT A OF THE FINAL ENVELOPE - SHEAR, BENDING AND TORSION ++++++++++ ")
% format = ' %6.6f          %6.6f          %6.6f\n';
% label  = ' Shear [daN]           Bending [daN * m]             Torsion [daN * m]\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % Subplots with Shear, Bending moment and Torsion
% 
% disp(" ")
% disp(" ++++ FIGURE 13 - POINT A SHEAR, BENDING, TORSION");
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_BendMom_diagram.value = Shear_Bending_Torsion_diag(obj2, flip(Aircraft.Geometry.Wing.y.value), ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Bend_mom_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Tors_mom_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.point_name.value);
% 
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointA.pdf', 'ContentType', 'vector')
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointA.png', 'ContentType', 'vector')
% % Saving figures inside correct folder
% fprintf('Saving ShearBendingTorsionDiagramPointA.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile ShearBendingTorsionDiagramPointA.pdf Output
% movefile ShearBendingTorsionDiagramPointA.png Output
% 
% %% POINT C CALCULATIONS                 
% % Lift coefficient distribution along the span at the Point C
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cl_C.value = abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cl_C.Attributes.unit = "Non dimensional";
% 
% % Drag coefficient ditribution along the span at the Point C (close to
% % stall)
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cd_C.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:)';
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cd_C.Attributes.unit = "Non dimensional";
% 
% % Pitching moment coefficient distribution along the span at Point C
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.value =  Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value(7,:)';
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.Attributes.unit = "Non dimensional";
% 
% %% IMPORT DATA FROM THE CORRECT SOURCE
% 
% %% PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
% % In this section of the code two vectors are defined to store the product 
% % c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCl_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cl_C.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCl_distr.Attributes.unit = "m";
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCd_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cd_C.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCd_distr.Attributes.unit = "m";
% 
% %% ALPHA VECTOR ALONG THE MAIN WING SPAN 
% cd .. 
% cd utilities
% % The 'dir' variable contains working directory path saved as a
% % char value
% dir = pwd;
% % Store working directory inside the log file
% fprintf('-----------------');
% fprintf('\n');
% fprintf('### Current directory ###');
% fprintf('\n');
% fprintf('%s\n', dir);
% 
% % Calling the function alpha_calc
% % Angle of attack = ALPHA(obj, CL, CL0, CL_star, CLalpha, a, b, c)
% % 
% %   This function is able to evaluate the angle of attack of
% %   the aircraft in the prescribed flight condition. It must be
% %   noticed that. Search the complete description inside the class
% %   file aero_model.m
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_alongspan_deg.value = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cl_C.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% %                          alpha_calc_lin(obj1, ...
% %                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cl_C.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_alongspan_deg.Attributes.unit = "degrees"; 
% 
% % Convert to radians 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_alongspan_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_alongspan_deg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_alongspan_rad.Attributes.unit = "radians";
% 
% % Changing folder 
% cd .. 
% cd csvla
% % The 'dir' variable contains working directory path saved as a
% % char value
% dir = pwd;
% % Store working directory inside the log file
% fprintf('-----------------');
% fprintf('\n');
% fprintf('### Current directory ###');
% fprintf('\n');
% fprintf('%s\n', dir);
% 
% %% TOTAL ANGLE OF ATTACK DATA STORAGE
% 
% % AoA_Tot = AoA + Twist_angle of the main wing
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_deg.value = Aircraft.Geometry.Wing.twist_angle.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alpha_C.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_deg.Attributes.unit = "Degrees"; 
% 
% % Convert in radians 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_deg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_rad.Attributes.unit = "Radians";
% 
% %% PROJECTION OF FORCES ALONG WING AXES - NORMAL AND AXIAL FORCES 
% 
% % Calculation of the normal force coefficient
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCz.value = calc_normal_force(obj2, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_deg.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCl_distr.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCd_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCz.Attributes.unit = "m"; 
% 
% % Calculation of the axial force coefficient 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCa.value = calc_axial_force(obj2, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_deg.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCl_distr.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCd_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCa.Attributes.unit = "m"; 
% 
% % Normal force 
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Normal_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCz.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Normal_force.Attributes.unit = "N/m";
% 
% % Axial force 
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Axial_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCa.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Axial_force.Attributes.unit = "N/m";
% 
% %% SHEAR FORCE CALCULATION 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_distr.value = calc_shear_force(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Normal_force.value)*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_distr.Attributes.unit = "daN";
% 
% %% BENDING MOMENT CALCULATION 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Bend_mom_distr.value = calc_bend_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Bend_mom_distr.Attributes.unit = "daN*m";
% 
% %% SHEAR, BENDING MOMENT AND TORSION DIAGRAM
% 
% % Pitching moment per unit length
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.m_distr.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.m_distr.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.value(i)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value*((Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value(i))^2);
% end
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.m_distr.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.value.*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value.^2);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.m_distr.Attributes.unit = "N";
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Tors_mom_distr.value = calc_tors_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.m_distr.value)*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Tors_mom_distr.Attributes.unit = "daN*m";
% 
% disp(" ")
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_distr.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Bend_mom_distr.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Tors_mom_distr.value];
% disp(" ++++++++++ POINT C OF THE FINAL ENVELOPE - SHEAR, BENDING AND TORSION ++++++++++ ")
% format = ' %6.6f          %6.6f          %6.6f\n';
% label  = ' Shear [daN]           Bending [daN * m]             Torsion [daN * m]\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% disp(" ")
% % Subplots with shear, bending moment and torsion diagram
% disp(" ++++ FIGURE 14 - POINT C SHEAR, BENDING, TORSION");
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_BendMom_diagram.value = Shear_Bending_Torsion_diag(obj2, flip(Aircraft.Geometry.Wing.y.value), ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Bend_mom_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Tors_mom_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value);
% 
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointC.pdf', 'ContentType', 'vector')
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointC.png', 'ContentType', 'vector')
% % Saving figures inside correct folder
% fprintf('Saving ShearBendingTorsionDiagramPointC.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile ShearBendingTorsionDiagramPointC.pdf Output
% movefile ShearBendingTorsionDiagramPointC.png Output
% 
% %% POINT D CALCULATIONS                 
% % Lift coefficient distribution along the span at the Point D
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cl_D.value = abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cl_D.Attributes.unit = "Non dimensional";
% 
% % Interpolation to obtain Cd at the desired global CL
% x   = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:));
% y   = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(:,1));
% xi  = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:));
% yi  = 1:0.1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(:,1));
% [XX,YY] = meshgrid(x,y);
% [XI,YI] = meshgrid(xi,yi);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cd.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value, ...
%                                                                         XI, YI, 'spline');
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cd.Attributes.unit = 'Non dimensional';
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_Interpolation_Graph.value = cd_interpolation_graph(x, y, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value, ... 
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cd.value, XI, YI);            
%             
% % Pitching moment interpolation along the span
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cm.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value, ...
%                                                                         XI, YI, 'spline'); 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cm.Attributes.unit = 'Non dimensional';
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CM_Interpolation_Graph.value = cm_interpolation_graph(x, y, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value, ... 
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cm.value, XI, YI);
%             
% % Export the interpolation carpet plot
% % Cd interpolation and saving diagram
% 
% disp(" ")
% disp(" ++++ FIGURE 15 - DRAG SPANWISE DISTRIBUTION INTERPOLATION ++++ ");
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_Interpolation_Graph.value, 'CdInterpolationdiag.pdf', 'ContentType', 'vector')
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_Interpolation_Graph.value, 'CdInterpolationdiag.png', 'ContentType', 'vector')
% % Saving figures inside correct folder
% fprintf('Saving CdInterpolationdiag.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile CdInterpolationdiag.pdf Output
% movefile CdInterpolationdiag.png Output
% 
% % Cd interpolation and saving diagram
% 
% disp(" ++++ FIGURE 16 - PITCH MOM. SPANWISE DISTRIBUTION INTERPOLATION ++++ ");
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CM_Interpolation_Graph.value, 'CmInterpolationdiag.pdf', 'ContentType', 'vector')
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CM_Interpolation_Graph.value, 'CmInterpolationdiag.png', 'ContentType', 'vector')
% % Saving figures inside correct folder
% fprintf('Saving CmInterpolationdiag.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile CmInterpolationdiag.pdf Output
% movefile CmInterpolationdiag.png Output
%             
% % Integrate along semi-span and assign drag and pitching moment along span
% Aircraft.Geometry.Wing.half_span_y.value = linspace(0, ...
%                                           Aircraft.Geometry.Wing.b.value*0.5, ...
%                                           length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value));
% Aircraft.Geometry.Wing.half_span_y.Attributes.unit = 'm';
% 
% % Selection of the interpolated distribution of CD and CM
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Global_CD.value = zeros(length(yi), 1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Global_CD.Attributes.unit = "Non dimensional";
% for i = 1:length(yi)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Global_CD.value(i) = trapz(Aircraft.Geometry.Wing.half_span_y.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cd.value(i,:));
%     if abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Global_CD.value(i) - Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_D.value) < 1e-2
%        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cd_D.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cd.value(i,:)';
%        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cm.value(i,:)';
%     end
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cd_D.Attributes.unit = "Non dimensional";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.Attributes.unit = "Non dimensional";
% 
% %% IMPORT DATA FROM THE CORRECT SOURCE
% 
% %% PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
% % In this section of the code two vectors are defined to store the product 
% % c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCl_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cl_D.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCl_distr.Attributes.unit = "m";
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCd_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cd_D.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCd_distr.Attributes.unit = "m";
% 
% %% ALPHA VECTOR ALONG THE MAIN WING SPAN 
% cd .. 
% cd utilities
% % The 'dir' variable contains working directory path saved as a
% % char value
% dir = pwd;
% % Store working directory inside the log file
% fprintf('-----------------');
% fprintf('\n');
% fprintf('### Current directory ###');
% fprintf('\n');
% fprintf('%s\n', dir);
% 
% % Calling the function alpha_calc
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_alongspan_deg.value = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cl_D.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% %                          alpha_calc_lin(obj1, ...
% %                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cl_D.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_alongspan_deg.Attributes.unit = "degrees"; 
% 
% % Convert to radians 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_alongspan_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_alongspan_deg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_alongspan_rad.Attributes.unit = "radians";
% 
% % Changing folder 
% cd .. 
% cd csvla
% % The 'dir' variable contains working directory path saved as a
% % char value
% dir = pwd;
% % Store working directory inside the log file
% fprintf('-----------------');
% fprintf('\n');
% fprintf('### Current directory ###');
% fprintf('\n');
% fprintf('%s\n', dir);
% 
% %% TOTAL ANGLE OF ATTACK DATA STORAGE
% 
% % AoA_Tot = AoA + Twist_angle of the main wing
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_deg.value = Aircraft.Geometry.Wing.twist_angle.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alpha_D.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_deg.Attributes.unit = "Degrees"; 
% 
% % Convert in radians 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_deg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_rad.Attributes.unit = "Radians";
% 
% %% PROJECTION OF FORCES ALONG WING AXES - NORMAL AND AXIAL FORCES 
% 
% % Calculation of the normal force coefficient
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCz.value = calc_normal_force(obj2, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_deg.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCl_distr.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCd_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCz.Attributes.unit = "m"; 
% 
% % Calculation of the axial force coefficient 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCa.value = calc_axial_force(obj2, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_deg.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCl_distr.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCd_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCa.Attributes.unit = "m"; 
% 
% % Normal force 
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Normal_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCz.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Normal_force.Attributes.unit = "N/m";
% 
% % Axial force 
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Axial_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCa.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Axial_force.Attributes.unit = "N/m";
% 
% %% SHEAR FORCE CALCULATION 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_distr.value = calc_shear_force(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Normal_force.value)*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_distr.Attributes.unit = "daN";
% 
% %% BENDING MOMENT CALCULATION 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Bend_mom_distr.value = calc_bend_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Bend_mom_distr.Attributes.unit = "daN*m";
% 
% %% SHEAR, BENDING MOMENT AND TORSION DIAGRAM
% 
% % Pitching moment per unit length
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.m_distr.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.m_distr.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.value(i)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value*((Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value(i))^2);
% end
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.m_distr.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.value.*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value.^2);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.m_distr.Attributes.unit = "N";
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Tors_mom_distr.value = calc_tors_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.m_distr.value)*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Tors_mom_distr.Attributes.unit = "daN*m";
% 
% disp(" ")
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_distr.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Bend_mom_distr.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Tors_mom_distr.value];
% disp(" ++++++++++ POINT D OF THE FINAL ENVELOPE - SHEAR, BENDING AND TORSION ++++++++++ ")
% format = ' %6.6f          %6.6f          %6.6f\n';
% label  = ' Shear [daN]           Bending [daN * m]             Torsion [daN * m]\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % Subplots with shear, bending moment and torsion 
% 
% disp(" ")
% disp(" ++++ FIGURE 17 - POINT D SHEAR, BENDING, TORSION ++++ ");
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_BendMom_diagram.value = Shear_Bending_Torsion_diag(obj2, flip(Aircraft.Geometry.Wing.y.value), ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Bend_mom_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Tors_mom_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value);
% 
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointD.pdf', 'ContentType', 'vector')
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointD.png', 'ContentType', 'vector')
% 
% % Saving figures inside correct folder
% fprintf('Saving ShearBendingTorsionDiagramPointD.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile ShearBendingTorsionDiagramPointD.pdf Output
% movefile ShearBendingTorsionDiagramPointD.png Output
% 
% %% POINT F CALCULATIONS                 
% % Lift coefficient distribution along the span at the Point F
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cl_F.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cl_F.Attributes.unit = "Non dimensional";
% 
% % Interpolation to obtain Cd at the desired global CL
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Cd.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value, ...
%                                                                         XI, YI, 'spline');
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Cd.Attributes.unit = 'Non dimensional';
% 
% % Pitching moment interpolation along the span
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Cm.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value, ...
%                                                                         XI, YI, 'spline'); 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Cm.Attributes.unit = 'Non dimensional';
%           
% % Integrate along semi-span and assign drag and pitching moment along span
% Aircraft.Geometry.Wing.half_span_y.value = linspace(0, ...
%                                           Aircraft.Geometry.Wing.b.value*0.5, ...
%                                           length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value));
% Aircraft.Geometry.Wing.half_span_y.Attributes.unit = 'm';
% 
% % Selection of CD and CM distribution for Point F
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Global_CD.value = zeros(length(yi), 1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Global_CD.Attributes.unit = "Non dimensional";
% for i = 1:length(yi)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Global_CD.value(i) = trapz(Aircraft.Geometry.Wing.half_span_y.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Cd.value(i,:));
%     if abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Global_CD.value(i) - Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.value) < 1e-2
%        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cd_F.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Cd.value(i,:)';
%        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cm_F.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Cm.value(i,:)';
%     end
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cd_F.Attributes.unit = "Non dimensional";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cm_F.Attributes.unit = "Non dimensional";
% 
% %% IMPORT DATA FROM THE CORRECT SOURCE
% 
% %% PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
% % In this section of the code two vectors are defined to store the product 
% % c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCl_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cl_F.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCl_distr.Attributes.unit = "m";
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCd_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cd_F.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCd_distr.Attributes.unit = "m";
% 
% %% ALPHA VECTOR ALONG THE MAIN WING SPAN 
% cd .. 
% cd utilities
% % The 'dir' variable contains working directory path saved as a
% % char value
% dir = pwd;
% % Store working directory inside the log file
% fprintf('-----------------');
% fprintf('\n');
% fprintf('### Current directory ###');
% fprintf('\n');
% fprintf('%s\n', dir);
% 
% % Calling the function alpha_calc
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_alongspan_deg.value = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cl_F.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p) - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_alongspan_deg.Attributes.unit = "degrees"; 
% 
% % Convert to radians 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_alongspan_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_alongspan_deg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_alongspan_rad.Attributes.unit = "radians";
% 
% % Changing folder 
% cd .. 
% cd csvla
% % The 'dir' variable contains working directory path saved as a
% % char value
% dir = pwd;
% % Store working directory inside the log file
% fprintf('-----------------');
% fprintf('\n');
% fprintf('### Current directory ###');
% fprintf('\n');
% fprintf('%s\n', dir);
% 
% %% TOTAL ANGLE OF ATTACK DATA STORAGE
% 
% % AoA_Tot = AoA + Twist_angle of the main wing
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_Tot_deg.value = Aircraft.Geometry.Wing.twist_angle.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alpha_F.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_Tot_deg.Attributes.unit = "Degrees"; 
% 
% % Convert in radians 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_Tot_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_Tot_deg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_Tot_rad.Attributes.unit = "Radians";
% 
% %% PROJECTION OF FORCES ALONG WING AXES - NORMAL AND AXIAL FORCES 
% 
% % Calculation of the normal force coefficient
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCz.value = calc_normal_force(obj2, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_Tot_deg.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCl_distr.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCd_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCz.Attributes.unit = "m"; 
% 
% % Calculation of the axial force coefficient 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCa.value = calc_axial_force(obj2, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_Tot_deg.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCl_distr.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCd_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCa.Attributes.unit = "m"; 
% 
% % Normal force 
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Normal_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCz.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Normal_force.Attributes.unit = "N/m";
% 
% % Axial force 
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Axial_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCa.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Axial_force.Attributes.unit = "N/m";
% 
% %% SHEAR FORCE CALCULATION 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Shear_distr.value = calc_shear_force(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Normal_force.value)*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Shear_distr.Attributes.unit = "daN";
% 
% %% BENDING MOMENT CALCULATION 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Bend_mom_distr.value = calc_bend_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Shear_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Bend_mom_distr.Attributes.unit = "daN*m";
% 
% %% SHEAR, BENDING MOMENT AND TORSION DIAGRAM
% 
% % Pitching moment per unit length
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.m_distr.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cm_F.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cm_F.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.m_distr.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cm_F.value(i)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value*((Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value(i))^2);
% end
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.m_distr.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cm_F.value.*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value.^2);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.m_distr.Attributes.unit = "N";
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Tors_mom_distr.value = calc_tors_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.m_distr.value)*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Tors_mom_distr.Attributes.unit = "daN*m";
% 
% disp(" ")
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Shear_distr.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Bend_mom_distr.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Tors_mom_distr.value];
% disp(" ++++++++++ POINT F OF THE FINAL ENVELOPE - SHEAR, BENDING AND TORSION ++++++++++ ")
% format = ' %6.6f          %6.6f          %6.6f\n';
% label  = ' Shear [daN]           Bending [daN * m]             Torsion [daN * m]\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % Subplots with shear, bending moment and torsion
% 
% disp(" ")
% disp(" ++++ FIGURE 18 - POINT F SHEAR, BENDING, TORSION ++++ ");
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Shear_BendMom_diagram.value = Shear_Bending_Torsion_diag(obj2, flip(Aircraft.Geometry.Wing.y.value), ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Shear_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Bend_mom_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Tors_mom_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value);
% 
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointF.pdf', 'ContentType', 'vector')
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointF.png', 'ContentType', 'vector')
% 
% % Saving figures inside correct folder
% fprintf('Saving ShearBendingTorsionDiagramPointF.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile ShearBendingTorsionDiagramPointF.pdf Output
% movefile ShearBendingTorsionDiagramPointF.png Output
% 
% %% POINT G CALCULATIONS                 
% % Lift coefficient distribution along the span at the Point G
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cl_G.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cl_G.Attributes.unit = "Non dimensional";
% 
% % Interpolation to obtain Cd at the desired global CL
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Cd.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value, ...
%                                                                         XI, YI, 'spline');
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Cd.Attributes.unit = 'Non dimensional';
% 
% % Pitching moment interpolation along the span
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Cm.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value, ...
%                                                                         XI, YI, 'spline'); 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Cm.Attributes.unit = 'Non dimensional';
%            
% % Integrate along semi-span and assign drag and pitching moment along span
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Global_CD.value = zeros(length(yi), 1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Global_CD.Attributes.unit = "Non dimensional";
% for i = 1:length(yi)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Global_CD.value(i) = trapz(Aircraft.Geometry.Wing.half_span_y.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Cd.value(i,:));
%     if abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Global_CD.value(i) - Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.value) < 1e-2
%        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cd_G.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Cd.value(i,:)';
%        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cm_G.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Cm.value(i,:)';
%     end
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cd_G.Attributes.unit = "Non dimensional";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cm_G.Attributes.unit = "Non dimensional";
% 
% %% IMPORT DATA FROM THE CORRECT SOURCE
% 
% %% PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
% % In this section of the code two vectors are defined to store the product 
% % c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCl_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cl_G.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCl_distr.Attributes.unit = "m";
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCd_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cd_G.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCd_distr.Attributes.unit = "m";
% 
% %% ALPHA VECTOR ALONG THE MAIN WING SPAN 
% cd .. 
% cd utilities
% % The 'dir' variable contains working directory path saved as a
% % char value
% dir = pwd;
% % Store working directory inside the log file
% fprintf('-----------------');
% fprintf('\n');
% fprintf('### Current directory ###');
% fprintf('\n');
% fprintf('%s\n', dir);
% 
% % Calling the function alpha_calc
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_alongspan_deg.value = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cl_G.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p) - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_alongspan_deg.Attributes.unit = "degrees"; 
% 
% % Convert to radians 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_alongspan_rad.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_alongspan_deg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_alongspan_rad.Attributes.unit = "Degrees";
% 
% % Changing folder 
% cd .. 
% cd csvla
% % The 'dir' variable contains working directory path saved as a
% % char value
% dir = pwd;
% % Store working directory inside the log file
% fprintf('-----------------');
% fprintf('\n');
% fprintf('### Current directory ###');
% fprintf('\n');
% fprintf('%s\n', dir);
% 
% %% TOTAL ANGLE OF ATTACK DATA STORAGE
% 
% % AoA_Tot = AoA + Twist_angle of the main wing
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_Tot_deg.value = Aircraft.Geometry.Wing.twist_angle.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alpha_G.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_Tot_deg.Attributes.unit = "Degrees"; 
% 
% % Convert in radians 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_Tot_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_Tot_deg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_Tot_rad.Attributes.unit = "Radians";
% 
% %% PROJECTION OF FORCES ALONG WING AXES - NORMAL AND AXIAL FORCES 
% 
% % Calculation of the normal force coefficient
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCz.value = calc_normal_force(obj2, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_Tot_deg.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCl_distr.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCd_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCz.Attributes.unit = "m"; 
% 
% % Calculation of the axial force coefficient 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCa.value = calc_axial_force(obj2, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_Tot_deg.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCl_distr.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCd_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCa.Attributes.unit = "m"; 
% 
% % Normal force 
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Normal_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCz.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Normal_force.Attributes.unit = "N/m";
% 
% % Axial force 
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Axial_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCa.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Axial_force.Attributes.unit = "N/m";
% 
% %% SHEAR FORCE CALCULATION 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Shear_distr.value = calc_shear_force(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Normal_force.value)*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Shear_distr.Attributes.unit = "daN";
% 
% %% BENDING MOMENT CALCULATION 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Bend_mom_distr.value = calc_bend_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Shear_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Bend_mom_distr.Attributes.unit = "daN*m";
% 
% %% SHEAR, BENDING MOMENT AND TORSION DIAGRAM
% 
% % Pitching moment per unit length
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.m_distr.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cm_G.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cm_G.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.m_distr.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cm_G.value(i)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value*((Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value(i))^2);
% end
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.m_distr.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cm_G.value.*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value.^2);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.m_distr.Attributes.unit = "N";
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Tors_mom_distr.value = calc_tors_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.m_distr.value)*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Tors_mom_distr.Attributes.unit = "daN*m";
% 
% disp(" ")
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Shear_distr.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Bend_mom_distr.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Tors_mom_distr.value];
% disp(" ++++++++++ POINT G OF THE FINAL ENVELOPE - SHEAR, BENDING AND TORSION ++++++++++ ")
% format = ' %6.6f          %6.6f          %6.6f\n';
% label  = ' Shear [daN]           Bending [daN * m]             Torsion [daN * m]\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % Subplots with shear, bending moment and torsion diagram
% 
% disp(" ")
% disp(" ++++ FIGURE 19 - POINT G SHEAR, BENDING, TORSION ++++ ");
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Shear_BendMom_diagram.value = Shear_Bending_Torsion_diag(obj2, flip(Aircraft.Geometry.Wing.y.value), ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Shear_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Bend_mom_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Tors_mom_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.point_name.value);
% 
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointG.pdf', 'ContentType', 'vector')
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointG.png', 'ContentType', 'vector')
% 
% % Saving figures inside correct folder
% fprintf('Saving ShearBendingTorsionDiagramPointG.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile ShearBendingTorsionDiagramPointG.pdf Output
% movefile ShearBendingTorsionDiagramPointG.png Output
% 
% %% POINT E CALCULATIONS                 
% % Lift coefficient distribution along the span at the Point E
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cl_E.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cl_E.Attributes.unit = "Non dimensional";
% 
% % Interpolation to obtain Cd at the desired global CL
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Cd.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value, ...
%                                                                         XI, YI, 'spline');
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Cd.Attributes.unit = 'Non dimensional';
% 
% % Pitching moment interpolation along the span
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Cm.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value, ...
%                                                                         XI, YI, 'spline'); 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Cm.Attributes.unit = 'Non dimensional';
%            
% % Integrate along semi-span and assign drag and pitching moment along span
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Global_CD.value = zeros(length(yi), 1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Global_CD.Attributes.unit = "Non dimensional";
% for i = 1:length(yi)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Global_CD.value(i) = trapz(Aircraft.Geometry.Wing.half_span_y.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Cd.value(i,:));
%     if abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Global_CD.value(i) - Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CD_E.value) < 1e-2
%        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cd_E.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Cd.value(i,:)';
%        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cm_E.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Cm.value(i,:)';
%     end
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cd_E.Attributes.unit = "Non dimensional";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cm_E.Attributes.unit = "Non dimensional";
% 
% %% IMPORT DATA FROM THE CORRECT SOURCE
% 
% %% PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
% % In this section of the code two vectors are defined to store the product 
% % c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCl_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cl_E.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCl_distr.Attributes.unit = "m";
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCd_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cd_E.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCd_distr.Attributes.unit = "m";
% 
% %% ALPHA VECTOR ALONG THE MAIN WING SPAN 
% cd .. 
% cd utilities
% % The 'dir' variable contains working directory path saved as a
% % char value
% dir = pwd;
% % Store working directory inside the log file
% fprintf('-----------------');
% fprintf('\n');
% fprintf('### Current directory ###');
% fprintf('\n');
% fprintf('%s\n', dir);
% 
% % Calling the function 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_alongspan_deg.value = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cl_E.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p) - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
% %                          alpha_calc_lin(obj1, ...
% %                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cl_E.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value) - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_alongspan_deg.Attributes.unit = "degrees"; 
% 
% % Convert to radians 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_alongspan_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_alongspan_deg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_alongspan_rad.Attributes.unit = "radians";
% 
% % Changing folder 
% cd .. 
% cd csvla
% % The 'dir' variable contains working directory path saved as a
% % char value
% dir = pwd;
% % Store working directory inside the log file
% fprintf('-----------------');
% fprintf('\n');
% fprintf('### Current directory ###');
% fprintf('\n');
% fprintf('%s\n', dir);
% 
% %% TOTAL ANGLE OF ATTACK DATA STORAGE
% 
% % AoA_Tot = AoA + Twist_angle of the main wing
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_Tot_deg.value = Aircraft.Geometry.Wing.twist_angle.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alpha_E.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_Tot_deg.Attributes.unit = "Degrees"; 
% 
% % Convert in radians 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_Tot_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_Tot_deg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_Tot_rad.Attributes.unit = "Radians";
% 
% %% PROJECTION OF FORCES ALONG WING AXES - NORMAL AND AXIAL FORCES 
% 
% % Calculation of the normal force coefficient
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCz.value = calc_normal_force(obj2, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_Tot_deg.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCl_distr.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCd_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCz.Attributes.unit = "m"; 
% 
% % Calculation of the axial force coefficient 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCa.value = calc_axial_force(obj2, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_Tot_deg.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCl_distr.value, ...
%                                                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCd_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCa.Attributes.unit = "m"; 
% 
% % Normal force 
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Normal_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCz.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Normal_force.Attributes.unit = "N/m";
% 
% % Axial force 
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Axial_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCa.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Axial_force.Attributes.unit = "N/m";
% 
% %% SHEAR FORCE CALCULATION 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Shear_distr.value = calc_shear_force(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Normal_force.value)*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Shear_distr.Attributes.unit = "daN";
% 
% %% BENDING MOMENT CALCULATION 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Bend_mom_distr.value = calc_bend_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Shear_distr.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Bend_mom_distr.Attributes.unit = "daN*m";
% 
% %% SHEAR, BENDING MOMENT AND TORSION DIAGRAM
% 
% % Pitching moment per unit length
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.m_distr.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cm_E.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cm_E.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.m_distr.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cm_E.value(i)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value*((Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value(i))^2);
% end
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.m_distr.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cm_E.value.*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value.^2);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.m_distr.Attributes.unit = "N";
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Tors_mom_distr.value = calc_tors_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
%                                                                                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.m_distr.value)*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Tors_mom_distr.Attributes.unit = "daN*m";
% 
% disp(" ")
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Shear_distr.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Bend_mom_distr.value, ...
%              Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Tors_mom_distr.value];
% 
% disp(" ")
% disp(" ++++++++++ POINT E OF THE FINAL ENVELOPE - SHEAR, BENDING AND TORSION ++++++++++ ")
% format = ' %6.6f          %6.6f          %6.6f\n';
% label  = ' Shear [daN]           Bending [daN * m]             Torsion [daN * m]\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % Subplots with shear, bending moment and torsion diagram
% 
% disp(" ")
% disp(" ++++ FIGURE 20 - POINT E SHEAR, BENDING, TORSION ++++ ");
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Shear_BendMom_diagram.value = Shear_Bending_Torsion_diag(obj2, flip(Aircraft.Geometry.Wing.y.value), ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Shear_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Bend_mom_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Tors_mom_distr.value, ...
%                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.point_name.value);
% 
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointE.pdf', 'ContentType', 'vector')
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointE.png', 'ContentType', 'vector')
% 
% % Saving figures inside correct folder
% fprintf('Saving ShearBendingTorsionDiagramPointE.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile ShearBendingTorsionDiagramPointE.pdf Output
% movefile ShearBendingTorsionDiagramPointE.png Output
% 
% %% COMPARING SHEAR CURVES 
% 
% disp(" ")
% disp(" ++++ FIGURE 21 - SHEAR COMPARISON ++++ ");
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Shear_Comparison.value = Compare_Shear_curves(flip(Aircraft.Geometry.Wing.y.value), ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Shear_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Shear_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Shear_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.point_name.value);
% % Saving
% 
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Shear_Comparison.value, 'ShearComparison.pdf', 'ContentType', 'vector')
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Shear_Comparison.value, 'ShearComparison.png', 'ContentType', 'vector')
% 
% % Saving figures inside correct folder
% fprintf('Saving ShearComparison.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile ShearComparison.pdf Output
% movefile ShearComparison.png Output
% 
% %% COMPARING BENDING CURVES 
% 
% disp(" ")
% disp(" ++++ FIGURE 22 - BENDING COMPARISON ++++ ");
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Bending_Comparison.value = Compare_Bending_curves(flip(Aircraft.Geometry.Wing.y.value), ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Bend_mom_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Bend_mom_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Bend_mom_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Bend_mom_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Bend_mom_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Bend_mom_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Bend_mom_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.point_name.value);
% % Saving
% 
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Bending_Comparison.value, 'BendingComparison.pdf', 'ContentType', 'vector')
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Bending_Comparison.value, 'BendingComparison.png', 'ContentType', 'vector')
% 
% % Saving figures inside correct folder
% fprintf('Saving BendingComparison.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile BendingComparison.pdf Output
% movefile BendingComparison.png Output
% 
% %% COMPARING TORSION CURVES 
% 
% disp(" ")
% disp(" ++++ FIGURE 23 - TORSION COMPARISON ++++ ");
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Torsion_Comparison.value = Compare_Torsion_curves(flip(Aircraft.Geometry.Wing.y.value), ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Tors_mom_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Tors_mom_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Tors_mom_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Tors_mom_distr.value, ...
% - Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Tors_mom_distr.value, ...
% - Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Tors_mom_distr.value, ...
% - Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Tors_mom_distr.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.point_name.value);
% % Saving
% 
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Torsion_Comparison.value, 'TorsionComparison.pdf', 'ContentType', 'vector')
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Torsion_Comparison.value, 'TorsionComparison.png', 'ContentType', 'vector')
% 
% % Saving figures inside correct folder
% fprintf('Saving TorsionComparison.pdf Output in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile TorsionComparison.pdf Output
% movefile TorsionComparison.png Output
% 
% %% COMPARING DISTRIBUTION ALONG THE SPAN OF CL, CD, CM
% % =========================================================================================================================================================
% 
% disp(" ")
% disp(" ++++ FIGURE 24 - SPANWISE LIFT COMPARISON ++++ ");
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.cl_Comparison.value = Compare_cl_curves(Aircraft.Geometry.Wing.y.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cl_S.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cl_A.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cl_C.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cl_D.value, ...
% abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cl_F.value), ...
% abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cl_G.value), ...
% abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cl_E.value), ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.point_name.value);
% % Saving
% 
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.cl_Comparison.value, 'clComparison.pdf', 'ContentType', 'vector')
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.cl_Comparison.value, 'clComparison.png', 'ContentType', 'vector')
% 
% % Saving figures inside correct folder
% fprintf('Saving clComparison.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile clComparison.pdf Output
% movefile clComparison.png Output
% % =========================================================================================================================================================
% 
% disp(" ")
% disp(" ++++ FIGURE 25 - SPANWISE DRAG COMPARISON ++++ ");
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.cd_Comparison.value = Compare_cd_curves(Aircraft.Geometry.Wing.y.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cd_S.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cd_A.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cd_C.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cd_D.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cd_F.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cd_G.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cd_E.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.point_name.value);
% % Saving
% 
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.cd_Comparison.value, 'cdComparison.pdf', 'ContentType', 'vector')
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.cd_Comparison.value, 'cdComparison.png', 'ContentType', 'vector')
% 
% % Saving figures inside correct folder
% fprintf('Saving cdComparison.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile cdComparison.pdf Output
% movefile cdComparison.png Output
% 
% % =========================================================================================================================================================
% 
% disp(" ")
% disp(" ++++ FIGURE 26 - SPANWISE PITCH MOMENT COMPARISON ++++ ");
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.cm_Comparison.value = Compare_cm_curves(Aircraft.Geometry.Wing.y.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cm_S.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cm_A.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cm_F.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cm_G.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cm_E.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.point_name.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.point_name.value);
% % Saving
% 
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.cm_Comparison.value, 'cmComparison.pdf', 'ContentType', 'vector')
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.cm_Comparison.value, 'cmComparison.png', 'ContentType', 'vector')
% % Saving figures inside correct folder
% fprintf('Saving cmComparison.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile cmComparison.pdf Output
% movefile cmComparison.png Output
% % =========================================================================================================================================================