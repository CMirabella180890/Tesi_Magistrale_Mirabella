% % ---------------------------------------------------------------------------
% % INPUT SOURCE VARIABLE 
% % ---------------------------------------------------------------------------
% InputSource = "From File"; 
% 
% %% Script to evaluate balancing horizontal tail loads
% %   DESCRIPTION
% %    In this script, all the methods relative to the calculation of
% %    aerodynamic coefficients, necessary to the correct estimation of the
% %    airloads acting on the horizontal empennage. 
% 
% %% RETURN INSIDE CSVLA
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
% %% CLASS INSTANTIATION 
% obj1 = aero_model; 
% % 
% %% LIFT CHARACTERISTIC
% % Before to start the balancing loads analysis, it is required to plot lift
% % coefficient values calculated with a non-linear formulation. Then, linear
% % and non-linear calculation are compared with raw data available. The
% % selected range is -2.0<CL<2.0. To find a complete documentation of the
% % function 'CL_Non_linear_model(...)' search inside aero_model.m file. 
% 
% % NUMBER OF ELEMENTS
% numb = 1e3;
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
% if InputSource == "From File"
%     % -------------------------------------------------------------------------
%     % ALFA FUNCTION
%     alfa_func = @(rho, S, V, WS, n, CLalfa, alfa_0lift) (2 / rho) * (1 / V.^2) * (1/CLalfa) * (WS) * n + alfa_0lift;
%     % =========================================================================
%     CL_wb    = str2num(Aircraft.Certification.Aerodynamic_data.CL.value);
%     CD_wb    = str2num(Aircraft.Certification.Aerodynamic_data.CD.value); 
%     CM_wb    = str2num(Aircraft.Certification.Aerodynamic_data.CM.value);
%     alpha_wb = str2num(Aircraft.Certification.Aerodynamic_data.alpha.value);
%     % =========================================================================
%     % OSWALDT EFFICIENCY FACTOR
%     e   = Aircraft.Certification.Aerodynamic_data.e.value;
%     % ASPECT RATIO 
%     AR  = Aircraft.Geometry.Wing.AR.value;
%     % ZERO-LIFT DRAG COEFFICIENT
%     CD0 = Aircraft.Certification.Aerodynamic_data.CD0.value;
%     % ENDING OF LINEAR PART OF LIFT CURVE
%     CL_star = Aircraft.Certification.Aerodynamic_data.CL_star.value;
%     % CLALFA IN DEG^-1
%     CLalfa  = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;
%     % ZERO LIFT COEFFICIENT
%     CL0         = Aircraft.Certification.Aerodynamic_data.CL0.value;
%     CL_alfa_deg = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;
%     alfa_0l     = -CL0 / CL_alfa_deg;
%     Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value           = alfa_0l;
%     Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.Attributes.unit = "deg";
%     % =========================================================================
%     % -------------------------------------------------------------------------
%     % CL INTERPOLATION
%     % -------------------------------------------------------------------------
%     alpha_interp  = linspace(alpha_wb(1), alpha_wb(end), numb)';
%     % -------------------------------------------------------------------------
%     % =========================================================================
%     % Alpha_star 
%     % =========================================================================
%     p_alfa_star   = polyfit(CL_wb(4:5), alpha_wb(4:5), 2);
%     alfa_star     = polyval(p_alfa_star, CL_star);
%     Aircraft.Certification.Aerodynamic_data.alpha_star.value = alfa_star;
%     Aircraft.Certification.Aerodynamic_data.Attributes.unit  = "deg";
%     % -------------------------------------------------------------------------
%     CL_wb_interp  = zeros(length(alpha_interp), 1);
%     n_pol1        = 1;
%     n_pol2        = 2;
%     p_CL_wb1      = polyfit([alpha_wb(1:3); alfa_star],   [CL_wb(1:3); CL_star],   n_pol1);
%     p_CL_wb2      = polyfit([alfa_star; alpha_wb(4:end)], [CL_star; CL_wb(4:end)], n_pol2);
%     Aircraft.Certification.Aerodynamic_data.p_CL_wb1.value = p_CL_wb1; 
%     Aircraft.Certification.Aerodynamic_data.p_CL_wb1.Attributes.unit = "1/deg";
%     Aircraft.Certification.Aerodynamic_data.p_CL_wb2.value = p_CL_wb2; 
%     Aircraft.Certification.Aerodynamic_data.p_CL_wb2.Attributes.unit = "1/deg";
%     CL_wb_interp1 = polyval(p_CL_wb1, alpha_interp);
%     CL_wb_interp2 = polyval(p_CL_wb2, alpha_interp);
%     for i = 1:length(alpha_interp)
%         if (CL_wb_interp1(i) < CL_star) && (CL_wb_interp1(i) ~= CL_wb_interp2(i)) 
%             CL_wb_interp(i) = CL_wb_interp1(i);
%         elseif CL_wb_interp1(i) > CL_star 
%             CL_wb_interp(i) = CL_wb_interp2(i);
%         end
%     end
%     % -------------------------------------------------------------------------
%     Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.value           = alpha_interp;
%     Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.Attributes.unit = "deg";
%     Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value                = CL_wb_interp;
%     Aircraft.Certification.Aerodynamic_data.CL_fullmodel.Attributes.unit      = "Non dimensional";
%     Aircraft.Certification.Aerodynamic_data.CL_max_fullmodel.value            = max(Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value);
%     Aircraft.Certification.Aerodynamic_data.CL_max_fullmodel.Attributes.unit  = "Non dimensional";
%     CL_max                                                                    = Aircraft.Certification.Aerodynamic_data.CL_max_fullmodel.value;
%     % -------------------------------------------------------------------------
%     disp(" ")
%     disp(" ++++ FIGURE 6 - LIFT CURVE INTERPOLATION ++++ ");
%     Aircraft.Certification.Aerodynamic_data.CL_fullmodel_diagram.value = Lift_fullmodel_curve(Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.value, ...
%         Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value, ...
%         CL_wb, ...
%         alpha_wb);    
%     
%     % SAVING FIGURES
%     exportgraphics(Aircraft.Certification.Aerodynamic_data.CL_fullmodel_diagram.value, 'FullLiftModelInterpolation.pdf', 'ContentType', 'vector');
%     exportgraphics(Aircraft.Certification.Aerodynamic_data.CL_fullmodel_diagram.value, 'FullLiftModelInterpolation.png', 'ContentType', 'vector');
%     
%     % Saving figures inside correct folder
%     fprintf('Saving FullLiftModelInterpolation.pdf in: ');
%     fprintf('\n'); 
%     fprintf('%s\n', SaveFolder);
%     % Moving file inside correct folder
%     movefile FullLiftModelInterpolation.pdf Output
%     movefile FullLiftModelInterpolation.png Output    
%     % -------------------------------------------------------------------------
%     % -------------------------------------------------------------------------
%     a1 = p_CL_wb1(1);
%     b1 = p_CL_wb1(2); 
%     Aircraft.Certification.Aerodynamic_data.CLlin_PolCoeff_a.value = a1;
%     Aircraft.Certification.Aerodynamic_data.CLlin_PolCoeff_b.value = b1;
%    
%     a2 = p_CL_wb2(1);
%     b2 = p_CL_wb2(2); 
%     c2 = p_CL_wb2(3);
%     Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a.value = a2;
%     Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b.value = b2;
%     Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c.value = c2;
%     % -------------------------------------------------------------------------
%     CL0          = p_CL_wb1(2);
%     CLalfa       = p_CL_wb1(1);
%     Aircraft.Certification.Aerodynamic_data.CL0.value = CL0;
%     Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value = CLalfa;
%     alfa_zero_lift = alpha_fullmodel(0, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
%     Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value           = alfa_zero_lift;
%     Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.Attributes.unit = "deg";
%     % -------------------------------------------------------------------------
%     % CD INTERPOLATION
%     % -------------------------------------------------------------------------
%     p_cd_star   = polyfit(alpha_wb(4:5), CD_wb(4:5), 2);
%     CD_star     = polyval(p_cd_star, alfa_star);
%     Aircraft.Certification.Aerodynamic_data.CD_star.value           = CD_star;
%     Aircraft.Certification.Aerodynamic_data.CD_star.Attributes.unit = "Non dimensional";
%     % -------------------------------------------------------------------------
%     n_pol1        = 2;
%     p_CD_wb1      = polyfit([alpha_wb(1:3); alfa_star; alpha_wb(4:end)], ...
%                             [CD_wb(1:3); CD_star; CD_wb(4:end)],   n_pol1);
%     CD_wb_interp  = polyval(p_CD_wb1, alpha_interp);
%     Aircraft.Certification.Aerodynamic_data.CD_PolCoeff_p.value           = p_CD_wb1;
%     Aircraft.Certification.Aerodynamic_data.CD_PolCoeff_p.Attributes.unit = "Non dimensional";
%     % -------------------------------------------------------------------------
%     CDalfa = figure(150);
%     hold on
%     grid on 
%     grid minor
%     plot(alpha_interp, CD_wb_interp, '-r', 'LineWidth', 1.5)
%     plot(alpha_wb, CD_wb, 'k.', 'MarkerSize', 10)
%     xlim 'padded' ;
%     ylim 'padded' ;
%     xlabel("Angle of attack - $\alpha$ $(deg)$", "Interpreter", "latex")
%     ylabel("Drag coefficient - $C_{D_{wb}}$", "Interpreter", "latex")
%     title("Drag model", "Interpreter", "latex")
%     legend({'Model','Data points'}, 'Interpreter', 'latex', 'Location', 'southeast')
%     
%     % SAVING FIGURES
%     exportgraphics(CDalfa, 'FullDragModelInterpolation.pdf', 'ContentType', 'vector');
%     exportgraphics(CDalfa, 'FullDragModelInterpolation.png', 'ContentType', 'vector');
%     
%     % Saving figures inside correct folder
%     fprintf('Saving FullDragModelInterpolation.pdf in: ');
%     fprintf('\n'); 
%     fprintf('%s\n', SaveFolder);
%     % Moving file inside correct folder
%     movefile FullDragModelInterpolation.pdf Output
%     movefile FullDragModelInterpolation.png Output        
%     % -------------------------------------------------------------------------
%     % CM INTERPOLATION
%     % -------------------------------------------------------------------------
%     n_pol     = 3;
%     p_CM_wb   = polyfit(alpha_wb, CM_wb, n_pol);
%     CM_interp = polyval(p_CM_wb, alpha_interp);
%     Aircraft.Certification.Aerodynamic_data.CM_fullmodel.value = CM_interp;
%     Aircraft.Certification.Aerodynamic_data.CM_fullmodel.Attributes.unit = "Non dimensional";
%     % -------------------------------------------------------------------------   
%     CMALFA = figure(151);
%     hold on
%     grid on 
%     grid minor
%     plot(alpha_interp, CM_interp, '-r', 'LineWidth', 1.5)
%     plot(alpha_wb, CM_wb, 'k.', 'MarkerSize', 10)
%     xlim 'padded' ;
%     ylim 'padded' ;
%     xlabel("Angle of attack - $\alpha$ $(deg)$", "Interpreter", "latex")
%     ylabel("Pitch. mom. coefficient - $C_{M_{wb}}$", "Interpreter", "latex")
%     title("Pitch. mom. model", "Interpreter", "latex")
%     legend({'Model','Data points'}, 'Interpreter', 'latex', 'Location', 'southeast')
%     
%     % SAVING FIGURES
%     exportgraphics(CDalfa, 'FullPitchMomModelInterpolation.pdf', 'ContentType', 'vector');
%     exportgraphics(CDalfa, 'FullPitchMomModelInterpolation.png', 'ContentType', 'vector');
%     
%     % Saving figures inside correct folder
%     fprintf('Saving FullPitchMomModelInterpolation.pdf in: ');
%     fprintf('\n'); 
%     fprintf('%s\n', SaveFolder);
%     % Moving file inside correct folder
%     movefile FullPitchMomModelInterpolation.pdf Output
%     movefile FullPitchMomModelInterpolation.png Output    
%     % -------------------------------------------------------------------------
%     clear CL_wb CD_wb CM_wb alpha_wb alpha_interp CM_interp CD_wb_interp
%     % -------------------------------------------------------------------------
%     elseif  InputSource == "From User"
%         disp("No data file. Insert data.")
% end




% ---------------------------------------------------------------------------
% INPUT SOURCE VARIABLE 
% ---------------------------------------------------------------------------
InputSource = "From File"; 

%% Script to evaluate balancing horizontal tail loads
%   DESCRIPTION
%    In this script, all the methods relative to the calculation of
%    aerodynamic coefficients, necessary to the correct estimation of the
%    airloads acting on the horizontal empennage. 

%% RETURN INSIDE CSVLA
cd .. 
cd csvla
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

%% CLASS INSTANTIATION 
obj1 = aero_model; 
% 
%% LIFT CHARACTERISTIC
% Before to start the balancing loads analysis, it is required to plot lift
% coefficient values calculated with a non-linear formulation. Then, linear
% and non-linear calculation are compared with raw data available. The
% selected range is -2.0<CL<2.0. To find a complete documentation of the
% function 'CL_Non_linear_model(...)' search inside aero_model.m file. 

% NUMBER OF ELEMENTS
numb = 1e3;
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if InputSource == "From File"
    % -------------------------------------------------------------------------
    % ALFA FUNCTION
    alfa_func = @(rho, S, V, WS, n, CLalfa, alfa_0lift) (2 / rho) * (1 / V.^2) * (1/CLalfa) * (WS) * n + alfa_0lift;
    % =========================================================================
    CL_wb    = str2num(Aircraft.Certification.Aerodynamic_data.CL.value);
    CD_wb    = str2num(Aircraft.Certification.Aerodynamic_data.CD.value); 
    CM_wb    = str2num(Aircraft.Certification.Aerodynamic_data.CM.value);
    alpha_wb = str2num(Aircraft.Certification.Aerodynamic_data.alpha.value);
%     CL_wb    = Aircraft.Certification.Aerodynamic_data.CL.value;
%     CD_wb    = Aircraft.Certification.Aerodynamic_data.CD.value; 
%     CM_wb    = Aircraft.Certification.Aerodynamic_data.CM.value;
%     alpha_wb = Aircraft.Certification.Aerodynamic_data.alpha.value;
    % =========================================================================
    length_alpha = length(alpha_wb)
    % =========================================================================
    aid_to_interp = figure(1400); 
    hold on; grid on; grid minor; 
    plot(alpha_wb, CL_wb, '.k', 'MarkerSize', 10)
    xlim 'padded' ;
    ylim 'padded' ;
    xlabel("Angle of attack - $\alpha$ $(deg)$", "Interpreter", "latex")
    ylabel("Lift coefficient - $C_{L_{wb}}$", "Interpreter", "latex")
    title("Aid to interpolation", "Interpreter", "latex")
    %legend({'Model','Data points'}, 'Interpreter', 'latex', 'Location', 'southeast')
    % =========================================================================
    
    % =========================================================================
    % OSWALDT EFFICIENCY FACTOR
    e   = Aircraft.Certification.Aerodynamic_data.e.value;
    % ASPECT RATIO 
    AR  = Aircraft.Geometry.Wing.AR.value;
    % ZERO-LIFT DRAG COEFFICIENT
    CD0 = Aircraft.Certification.Aerodynamic_data.CD0.value;
    % ENDING OF LINEAR PART OF LIFT CURVE
    CL_star = Aircraft.Certification.Aerodynamic_data.CL_star.value;
    % CLALFA IN DEG^-1
    CLalfa  = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;
    % ZERO LIFT COEFFICIENT
    CL0         = Aircraft.Certification.Aerodynamic_data.CL0.value;
    CL_alfa_deg = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;
    alfa_0l     = -CL0 / CL_alfa_deg;
    Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value           = alfa_0l;
    Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.Attributes.unit = "deg";
    % =========================================================================
    % -------------------------------------------------------------------------
    % CL INTERPOLATION
    % -------------------------------------------------------------------------
    alpha_interp  = linspace(alpha_wb(1), alpha_wb(end), numb)';
    % -------------------------------------------------------------------------
    % =========================================================================
    % Alpha_star 
    % =========================================================================
    % CL_difference = CL_wb-CL_star
    index_cl_star  = find((abs(CL_wb-CL_star))<1e-1);
    index_aoa_star = index_cl_star(1);
    p_alfa_star    = polyfit(CL_wb(index_aoa_star-1:index_aoa_star+1), alpha_wb(index_aoa_star-1:index_aoa_star+1), 2);
    alfa_star      = polyval(p_alfa_star, CL_star);
    Aircraft.Certification.Aerodynamic_data.alpha_star.value = alfa_star;
    Aircraft.Certification.Aerodynamic_data.Attributes.unit  = "deg";
    % -------------------------------------------------------------------------
    CL_wb_interp  = zeros(length(alpha_interp), 1);
    n_pol1        = 1;
    n_pol2        = 2;
    p_CL_wb1      = polyfit([alpha_wb(1:index_aoa_star-1); alfa_star],   [CL_wb(1:index_aoa_star-1); CL_star],   n_pol1);
    p_CL_wb2      = polyfit([alfa_star; alpha_wb(index_aoa_star:end)], [CL_star; CL_wb(index_aoa_star:end)], n_pol2);
    Aircraft.Certification.Aerodynamic_data.p_CL_wb1.value = p_CL_wb1; 
    Aircraft.Certification.Aerodynamic_data.p_CL_wb1.Attributes.unit = "1/deg";
    Aircraft.Certification.Aerodynamic_data.p_CL_wb2.value = p_CL_wb2; 
    Aircraft.Certification.Aerodynamic_data.p_CL_wb2.Attributes.unit = "1/deg";
    CL_wb_interp1 = polyval(p_CL_wb1, alpha_interp);
    CL_wb_interp2 = polyval(p_CL_wb2, alpha_interp);
    for i = 1:length(alpha_interp)
        if (CL_wb_interp1(i) < CL_star) && (CL_wb_interp1(i) ~= CL_wb_interp2(i)) 
            CL_wb_interp(i) = CL_wb_interp1(i);
        elseif CL_wb_interp1(i) > CL_star 
            CL_wb_interp(i) = CL_wb_interp2(i);
        end
    end
    % -------------------------------------------------------------------------
    Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.value           = alpha_interp;
    Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.Attributes.unit = "deg";
    Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value                = CL_wb_interp;
    Aircraft.Certification.Aerodynamic_data.CL_fullmodel.Attributes.unit      = "Non dimensional";
    Aircraft.Certification.Aerodynamic_data.CL_max_fullmodel.value            = max(Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value);
    Aircraft.Certification.Aerodynamic_data.CL_max_fullmodel.Attributes.unit  = "Non dimensional";
    CL_max                                                                    = Aircraft.Certification.Aerodynamic_data.CL_max_fullmodel.value;
    % -------------------------------------------------------------------------
    disp(" ")
    disp(" ++++ FIGURE 6 - LIFT CURVE INTERPOLATION ++++ ");
    Aircraft.Certification.Aerodynamic_data.CL_fullmodel_diagram.value = Lift_fullmodel_curve(Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.value, ...
        Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value, ...
        CL_wb, ...
        alpha_wb);    
    
    % SAVING FIGURES
    exportgraphics(Aircraft.Certification.Aerodynamic_data.CL_fullmodel_diagram.value, 'FullLiftModelInterpolation.pdf', 'ContentType', 'vector');
    exportgraphics(Aircraft.Certification.Aerodynamic_data.CL_fullmodel_diagram.value, 'FullLiftModelInterpolation.png', 'ContentType', 'vector');
    
    % Saving figures inside correct folder
    fprintf('Saving FullLiftModelInterpolation.pdf in: ');
    fprintf('\n'); 
    fprintf('%s\n', SaveFolder);
    % Moving file inside correct folder
    movefile FullLiftModelInterpolation.pdf Output
    movefile FullLiftModelInterpolation.png Output    
    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    a1 = p_CL_wb1(1);
    b1 = p_CL_wb1(2); 
    Aircraft.Certification.Aerodynamic_data.CLlin_PolCoeff_a.value = a1;
    Aircraft.Certification.Aerodynamic_data.CLlin_PolCoeff_b.value = b1;
   
    a2 = p_CL_wb2(1);
    b2 = p_CL_wb2(2); 
    c2 = p_CL_wb2(3);
    Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a.value = a2;
    Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b.value = b2;
    Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c.value = c2;
    % -------------------------------------------------------------------------
    CL0          = p_CL_wb1(2);
    CLalfa       = p_CL_wb1(1);
    Aircraft.Certification.Aerodynamic_data.CL0.value = CL0;
    Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value = CLalfa;
    alfa_zero_lift = alpha_fullmodel(0, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
    Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value           = alfa_zero_lift;
    Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.Attributes.unit = "deg";
    % -------------------------------------------------------------------------
    % CD INTERPOLATION
    % -------------------------------------------------------------------------
    p_cd_star   = polyfit(alpha_wb(index_aoa_star-1:index_aoa_star+1), CD_wb(index_aoa_star-1:index_aoa_star+1), 2);
    CD_star     = polyval(p_cd_star, alfa_star);
    Aircraft.Certification.Aerodynamic_data.CD_star.value           = CD_star;
    Aircraft.Certification.Aerodynamic_data.CD_star.Attributes.unit = "Non dimensional";
    % -------------------------------------------------------------------------
    n_pol1        = 2;
    p_CD_wb1      = polyfit([alpha_wb(1:index_aoa_star-1); alfa_star; alpha_wb(index_aoa_star:end)], ...
                            [CD_wb(1:index_aoa_star-1); CD_star; CD_wb(index_aoa_star:end)],   n_pol1);
    CD_wb_interp  = polyval(p_CD_wb1, alpha_interp);
    Aircraft.Certification.Aerodynamic_data.CD_PolCoeff_p.value           = p_CD_wb1;
    Aircraft.Certification.Aerodynamic_data.CD_PolCoeff_p.Attributes.unit = "Non dimensional";
    % -------------------------------------------------------------------------
    CDalfa = figure(150);
    hold on
    grid on 
    grid minor
    plot(alpha_interp, CD_wb_interp, '-r', 'LineWidth', 1.5)
    plot(alpha_wb, CD_wb, 'k.', 'MarkerSize', 10)
    xlim 'padded' ;
    ylim 'padded' ;
    xlabel("Angle of attack - $\alpha$ $(deg)$", "Interpreter", "latex")
    ylabel("Drag coefficient - $C_{D_{wb}}$", "Interpreter", "latex")
    title("Drag model", "Interpreter", "latex")
    legend({'Model','Data points'}, 'Interpreter', 'latex', 'Location', 'southeast')
    
    % SAVING FIGURES
    exportgraphics(CDalfa, 'FullDragModelInterpolation.pdf', 'ContentType', 'vector');
    exportgraphics(CDalfa, 'FullDragModelInterpolation.png', 'ContentType', 'vector');
    
    % Saving figures inside correct folder
    fprintf('Saving FullDragModelInterpolation.pdf in: ');
    fprintf('\n'); 
    fprintf('%s\n', SaveFolder);
    % Moving file inside correct folder
    movefile FullDragModelInterpolation.pdf Output
    movefile FullDragModelInterpolation.png Output        
    % -------------------------------------------------------------------------
    % CM INTERPOLATION
    % -------------------------------------------------------------------------
    n_pol     = 3;
    p_CM_wb   = polyfit(alpha_wb, CM_wb, n_pol);
    CM_interp = polyval(p_CM_wb, alpha_interp);
    Aircraft.Certification.Aerodynamic_data.CM_fullmodel.value = CM_interp;
    Aircraft.Certification.Aerodynamic_data.CM_fullmodel.Attributes.unit = "Non dimensional";
    % -------------------------------------------------------------------------   
    CMALFA = figure(151);
    hold on
    grid on 
    grid minor
    plot(alpha_interp, CM_interp, '-r', 'LineWidth', 1.5)
    plot(alpha_wb, CM_wb, 'k.', 'MarkerSize', 10)
    xlim 'padded' ;
    ylim 'padded' ;
    xlabel("Angle of attack - $\alpha$ $(deg)$", "Interpreter", "latex")
    ylabel("Pitch. mom. coefficient - $C_{M_{wb}}$", "Interpreter", "latex")
    title("Pitch. mom. model", "Interpreter", "latex")
    legend({'Model','Data points'}, 'Interpreter', 'latex', 'Location', 'southeast')
    
    % SAVING FIGURES
    exportgraphics(CDalfa, 'FullPitchMomModelInterpolation.pdf', 'ContentType', 'vector');
    exportgraphics(CDalfa, 'FullPitchMomModelInterpolation.png', 'ContentType', 'vector');
    
    % Saving figures inside correct folder
    fprintf('Saving FullPitchMomModelInterpolation.pdf in: ');
    fprintf('\n'); 
    fprintf('%s\n', SaveFolder);
    % Moving file inside correct folder
    movefile FullPitchMomModelInterpolation.pdf Output
    movefile FullPitchMomModelInterpolation.png Output    
    % -------------------------------------------------------------------------
    clear CL_wb CD_wb CM_wb alpha_wb alpha_interp CM_interp CD_wb_interp
    % -------------------------------------------------------------------------
    elseif  InputSource == "From User"
        disp("No data file. Insert data.")
end
