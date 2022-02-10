% ---------------------------------------------------------------------------
% INPUT SOURCE VARIABLE 
% ---------------------------------------------------------------------------
InputSource = "From File"; 
% % =========================================================================
% % Alpha_star and Alpha_max 
% % =========================================================================
% a       = -0.021837866;
% b       =  0.436064773;
% c       = -0.56312855;
% p_model = [a b c];
% Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a.value = a;
% Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b.value = b;
% Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c.value = c;
% % =========================================================================
% % ZERO LIFT COEFFICIENT
% CL0         = Aircraft.Certification.Aerodynamic_data.CL0.value;
% CL_alfa_deg = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;
% % ==== USEFUL FUNCTION DEFINED LOCALLY ====
% % -------------------------------------------------------------------------
% % CLMAX FUNCTION
% CLmax_func = @(rho, V, WS, n) (2 / rho) * (1 / V.^2) * (WS) * n;
% % -------------------------------------------------------------------------
% % CLMAX LINEAR
% CLmax_lin = @(CL_alfa, alfa) (CL0 + 0.01) + CL_alfa * alfa;
% % -------------------------------------------------------------------------
% % CLMAX NON LINEAR
% CLmax_non_lin = @(alfa) a * alfa^2 + b * alfa + c;
% % -------------------------------------------------------------------------
% % ALFA FUNCTION
% alfa_func = @(rho, S, V, WS, n, CLalfa, alfa_0lift) (2 / rho) * (1 / V.^2) * (1/CLalfa) * (WS) * n + alfa_0lift;
% % -------------------------------------------------------------------------
% % GUST LOAD FACTOR - POSITIVE FLIGHT
% nGust  = @(rho, V, a, kG, Ude, WS) 1 + (0.5 * rho * V * a * kG * Ude)/(WS); 
% % -------------------------------------------------------------------------
% % GUST LOAD FACTOR - INVERTED FLIGHT
% nGust_inverted  = @(rho, V, a, kG, Ude, WS) 1 - (0.5 * rho * V * a * kG * Ude)/(WS); 
% % -------------------------------------------------------------------------
% % STALL SPEED FUNCTION
% Vstall = @(WS, rho, CLmax, n) sqrt(WS * (2/rho) * (1/CLmax).*n); 
% % -------------------------------------------------------------------------
% 
% % =========================================================================
% % DRAG POLYNOMIAL COEFFICIENTS DEFINITION 
% % =========================================================================
% Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value = 0.079;                       % Coefficient inside an expression for the CD in polynomial form
% Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.Attributes.unit = "Non dimensional";
% Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value = 0.365;                       % Coefficient inside an expressione for the CD in polynomial form
% Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.Attributes.unit = "Non dimensional"; 
% k1 = Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value;
% k2 = Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value;

% % =========================================================================
% % OSWALDT EFFICIENCY FACTOR
% e   = Aircraft.Certification.Aerodynamic_data.e.value;
% 
% % ASPECT RATIO 
% AR  = Aircraft.Geometry.Wing.AR.value;
% 
% % ZERO-LIFT DRAG COEFFICIENT
% CD0 = Aircraft.Certification.Aerodynamic_data.CD0.value;
% 
% % ENDING OF LINEAR PART OF LIFT CURVE
% CL_star = Aircraft.Certification.Aerodynamic_data.CL_star.value;
% 
% % CLALFA IN DEG^-1
% CLalfa = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;
% % =========================================================================

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
% 
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% CL_WB_model = @(alpha) a*alpha.^2 + b*alpha + c; 
% alpha_plus  = @(CL) (-b + sqrt(b^2 - 4*a*(c - CL)))/(2*a);
% alpha_meno  = @(CL) (-b - sqrt(b^2 - 4*a*(c - CL)))/(2*a);
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% % STORE INSIDE STRUCT VARIABLE
% Aircraft.Certification.Aerodynamic_data.alpha_star.value            = alpha_star;
% Aircraft.Certification.Aerodynamic_data.alpha_star.Attributes.unit  = "deg";
% Aircraft.Certification.Aerodynamic_data.alpha_max.value             = alpha_max;
% Aircraft.Certification.Aerodynamic_data.alpha_max.Attributes.unit   = "deg";
% 
% alpha_i = - 4.0;
% alpha_f = 13.0;
% alpha_interpolation_interval = linspace(alpha_i, alpha_f, numb)';
% 
% % Full lift model 
% CL_interpolation = zeros(length(alpha_interpolation_interval), 1);
% for i = 1:length(alpha_interpolation_interval)
%     CL_interpolation(i) = CL_fullmodel(alpha_interpolation_interval(i));
% end
% Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value = CL_interpolation;
% Aircraft.Certification.Aerodynamic_data.CL_fullmodel.Attributes.unit = "Non dimensional";
% Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.value = alpha_interpolation_interval;
% Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.Attributes.unit = "deg";
if InputSource == "From File"
    % -------------------------------------------------------------------------
    % ALFA FUNCTION
    alfa_func = @(rho, S, V, WS, n, CLalfa, alfa_0lift) (2 / rho) * (1 / V.^2) * (1/CLalfa) * (WS) * n + alfa_0lift;
    % =========================================================================
    CL_wb    = str2num(Aircraft.Certification.Aerodynamic_data.CL.value);
    CD_wb    = str2num(Aircraft.Certification.Aerodynamic_data.CD.value); 
    CM_wb    = str2num(Aircraft.Certification.Aerodynamic_data.CM.value);
    alpha_wb = str2num(Aircraft.Certification.Aerodynamic_data.alpha.value);
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
    Aircraft.Certification.Aerodynamic_data.alfa_zero_lift.value           = alfa_0l;
    Aircraft.Certification.Aerodynamic_data.alfa_zero_lift.Attributes.unit = "deg";
    % =========================================================================
    % -------------------------------------------------------------------------
    % CL INTERPOLATION
    % -------------------------------------------------------------------------
    alpha_interp  = linspace(alpha_wb(1), alpha_wb(end), numb)';
    % -------------------------------------------------------------------------
    p_alfa_star   = polyfit(CL_wb(4:5), alpha_wb(4:5), 2);
    alfa_star     = polyval(p_alfa_star, CL_star);
    Aircraft.Certification.Aerodynamic_data.alpha_star.value = alfa_star;
    Aircraft.Certification.Aerodynamic_data.Attributes.unit  = "deg";
    % -------------------------------------------------------------------------
    CL_wb_interp  = zeros(length(alpha_interp), 1);
    n_pol1        = 1;
    n_pol2        = 2;
    p_CL_wb1      = polyfit([alpha_wb(1:3); alfa_star],   [CL_wb(1:3); CL_star],   n_pol1);
    p_CL_wb2      = polyfit([alfa_star; alpha_wb(4:end)], [CL_star; CL_wb(4:end)], n_pol2);
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
    Aircraft.Certification.Aerodynamic_data.CL_max.value = Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value;
    Aircraft.Certification.Aerodynamic_data.CL_max.Attributes.unit = "Non dimensional";
    CL_max = Aircraft.Certification.Aerodynamic_data.CL_max.value;
    % -------------------------------------------------------------------------
    disp(" ")
    disp(" ++++ FIGURE 6 - LIFT CURVE INTERPOLATION ++++ ");
    Aircraft.Certification.Aerodynamic_data.CL_fullmodel_diagram.value = Lift_fullmodel_curve(Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.value, ...
        Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value, ...
        CL_wb, ...
        alpha_wb);
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
    % CD INTERPOLATION
    % -------------------------------------------------------------------------
    p_cd_star   = polyfit(alpha_wb(4:5), CD_wb(4:5), 2);
    CD_star     = polyval(p_cd_star, alfa_star);
    % -------------------------------------------------------------------------
    n_pol1        = 2;
    p_CD_wb1      = polyfit([alpha_wb(1:3); alfa_star; alpha_wb(4:end)], ...
                            [CD_wb(1:3); CD_star; CD_wb(4:end)],   n_pol1);
    CD_wb_interp  = polyval(p_CD_wb1, alpha_interp);
    % -------------------------------------------------------------------------
    CDalfa = figure(153);
    hold on
    grid on 
    grid minor
    plot(alpha_interp, CD_wb_interp, '-r', 'LineWidth', 1.5)
    plot(alpha_wb, CD_wb, 'k.', 'MarkerSize', 10)
    xlabel("Angle of attack - $\alpha$ $(deg)$", "Interpreter", "latex")
    ylabel("Drag coefficient - $C_{D_{wb}}$", "Interpreter", "latex")
    title("Drag full model", "Interpreter", "latex")
    legend({'Full model','Data points'}, 'Interpreter', 'latex', 'Location', 'southeast')
    % -------------------------------------------------------------------------    
    elseif  InputSource == "From User"
        disp("No data file. Insert data.")
end

% % -------------------------------------------------------------------------
% % CL INTERPOLATION
% % -------------------------------------------------------------------------
% alpha_interp = linspace(alpha_wb(1), alpha_wb(end), numb)';
% n_pol3       = 3;
% p_CL_wb3     = polyfit(alpha_wb, CL_wb, n_pol3);
% CL_wb_interp = polyval(p_CL_wb3, alpha_interp);
% Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.value = alpha_interp;
% Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.Attributes.unit = "deg";
% Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value = CL_wb_interp;
% Aircraft.Certification.Aerodynamic_data.CL_fullmodel.Attributes.unit = "Non dimensional";
% % -------------------------------------------------------------------------
% disp(" ")
% disp(" ++++ FIGURE 6 - LIFT CURVE INTERPOLATION ++++ ");
% Aircraft.Certification.Aerodynamic_data.CL_fullmodel_diagram.value = Lift_fullmodel_curve(Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.value, ...
%     Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value, ...
%     CL_wb, ...
%     alpha_wb);
% 
% % -------------------------------------------------------------------------
% % CD INTERPOLATION
% % -------------------------------------------------------------------------
% n_pol         = 3;
% p_CD_wb1      = polyfit(CL_wb(1:5),   CD_wb(1:5),   n_pol);
% p_CD_wb2      = polyfit(CL_wb(5:end), CD_wb(5:end), n_pol);
% CL_aux1       = linspace(CL_wb(1), CL_wb(5),   numb*0.5);
% CL_aux2       = linspace(CL_wb(5), CL_wb(end), numb*0.5);
% CL_aux        = [CL_aux1 CL_aux2];
% CD_CL_interp1 = polyval(p_CD_wb1, CL_aux1);
% CD_CL_interp2 = polyval(p_CD_wb2, CL_aux2);
% CD_CL_interp  = [CD_CL_interp1 CD_CL_interp2];
% % CD_CL_interp = polyval(p_CD_wb, CL_aux);
% Aircraft.Certification.Aerodynamic_data.CL_aux_fullmodel.value = CL_aux;
% Aircraft.Certification.Aerodynamic_data.CL_aux_fullmodel.Attributes.unit = "Non dimensional";
% Aircraft.Certification.Aerodynamic_data.CD_CLfullmodel.value = CD_CL_interp;
% Aircraft.Certification.Aerodynamic_data.CD_CLfullmodel.Attributes.unit = "Non dimensional";
% % -------------------------------------------------------------------------
% % -------------------------------------------------------------------------
% % CD INTERPOLATION
% % -------------------------------------------------------------------------
% n_pol        = 3;
% p_alpha_cd   = polyfit(alpha_wb, CD_wb, n_pol);
% alpha_interp = linspace(alpha_wb(1), alpha_wb(end), numb);
% CD_wb_interp = polyval(p_alpha_cd, alpha_interp);
% Aircraft.Certification.Aerodynamic_data.CD_wb_interp.value = CD_wb_interp;
% Aircraft.Certification.Aerodynamic_data.CD_wb_interp.Attributes.unit = "Non dimensional";
% % -------------------------------------------------------------------------
% % CM INTERPOLATION
% % -------------------------------------------------------------------------
% n_pol     = 3;
% p_CM_wb   = polyfit(alpha_wb, CM_wb, n_pol);
% CM_interp = polyval(p_CM_wb, alpha_interp);
% Aircraft.Certification.Aerodynamic_data.CM_fullmodel.value = CM_interp;
% Aircraft.Certification.Aerodynamic_data.CM_fullmodel.Attributes.unit = "Non dimensional";
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
% % Alpha zero lift
% Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value = alpha_calc_lin(obj1, ...
%                                                                                0.0, ...
%                                                                                Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                                                                Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
% Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.Attributes.unit = "deg";
% 
% alpha_zerol = Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
% 
% % alpha_star  = alpha_plus(Aircraft.Certification.Aerodynamic_data.CL_star.value);
% % alpha_max   = alpha_meno(Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value);
% 
% % -------------------------------------------------------------------------
% % ALPHA INTERPOLATION
% % -------------------------------------------------------------------------
% n_pol1           = 2;
% n_pol2           = 2;
% p_alpha_wb1      = polyfit(CL_wb(1:5),   alpha_wb(1:5), n_pol1);
% p_alpha_wb2      = polyfit(CL_wb(5:end), alpha_wb(5:end), n_pol);
% CL_aux1          = linspace(CL_wb(1), CL_wb(5),   numb*0.5);
% CL_aux2          = linspace(CL_wb(5), CL_wb(end), numb*0.5);
% CL_aux           = [CL_aux1 CL_aux2];
% alpha_interp1    = polyval(p_alpha_wb1, CL_aux1);
% alpha_interp2    = polyval(p_alpha_wb2, CL_aux2);
% alpha_interp_AoA = [alpha_interp1 alpha_interp2];
% % CD_CL_interp = polyval(p_CD_wb, CL_aux);
% Aircraft.Certification.Aerodynamic_data.CL_alpha_aux_fullmodel.value = CL_aux;
% Aircraft.Certification.Aerodynamic_data.CL_alpha_aux_fullmodel.Attributes.unit = "Non dimensional";
% Aircraft.Certification.Aerodynamic_data.alpha_interp_AoA.value = alpha_interp_AoA;
% Aircraft.Certification.Aerodynamic_data.alpha_interp_AoA.Attributes.unit = "Non dimensional";
% % -------------------------------------------------------------------------
% 
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % disp(" ")
% % disp(" ++++ FIGURE 6 - LIFT CURVE INTERPOLATION ++++ ");
% % Aircraft.Certification.Aerodynamic_data.CL_fullmodel_diagram.value = Lift_fullmodel_curve(Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.value, ...
% %     Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value, ...
% %     CL_wb, ...
% %     alpha_wb);
% 
% % % SAVING FIGURES
% % exportgraphics(Aircraft.Certification.Aerodynamic_data.CL_fullmodel_diagram.value, 'FullLiftModelInterpolation.pdf', 'ContentType', 'vector');
% % exportgraphics(Aircraft.Certification.Aerodynamic_data.CL_fullmodel_diagram.value, 'FullLiftModelInterpolation.png', 'ContentType', 'vector');
% % 
% % % Saving figures inside correct folder
% % fprintf('Saving FullLiftModelInterpolation.pdf in: ');
% % fprintf('\n'); 
% % fprintf('%s\n', SaveFolder);
% % % Moving file inside correct folder
% % movefile FullLiftModelInterpolation.pdf Output
% % movefile FullLiftModelInterpolation.png Output
% 
% %% MAX LIFT COEFFICIENT - WING BODY CONFIGURATION
% % -------------------------------------------------------------------------
% % Maximum lift coefficient for the wing - body configuration 
% % A complete documentation of the function CLMax(...) used here is inside
% % the class file aero_model.m, which can be found inside the 'utilities'
% % folder of this library.
% Aircraft.Certification.Aerodynamic_data.CLMAX_wingbody.value = max(Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value); 
% Aircraft.Certification.Aerodynamic_data.CLMAX_wingbody.Attributes.unit = "Non dimensional";
% % -------------------------------------------------------------------------
% CDCL = figure(150);
% hold on
% grid on 
% grid minor
% plot(CD_CL_interp, CL_aux, '-r', 'LineWidth', 1.5)
% % plot(AoA, CL_inverted, '-b', 'LineWidth', 1.5)
% plot(CD_wb, CL_wb, 'k.', 'MarkerSize', 10)
% % xlim([]);
% % ylim([0.0 1.80]);
% xlabel("Drag coefficient - $C_{D_{wb}}$", "Interpreter", "latex")
% ylabel("Lift coefficient - $C_{L_{wb}}$", "Interpreter", "latex")
% title("Drag full model", "Interpreter", "latex")
% legend({'Full model','Data points'}, 'Interpreter', 'latex', 'Location', 'southeast')
% % -------------------------------------------------------------------------
% % -------------------------------------------------------------------------
% CMalpha = figure(151);
% hold on
% grid on 
% grid minor
% plot(alpha_interp, CM_interp, '-r', 'LineWidth', 1.5)
% % plot(AoA, CL_inverted, '-b', 'LineWidth', 1.5)
% plot(alpha_wb, CM_wb, 'k.', 'MarkerSize', 10)
% % xlim([]);
% % ylim([0.0 1.80]);
% xlabel("Angle of attack - $\alpha$ $(deg)$", "Interpreter", "latex")
% ylabel("Pitching mom. coefficient - $C_{M_{wb}}$", "Interpreter", "latex")
% title("Pitching mom. full model", "Interpreter", "latex")
% legend({'Full model','Data points'}, 'Interpreter', 'latex', 'Location', 'southeast')
% % -------------------------------------------------------------------------
% % -------------------------------------------------------------------------
% alphamodel = figure(152);
% hold on
% grid on 
% grid minor
% plot(Aircraft.Certification.Aerodynamic_data.CL_alpha_aux_fullmodel.value, ...
%      Aircraft.Certification.Aerodynamic_data.alpha_interp_AoA.value, '-r', 'LineWidth', 1.5)
% % plot(AoA, CL_inverted, '-b', 'LineWidth', 1.5)
% plot(CL_wb, alpha_wb, 'k.', 'MarkerSize', 10)
% % xlim([]);
% % ylim([0.0 1.80]);
% xlabel("Lift coefficient - $C_{L_{wb}}$", "Interpreter", "latex")
% ylabel("Angle of attack - $\alpha$ $(deg)$", "Interpreter", "latex")
% title("Test", "Interpreter", "latex")
% legend({'Full model','Data points'}, 'Interpreter', 'latex', 'Location', 'southeast')
% % -------------------------------------------------------------------------
% % -------------------------------------------------------------------------
% CDalfa = figure(153);
% hold on
% grid on 
% grid minor
% plot(alpha_interp, CD_wb_interp, '-r', 'LineWidth', 1.5)
% % plot(AoA, CL_inverted, '-b', 'LineWidth', 1.5)
% plot(alpha_wb, CD_wb, 'k.', 'MarkerSize', 10)
% % xlim([]);
% % ylim([0.0 1.80]);
% xlabel("Angle of attack - $\alpha$ $(deg)$", "Interpreter", "latex")
% ylabel("Drag coefficient - $C_{D_{wb}}$", "Interpreter", "latex")
% title("Drag full model", "Interpreter", "latex")
% legend({'Full model','Data points'}, 'Interpreter', 'latex', 'Location', 'southeast')
% % -------------------------------------------------------------------------
% 
% % INSTANTIATION OF SOME CONSTANTS IMPORTANT TO MOMENT CONTRIB. CALCULATIONS
% XAC         = Aircraft.Geometry.General.XAC_nondim.value;
% XCG         = Aircraft.Geometry.General.XAC_nondim.value;
% bCG         = Aircraft.Geometry.General.bcg.value;
% MAC         = Aircraft.Geometry.Wing.mac.value;
% Thrust_axes = Aircraft.Geometry.Engine.Primary.Thrust_axes.value;
% l_ht        = Aircraft.Geometry.Horizontal.l.value;
% 
% % OTHER AERODYNAMIC COEFFICIENTS
% CM_landing_gear = Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value;
% CM0             = Aircraft.Certification.Aerodynamic_data.CM0.value;
% CM_slope        = Aircraft.Certification.Aerodynamic_data.CMCL.value;
% 
% % %% TAIL BALANCING LOADS - N = 1.0 
% % 
% % V_unit_load_factor = zeros(length(Positive_stall_speed), 1);
% % VS = Positive_stall_speed(1);
% % V_unit_load_factor(1) = VS;
% % VD = V_fromfgtoD(end);
% % for i = 2:length(V_unit_load_factor)
% %         V_unit_load_factor(i) = V_unit_load_factor(i-1) + (VD - VS)*(1/length(V_unit_load_factor));
% % end
% % 
% % % DYNAMIC PRESSURE
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.value = 0.5*Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value.^2);
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.Attributes.unit = "Pa";
% % 
