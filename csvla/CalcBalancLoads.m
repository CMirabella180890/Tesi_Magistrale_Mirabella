% =========================================================================
% COEFFICIENT OF THE AERODYNAMIC MODEL
% =========================================================================
alpha_wb                                                       = str2num(Aircraft.Certification.Aerodynamic_data.alpha.value);
p_CD_wb1                                                       = Aircraft.Certification.Aerodynamic_data.CD_PolCoeff_p.value;
Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a.value = a2; 
Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b.value = b2; 
Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c.value = c2; 
Aircraft.Certification.Aerodynamic_data.CL_max.value           = CL_max;
CL_star      = Aircraft.Certification.Aerodynamic_data.CL_star.value;
alpha_interp = linspace(alpha_wb(1), alpha_wb(end), numb)'; 
p_CL_wb1     = [Aircraft.Certification.Aerodynamic_data.CLlin_PolCoeff_a.value; Aircraft.Certification.Aerodynamic_data.CLlin_PolCoeff_b.value];
CL0          = p_CL_wb1(2);
CLalfa       = p_CL_wb1(1);
clear alpha_wb
CL_max_fullmodel = Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value;
% CL_max_wingbody  = Aircraft.Certification.Aerodynamic_data.CL_max_fullmodel.value;
% CL_max_fullmodel = Aircraft.Certification.Aerodynamic_data.CL_max_fullmodel.value;
CL_max_inv      = Aircraft.Certification.Aerodynamic_data.Max_Inverted_Lift_Coefficient.value;
% =========================================================================
alfa_star = Aircraft.Certification.Aerodynamic_data.alpha_star.value;
p_CL1     = Aircraft.Certification.Aerodynamic_data.p_CL_wb1.value;
p_CL2     = Aircraft.Certification.Aerodynamic_data.p_CL_wb2.value;
% =========================================================================
% % a       = -0.021837866;
% % b       =  0.436064773;
% % c       = -0.56312855;
% % p_model = [a b c];
% % Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a.value = a;
% % Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b.value = b;
% % Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c.value = c;
% % =========================================================================
% % ZERO LIFT COEFFICIENT
% CL0         = Aircraft.Certification.Aerodynamic_data.CL0.value;
% CL_alfa_deg = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;
% % ==== USEFUL FUNCTION DEFINED LOCALLY ====
% % -------------------------------------------------------------------------
% CLMAX FUNCTION
CLmax_func = @(rho, V, WS, n) (2 / rho) * (1 / V.^2) * (WS) * n;
% % -------------------------------------------------------------------------
% % CLMAX LINEAR
% CLmax_lin = @(CL_alfa, alfa) (CL0 + 0.01) + CL_alfa * alfa;
% % -------------------------------------------------------------------------
% % CLMAX NON LINEAR
% CLmax_non_lin = @(alfa) a * alfa^2 + b * alfa + c;
% -------------------------------------------------------------------------
% ALFA FUNCTION
alfa_func = @(rho, S, V, WS, n, CLalfa, alfa_0lift) (2 / rho) * (1 / V.^2) * (1/CLalfa) * (WS) * n + alfa_0lift;
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
% % 
% % % =========================================================================
% % % DRAG POLYNOMIAL COEFFICIENTS DEFINITION 
% % % =========================================================================
% % Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value = 0.079;                       % Coefficient inside an expression for the CD in polynomial form
% % Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.Attributes.unit = "Non dimensional";
% % Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value = 0.365;                       % Coefficient inside an expressione for the CD in polynomial form
% % Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.Attributes.unit = "Non dimensional"; 
% % k1 = Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value;
% % k2 = Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value;
% 
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
% 
% %% Script to evaluate balancing horizontal tail loads
% %   DESCRIPTION
% %    In this script, all the methods relative to the calculation of
% %    aerodynamic coefficients, necessary to the correct estimation of the
% %    airloads acting on the horizontal empennage. 
% 
% %% RETURN INSIDE UTILITIES
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
% obj1 = aero_model; 
% % 
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
% % 
% % % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % CL_WB_model = @(alpha) a*alpha.^2 + b*alpha + c; 
% % alpha_plus  = @(CL) (-b + sqrt(b^2 - 4*a*(c - CL)))/(2*a);
% % alpha_meno  = @(CL) (-b - sqrt(b^2 - 4*a*(c - CL)))/(2*a);
% % % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
% % % STORE INSIDE STRUCT VARIABLE
% % Aircraft.Certification.Aerodynamic_data.alpha_star.value            = alpha_star;
% % Aircraft.Certification.Aerodynamic_data.alpha_star.Attributes.unit  = "deg";
% % Aircraft.Certification.Aerodynamic_data.alpha_max.value             = alpha_max;
% % Aircraft.Certification.Aerodynamic_data.alpha_max.Attributes.unit   = "deg";
% % 
% % alpha_i = - 4.0;
% % alpha_f = 13.0;
% % alpha_interpolation_interval = linspace(alpha_i, alpha_f, numb)';
% % 
% % % Full lift model 
% % CL_interpolation = zeros(length(alpha_interpolation_interval), 1);
% % for i = 1:length(alpha_interpolation_interval)
% %     CL_interpolation(i) = CL_fullmodel(alpha_interpolation_interval(i));
% % end
% % Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value = CL_interpolation;
% % Aircraft.Certification.Aerodynamic_data.CL_fullmodel.Attributes.unit = "Non dimensional";
% % Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.value = alpha_interpolation_interval;
% % Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.Attributes.unit = "deg";
% CL_wb    = str2num(Aircraft.Certification.Aerodynamic_data.CL.value);
% CD_wb    = str2num(Aircraft.Certification.Aerodynamic_data.CD.value); 
% CM_wb    = str2num(Aircraft.Certification.Aerodynamic_data.CM.value);
% alpha_wb = str2num(Aircraft.Certification.Aerodynamic_data.alpha.value);
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
% INSTANTIATION OF SOME CONSTANTS IMPORTANT TO MOMENT CONTRIB. CALCULATIONS
XAC         = Aircraft.Geometry.General.XAC_nondim.value;
XCG         = Aircraft.Geometry.General.XAC_nondim.value;
bCG         = Aircraft.Geometry.General.bcg.value;
MAC         = Aircraft.Geometry.Wing.mac.value;
Thrust_axes = Aircraft.Geometry.Engine.Primary.Thrust_axes.value;
l_ht        = Aircraft.Geometry.Horizontal.l.value;

% OTHER AERODYNAMIC COEFFICIENTS
CM_landing_gear = Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value;
CM0             = Aircraft.Certification.Aerodynamic_data.CM0.value;
CM_slope        = Aircraft.Certification.Aerodynamic_data.CMCL.value;

% %% TAIL BALANCING LOADS - N = 1.0 
% 
% V_unit_load_factor = zeros(length(Positive_stall_speed), 1);
% VS = Positive_stall_speed(1);
% V_unit_load_factor(1) = VS;
% VD = V_fromfgtoD(end);
% for i = 2:length(V_unit_load_factor)
%         V_unit_load_factor(i) = V_unit_load_factor(i-1) + (VD - VS)*(1/length(V_unit_load_factor));
% end
% 
% % DYNAMIC PRESSURE
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.value = 0.5*Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value.^2);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.Attributes.unit = "Pa";
% % 

%% BALANCING LOADS CALCULATIONS - POSITIVE LOAD FACTORS
% ------------------------------------------------------------------------- 
% Calculation of the CL for the wing - body configuration. It must be
% noticed that the calculations relative to lift coefficient in this
% particular range of airspeed V and load factor n are going to give a
% constant as a results, which is the maximum lift coefficient following
% the lift curve of the aircrat. The previously defined CLMAX is a little
% bit greater than the one obtained from V - n diagram flight condition to
% take into account the extra lift produced by the horizontal empennage.
% 
% ++ AIRWORTHINESS RULES ++
%% CS - VLA 421 Balancing loads 
%  (a) A horizontal tail balancing load is a load necessary to maintain
%      equilibrium in any specified flight condition with no pitching
%      acceleration.
%  (b) Horizontal tail surfaces must be designed for the balancing loads
%      occurring at any point on the limit manoeuvring envelope and in the
%      flap conditions specified in CS - VLA 345. The distribution in
%      figure B6 of Appendix B may be used.

Straight_flight_Case = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Straight_flight.value;

switch (Straight_flight_Case)
    % CASE 1: VA greater than the intercept
    case 'Case 1'
        % CL, ALFA, CD CALCULATIONS - POSITIVE LOAD FACTORS
        % ----------------------------------------------------------------- 
        % Calculation of the CL for the wing - body configuration. It must be
        % noticed that the calculations relative to lift coefficient in this
        % particular range of airspeed V and load factor n are going to give a
        % constant as a results, which is the maximum lift coefficient following
        % the lift curve of the aircrat. The previously defined CLMAX is a little
        % bit greater than the one obtained from V - n diagram flight condition to
        % take into account the extra lift produced by the horizontal empennage.
        % WING BODY LIFT CALCULATION
        % ----------------------------------------------------------------- 
        % In this section of the code the following formula is applied to evaluate
        % the wing-body lift with the formula 
        %        
        %   L = 0.5*rho*(V^2)*Sw*CL 
        %
        % using the previously calculated lift coefficients, coming from the 
        % various curves of the final envelope diagram.
        % =================================================================
        % FROM 0 TO S
        CL_from0toS   = zeros(length(V_from0toS), 1);
        alfa_from0toS = zeros(length(V_from0toS), 1);
        CD_from0toS   = zeros(length(V_from0toS), 1);
        q_from0toS    = zeros(length(V_from0toS), 1);
        WBL_from0toS  = zeros(length(V_from0toS), 1);
        WBL_from0toS_new  = zeros(length(V_from0toS), 1);
        CMCL_from0toS = zeros(length(V_from0toS), 1);
        CMCD_from0toS = zeros(length(V_from0toS), 1);
        CMCT_from0toS = zeros(length(V_from0toS), 1);
        CMCG_from0toS = zeros(length(V_from0toS), 1); 
        CLHT_from0toS = zeros(length(V_from0toS), 1);
        LHT_from0toS  = zeros(length(V_from0toS), 1);
        % =================================================================
        CL_from0toS_fullvehicle   = zeros(length(V_from0toS), 1);
        alfa_from0toS_fullvehicle = zeros(length(V_from0toS), 1);
        % =================================================================
        for i = 1:length(V_from0toS)
            CL_from0toS(i)   = CLmax_func(rho0, V_from0toS(i), WS, n_from0toS(i));
            if CL_from0toS(i) < CL_max_fullmodel
                CL_from0toS(i) = CL_from0toS(i);
            elseif CL_from0toS(i) > CL_max_fullmodel
                CL_from0toS(i) = CL_max_fullmodel;
            end
            alfa_from0toS(i) = alpha_fullmodel(CL_from0toS(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
            % alfa_from0toS(i) = alpha_fullmodel(CL_from0toS(i), a2, b2, c2, CL_max_wingbody, CL_star, CL0, CLalfa);
            CD_from0toS(i)   = polyval(p_CD_wb1, alfa_from0toS(i));
            q_from0toS(i)    = 0.5*rho0*(V_from0toS(i))^2;
            WBL_from0toS(i)  = q_from0toS(i)*S*CL_from0toS(i)*1e-1;
            CMCL_from0toS(i) = CLWB_contrib(obj1, CL_from0toS(i), deg2rad(alfa_from0toS(i)), XAC, XCG, bCG, MAC);
            CMCD_from0toS(i) = CDWB_contrib(obj1, CD_from0toS(i), deg2rad(alfa_from0toS(i)), XAC, XCG, bCG, MAC);
            CMCT_from0toS(i) = CT_contr(obj1, CD_from0toS(i), Thrust_axes, MAC);
            CMCG_from0toS(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_from0toS(i));
            % HORIZONTAL TAIL LIFT COEFFICIENT
            CLHT_from0toS(i) = CL_Tail(obj1, CMCL_from0toS(i), CMCD_from0toS(i), ...
                                             CMCT_from0toS(i), CMCG_from0toS(i), ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_from0toS(i)));
            % HORIZONTAL TAIL LIFT
            LHT_from0toS(i) = (0.5)*(V_from0toS(i)^2)*(S)*(rho0)*(CLHT_from0toS(i))*(1e-1);
            
            % DRAFT VERSION OF ITERATION 
            CL_tail         = CLHT_from0toS(i);
            CL_wb           = CL_from0toS(i);
            CL_new_from0toS = CL_wb - CL_tail;
            tol             = 1e-3; 
            n               = 1;
            while abs(CL_wb - CL_tail) > tol
               alfa_new_from0toS = alpha_fullmodel(CL_new_from0toS, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
               CD_wb             = polyval(p_CD_wb1, alfa_new_from0toS);
               CMCL_new          = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_from0toS), XAC, XCG, bCG, MAC);
               CMCD_new          = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_from0toS), XAC, XCG, bCG, MAC);
               CMCT_new          = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
               CMCG_new          = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
               CLHT_new          = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_new_from0toS));
               CL_tail           = CLHT_new;
               CL_new_from0toS   = CL_from0toS(i) - CL_tail;
               CL_wb             = CL_wb + (CL_new_from0toS - CL_wb) * 1e-1;
               n                 = n + 1;
               if n == 15
                   break
               end
            end
            CL_from0toS_fullvehicle(i)   = CL_new_from0toS; 
%             CL_from0toS(i)   = CL_new_from0toS; 
            CLHT_from0toS(i) = CL_tail;
            alfa_from0toS_fullvehicle(i) = alfa_new_from0toS;
%             alfa_from0toS(i) = alfa_new_from0toS;
            WBL_from0toS(i)  = q_from0toS(i) * S * CL_from0toS(i) * 1e-1;
            WBL_from0toS_new(i)  = q_from0toS(i) * S * CL_from0toS_fullvehicle(i) * 1e-1;
            LHT_from0toS(i)  = (0.5)*(V_from0toS(i)^2)*(S)*(rho0)*(CLHT_from0toS(i))*(1e-1);
            CMCL_from0toS(i) = CMCL_new;
            CMCD_from0toS(i) = CMCD_new;
            CMCG_from0toS(i) = CMCG_new;
        end 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS.value = CL_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS_new.value = CL_from0toS_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS.value = alfa_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS_rad.value = deg2rad(alfa_from0toS);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS_rad.Attributes.unit = "rad";    
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS_new.value = alfa_from0toS_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS_new_rad.value = deg2rad(alfa_from0toS_fullvehicle);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS_new_rad.Attributes.unit = "rad";     
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toS.value = CD_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toS.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toS.value = q_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toS.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toS.value = WBL_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toS.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.value = WBL_from0toS_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_from0toS.value = CMCL_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_from0toS.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_from0toS.value = CMCD_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_from0toS.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_from0toS.value = CMCT_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_from0toS.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_from0toS.value = CMCG_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_from0toS.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_from0toS.value = CLHT_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_from0toS.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_from0toS.value = LHT_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_from0toS.Attributes.unit = "daN";        
        % ================================================================= 
        % FROM S TO A1
        CL_fromStoA1   = zeros(length(V_fromStoA1), 1);
        alfa_fromStoA1 = zeros(length(V_fromStoA1), 1);
        CD_fromStoA1   = zeros(length(V_fromStoA1), 1);
        q_fromStoA1    = zeros(length(V_fromStoA1), 1);
        WBL_fromStoA1  = zeros(length(V_fromStoA1), 1);
        WBL_fromStoA1_new  = zeros(length(V_fromStoA1), 1);
        CMCL_fromStoA1 = zeros(length(V_fromStoA1), 1);
        CMCD_fromStoA1 = zeros(length(V_fromStoA1), 1);
        CMCT_fromStoA1 = zeros(length(V_fromStoA1), 1);
        CMCG_fromStoA1 = zeros(length(V_fromStoA1), 1); 
        CLHT_fromStoA1 = zeros(length(V_fromStoA1), 1);
        LHT_fromStoA1  = zeros(length(V_fromStoA1), 1);
        % =================================================================
        CL_fromStoA1_fullvehicle   = zeros(length(V_fromStoA1), 1);
        alfa_fromStoA1_fullvehicle = zeros(length(V_fromStoA1), 1);
        % =================================================================
        for i = 1:length(V_fromStoA1)
            CL_fromStoA1(i)   = CLmax_func(rho0, V_fromStoA1(i), WS, n_fromStoA1(i));
            if CL_fromStoA1(i) < CL_max_fullmodel
                CL_fromStoA1(i) = CL_fromStoA1(i);
            elseif CL_fromStoA1(i) > CL_max_fullmodel
                CL_fromStoA1(i) = CL_max_fullmodel;
            end
            alfa_fromStoA1(i) = alpha_fullmodel(CL_fromStoA1(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
            CD_fromStoA1(i)   = polyval(p_CD_wb1, alfa_fromStoA1(i));
            q_fromStoA1(i)    = 0.5*rho0*(V_fromStoA1(i))^2;
            WBL_fromStoA1(i)  = q_fromStoA1(i)*S*CL_fromStoA1(i)*1e-1;
            CMCL_fromStoA1(i) = CLWB_contrib(obj1, CL_fromStoA1(i), deg2rad(alfa_fromStoA1(i)), XAC, XCG, bCG, MAC);
            CMCD_fromStoA1(i) = CDWB_contrib(obj1, CL_fromStoA1(i), deg2rad(alfa_fromStoA1(i)), XAC, XCG, bCG, MAC);
            CMCT_fromStoA1(i) = CT_contr(obj1, CD_fromStoA1(i), Thrust_axes, MAC);
            CMCG_fromStoA1(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromStoA1(i));
            % HORIZONTAL TAIL LIFT COEFFICIENT
            CLHT_fromStoA1(i) = CL_Tail(obj1, CMCL_fromStoA1(i), CMCD_fromStoA1(i), ...
                                             CMCT_fromStoA1(i), CMCG_fromStoA1(i), ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_fromStoA1(i)));
            % HORIZONTAL TAIL LIFT
            LHT_fromStoA1(i) = (0.5)*(V_fromStoA1(i)^2)*(S)*(rho0)*(CLHT_fromStoA1(i))*(1e-1);
            
            % DRAFT VERSION OF ITERATION 
            CL_tail          = CLHT_fromStoA1(i);
            CL_wb            = CL_fromStoA1(i);
            CL_new_fromStoA1 = CL_wb - CL_tail;
            tol              = 1e-3; 
            n                = 1;
            while abs(CL_wb - CL_tail) > tol
               alfa_new_fromStoA1 = alpha_fullmodel(CL_new_fromStoA1, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
               CD_wb             = polyval(p_CD_wb1, alfa_new_fromStoA1);
               CMCL_new          = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromStoA1), XAC, XCG, bCG, MAC);
               CMCD_new          = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromStoA1), XAC, XCG, bCG, MAC);
               CMCT_new          = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
               CMCG_new          = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
               CLHT_new          = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromStoA1));
               CL_tail           = CLHT_new;
               CL_new_fromStoA1  = CL_fromStoA1(i) - CL_tail;
               CL_wb             = CL_wb + (CL_new_fromStoA1 - CL_wb) * 1e-1;
               n                 = n + 1;
               if n == 15
                   break
               end
            end
            CL_fromStoA1_fullvehicle(i)   = CL_new_fromStoA1; 
%             CL_fromStoA1(i)   = CL_new_fromStoA1; 
            CLHT_fromStoA1(i) = CL_tail;
            alfa_fromStoA1_fullvehicle(i) = alfa_new_fromStoA1;
%             alfa_fromStoA1(i) = alfa_new_fromStoA1;
            WBL_fromStoA1(i)  = q_fromStoA1(i) * S * CL_fromStoA1(i) * 1e-1;
            WBL_fromStoA1_new(i)  = q_fromStoA1(i) * S * CL_fromStoA1_fullvehicle(i) * 1e-1;
            LHT_fromStoA1(i)  = (0.5)*(V_fromStoA1(i)^2)*(S)*(rho0)*(CLHT_fromStoA1(i))*(1e-1);
            CMCL_fromStoA1(i) = CMCL_new;
            CMCD_fromStoA1(i) = CMCD_new;
            CMCG_fromStoA1(i) = CMCG_new;
            
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1.value = CL_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1_new.value = CL_fromStoA1_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1.value = alfa_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1_rad.value = deg2rad(alfa_fromStoA1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1_rad.Attributes.unit = "rad";    
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1_new.value = alfa_fromStoA1_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1_new_rad.value = deg2rad(alfa_fromStoA1_fullvehicle);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1_new_rad.Attributes.unit = "rad";      
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromStoA1.value = CD_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromStoA1.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromStoA1.value = q_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromStoA1.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromStoA1.value = WBL_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromStoA1.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromStoA1_new.value = WBL_fromStoA1_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromStoA1_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromStoA1.value = CMCL_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromStoA1.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromStoA1.value = CMCD_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromStoA1.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromStoA1.value = CMCT_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromStoA1.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromStoA1.value = CMCG_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromStoA1.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromStoA1.value = CLHT_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromStoA1.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromStoA1.value = LHT_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromStoA1.Attributes.unit = "daN";        
        % =================================================================
        if max(n_gust_cruise_plus) > nmax
        % =================================================================    
            % CL CALCULATIONS - POSITIVE LOAD FACTORS
            % ------------------------------------------------------------------------- 
            % Calculation of the CL for the wing - body configuration. It must be
            % noticed that the calculations relative to lift coefficient in this
            % particular range of airspeed V and load factor n are going to give a
            % constant as a results, which is the maximum lift coefficient following
            % the lift curve of the aircrat. The previously defined CLMAX is a little
            % bit greater than the one obtained from V - n diagram flight condition to
            % take into account the extra lift produced by the horizontal empennage.
            % =============================================================
            % FROM A1 TO C1
            CL_fromA1toC1   = zeros(length(V_fromA1toC1), 1);
            alfa_fromA1toC1 = zeros(length(V_fromA1toC1), 1);
            CD_fromA1toC1   = zeros(length(V_fromA1toC1), 1);
            q_fromA1toC1    = zeros(length(V_fromA1toC1), 1);
            WBL_fromA1toC1  = zeros(length(V_fromA1toC1), 1);
            WBL_fromA1toC1_new  = zeros(length(V_fromA1toC1), 1);
            CMCL_fromA1toC1 = zeros(length(V_fromA1toC1), 1);
            CMCD_fromA1toC1 = zeros(length(V_fromA1toC1), 1);
            CMCT_fromA1toC1 = zeros(length(V_fromA1toC1), 1);
            CMCG_fromA1toC1 = zeros(length(V_fromA1toC1), 1); 
            CLHT_fromA1toC1 = zeros(length(V_fromA1toC1), 1);
            LHT_fromA1toC1  = zeros(length(V_fromA1toC1), 1);
            % =================================================================
            CL_fromA1toC1_fullvehicle   = zeros(length(V_fromA1toC1), 1);
            alfa_fromA1toC1_fullvehicle = zeros(length(V_fromA1toC1), 1);
            % =================================================================
            for i = 1:length(V_fromA1toC1)
                CL_fromA1toC1(i)   = CLmax_func(rho0, V_fromA1toC1(i), WS, n_fromA1toC1(i)); 
                if CL_fromA1toC1(i) < CL_max_fullmodel
                    CL_fromA1toC1(i) = CL_fromA1toC1(i);
                elseif CL_fromA1toC1(i) > CL_max_fullmodel
                    CL_fromA1toC1(i) = CL_max_fullmodel;
                end
                alfa_fromA1toC1(i) = alpha_fullmodel(CL_fromA1toC1(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
                CD_fromA1toC1(i)   = polyval(p_CD_wb1, alfa_fromA1toC1(i));
                q_fromA1toC1(i)    = 0.5*rho0*(V_fromA1toC1(i))^2;
                WBL_fromA1toC1(i)  = q_fromA1toC1(i)*S*CL_fromA1toC1(i)*1e-1;
                CMCL_fromA1toC1(i) = CLWB_contrib(obj1, CL_fromA1toC1(i), deg2rad(alfa_fromA1toC1(i)), XAC, XCG, bCG, MAC);
                CMCD_fromA1toC1(i) = CDWB_contrib(obj1, CL_fromA1toC1(i), deg2rad(alfa_fromA1toC1(i)), XAC, XCG, bCG, MAC);
                CMCT_fromA1toC1(i) = CT_contr(obj1, CD_fromA1toC1(i), Thrust_axes, MAC);
                CMCG_fromA1toC1(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromA1toC1(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromA1toC1(i) = CL_Tail(obj1, CMCL_fromA1toC1(i), CMCD_fromA1toC1(i), ...
                                                 CMCT_fromA1toC1(i), CMCG_fromA1toC1(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromA1toC1(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromA1toC1(i) = (0.5)*(V_fromA1toC1(i)^2)*(S)*(rho0)*(CLHT_fromA1toC1(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail           = CLHT_fromA1toC1(i);
                CL_wb             = CL_fromA1toC1(i);
                CL_new_fromA1toC1 = CL_wb - CL_tail;
                tol               = 1e-3; 
                n                 = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromA1toC1 = alpha_fullmodel(CL_new_fromA1toC1, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
                   CD_wb               = polyval(p_CD_wb1, alfa_new_fromA1toC1);
                   CMCL_new            = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromA1toC1), XAC, XCG, bCG, MAC);
                   CMCD_new            = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromA1toC1), XAC, XCG, bCG, MAC);
                   CMCT_new            = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new            = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new            = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromA1toC1));
                   CL_tail             = CLHT_new;
                   CL_new_fromA1toC1   = CL_fromA1toC1(i) - CL_tail;
                   CL_wb               = CL_wb + (CL_new_fromA1toC1 - CL_wb) * 1e-1;
                   n                   = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromA1toC1_fullvehicle(i)   = CL_new_fromA1toC1; 
%                 CL_fromA1toC1(i)   = CL_new_fromA1toC1; 
                CLHT_fromA1toC1(i) = CL_tail;
                alfa_fromA1toC1_fullvehicle(i) = alfa_new_fromA1toC1;
%                 alfa_fromA1toC1(i) = alfa_new_fromA1toC1;
                WBL_fromA1toC1(i)  = q_fromA1toC1(i) * S * CL_fromA1toC1(i) * 1e-1;
                WBL_fromA1toC1_new(i)  = q_fromA1toC1(i) * S * CL_fromA1toC1_fullvehicle(i) * 1e-1;
                LHT_fromA1toC1(i)  = (0.5)*(V_fromA1toC1(i)^2)*(S)*(rho0)*(CLHT_fromA1toC1(i))*(1e-1);
                CMCL_fromA1toC1(i) = CMCL_new;
                CMCD_fromA1toC1(i) = CMCD_new;
                CMCG_fromA1toC1(i) = CMCG_new;
                
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC1.value = CL_fromA1toC1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC1.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC1_new.value = CL_fromA1toC1_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC1_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC1.value = alfa_fromA1toC1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC1.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC1_rad.value = deg2rad(alfa_fromA1toC1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC1_rad.Attributes.unit = "rad";     
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC1_new.value = alfa_fromA1toC1_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC1_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC1_new_rad.value = deg2rad(alfa_fromA1toC1_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC1_new_rad.Attributes.unit = "rad";    
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromA1toC1.value = CD_fromA1toC1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromA1toC1.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromA1toC1.value = q_fromA1toC1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromA1toC1.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromA1toC1.value = WBL_fromA1toC1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromA1toC1.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC1_new.value = WBL_fromA1toC1_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC1_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromA1toC1.value = CMCL_fromA1toC1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromA1toC1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromA1toC1.value = CMCD_fromA1toC1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromA1toC1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromA1toC1.value = CMCT_fromA1toC1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromA1toC1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromA1toC1.value = CMCG_fromA1toC1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromA1toC1.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromA1toC1.value = CLHT_fromA1toC1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromA1toC1.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromA1toC1.value = LHT_fromA1toC1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromA1toC1.Attributes.unit = "daN";  
            % =============================================================
            % FROM C1 TO C
            CL_fromC1toC   = zeros(length(V_fromC1toC), 1);
            alfa_fromC1toC = zeros(length(V_fromC1toC), 1);
            CD_fromC1toC   = zeros(length(V_fromC1toC), 1);
            q_fromC1toC    = zeros(length(V_fromC1toC), 1);
            WBL_fromC1toC  = zeros(length(V_fromC1toC), 1);
            WBL_fromC1toC_new  = zeros(length(V_fromC1toC), 1);
            CMCL_fromC1toC = zeros(length(V_fromC1toC), 1);
            CMCD_fromC1toC = zeros(length(V_fromC1toC), 1);
            CMCT_fromC1toC = zeros(length(V_fromC1toC), 1);
            CMCG_fromC1toC = zeros(length(V_fromC1toC), 1); 
            CLHT_fromC1toC = zeros(length(V_fromC1toC), 1);
            LHT_fromC1toC  = zeros(length(V_fromC1toC), 1);
            % =================================================================
            CL_fromC1toC_fullvehicle   = zeros(length(V_fromC1toC), 1);
            alfa_fromC1toC_fullvehicle = zeros(length(V_fromC1toC), 1);
            % =================================================================
            for i = 1:length(V_fromC1toC)
                CL_fromC1toC(i)   = CLmax_func(rho0, V_fromC1toC(i), WS, n_fromC1toC(i));
                if CL_fromC1toC(i) < CL_max_fullmodel
                    CL_fromC1toC(i) = CL_fromC1toC(i);
                elseif CL_fromC1toC(i) > CL_max_fullmodel
                    CL_fromC1toC(i) = CL_max_fullmodel;
                end
                alfa_fromC1toC(i) = alpha_fullmodel(CL_fromC1toC(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
                CD_fromC1toC(i)   = polyval(p_CD_wb1, alfa_fromC1toC(i));
                q_fromC1toC(i)    = 0.5*rho0*(V_fromC1toC(i))^2;
                WBL_fromC1toC(i)  = q_fromC1toC(i)*S*CL_fromC1toC(i)*1e-1;
                CMCL_fromC1toC(i) = CLWB_contrib(obj1, CL_fromC1toC(i), deg2rad(alfa_fromC1toC(i)), XAC, XCG, bCG, MAC);
                CMCD_fromC1toC(i) = CDWB_contrib(obj1, CL_fromC1toC(i), deg2rad(alfa_fromC1toC(i)), XAC, XCG, bCG, MAC);
                CMCT_fromC1toC(i) = CT_contr(obj1, CD_fromC1toC(i), Thrust_axes, MAC);
                CMCG_fromC1toC(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromC1toC(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromC1toC(i) = CL_Tail(obj1, CMCL_fromC1toC(i), CMCD_fromC1toC(i), ...
                                                 CMCT_fromC1toC(i), CMCG_fromC1toC(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromC1toC(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromC1toC(i) = (0.5)*(V_fromC1toC(i)^2)*(S)*(rho0)*(CLHT_fromC1toC(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail           = CLHT_fromC1toC(i);
                CL_wb             = CL_fromC1toC(i);
                CL_new_fromC1toC  = CL_wb - CL_tail;
                tol               = 1e-3; 
                n                 = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromC1toC = alpha_fullmodel(CL_new_fromC1toC, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
                   CD_wb             = polyval(p_CD_wb1, alfa_fromC1toC(i));
                   CMCL_new          = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromC1toC), XAC, XCG, bCG, MAC);
                   CMCD_new          = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromC1toC), XAC, XCG, bCG, MAC);
                   CMCT_new          = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new          = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new          = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromC1toC));
                   CL_tail           = CLHT_new;
                   CL_new_fromC1toC  = CL_fromC1toC(i) - CL_tail;
                   CL_wb             = CL_wb + (CL_new_fromC1toC - CL_wb) * 1e-1;
                   n                 = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromC1toC_fullvehicle(i)   = CL_new_fromC1toC; 
                CLHT_fromC1toC(i) = CL_tail;
                alfa_fromC1toC_fullvehicle(i) = alfa_new_fromC1toC;
%                 alfa_fromC1toC(i) = alfa_new_fromC1toC;
                WBL_fromC1toC(i)  = q_fromC1toC(i) * S * CL_fromC1toC(i) * 1e-1;
                WBL_fromC1toC_new(i)  = q_fromC1toC(i) * S * CL_fromC1toC_fullvehicle(i) * 1e-1;
                LHT_fromC1toC(i)  = (0.5)*(V_fromC1toC(i)^2)*(S)*(rho0)*(CLHT_fromC1toC(i))*(1e-1);
                CMCL_fromC1toC(i) = CMCL_new;
                CMCD_fromC1toC(i) = CMCD_new;
                CMCG_fromC1toC(i) = CMCG_new;
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC1toC.value = CL_fromC1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC1toC.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC1toC_new.value = CL_fromC1toC_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC1toC_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC1toC.value = alfa_fromC1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC1toC.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC1toC_rad.value = deg2rad(alfa_fromC1toC);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC1toC_rad.Attributes.unit = "rad";    
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC1toC_new.value = alfa_fromC1toC_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC1toC_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC1toC_new_rad.value = deg2rad(alfa_fromC1toC_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC1toC_new_rad.Attributes.unit = "rad";       
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromC1toC.value = CD_fromC1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromC1toC.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromC1toC.value = q_fromC1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromC1toC.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromC1toC.value = WBL_fromC1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromC1toC.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromC1toC_new.value = WBL_fromC1toC_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromC1toC_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromC1toC.value = CMCL_fromC1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromC1toC.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromC1toC.value = CMCD_fromC1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromC1toC.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromC1toC.value = CMCT_fromC1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromC1toC.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromC1toC.value = CMCG_fromC1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromC1toC.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromC1toC.value = CLHT_fromC1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromC1toC.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromC1toC.value = LHT_fromC1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromC1toC.Attributes.unit = "daN";
            % =============================================================
            % FROM C TO C2
            CL_fromCtoC2   = zeros(length(V_fromCtoC2), 1);
            alfa_fromCtoC2 = zeros(length(V_fromCtoC2), 1);
            CD_fromCtoC2   = zeros(length(V_fromCtoC2), 1);
            q_fromCtoC2    = zeros(length(V_fromCtoC2), 1);
            WBL_fromCtoC2  = zeros(length(V_fromCtoC2), 1);
            WBL_fromCtoC2_new  = zeros(length(V_fromCtoC2), 1);
            CMCL_fromCtoC2 = zeros(length(V_fromCtoC2), 1);
            CMCD_fromCtoC2 = zeros(length(V_fromCtoC2), 1);
            CMCT_fromCtoC2 = zeros(length(V_fromCtoC2), 1);
            CMCG_fromCtoC2 = zeros(length(V_fromCtoC2), 1); 
            CLHT_fromCtoC2 = zeros(length(V_fromCtoC2), 1);
            LHT_fromCtoC2  = zeros(length(V_fromCtoC2), 1);
            % =================================================================
            CL_fromCtoC2_fullvehicle   = zeros(length(V_fromCtoC2), 1);
            alfa_fromCtoC2_fullvehicle = zeros(length(V_fromCtoC2), 1);
            % =================================================================
            for i = 1:length(V_fromCtoC2)
                CL_fromCtoC2(i)   = CLmax_func(rho0, V_fromCtoC2(i), WS, n_fromCtoC2(i));
                if CL_fromCtoC2(i) < CL_max_fullmodel
                    CL_fromCtoC2(i) = CL_fromCtoC2(i);
                elseif CL_fromCtoC2(i) > CL_max_fullmodel
                    CL_fromCtoC2(i) = CL_max_fullmodel;
                end
                alfa_fromCtoC2(i) = alpha_fullmodel(CL_fromCtoC2(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
                CD_fromCtoC2(i)   = polyval(p_CD_wb1, alfa_fromCtoC2(i));
                q_fromCtoC2(i)    = 0.5*rho0*(V_fromCtoC2(i))^2;
                WBL_fromCtoC2(i)  = q_fromCtoC2(i)*S*CL_fromCtoC2(i)*1e-1;
                CMCL_fromCtoC2(i) = CLWB_contrib(obj1, CL_fromCtoC2(i), deg2rad(alfa_fromCtoC2(i)), XAC, XCG, bCG, MAC);
                CMCD_fromCtoC2(i) = CDWB_contrib(obj1, CL_fromCtoC2(i), deg2rad(alfa_fromCtoC2(i)), XAC, XCG, bCG, MAC);
                CMCT_fromCtoC2(i) = CT_contr(obj1, CD_fromCtoC2(i), Thrust_axes, MAC);
                CMCG_fromCtoC2(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromCtoC2(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromCtoC2(i) = CL_Tail(obj1, CMCL_fromCtoC2(i), CMCD_fromCtoC2(i), ...
                                                 CMCT_fromCtoC2(i), CMCG_fromCtoC2(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromCtoC2(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromCtoC2(i) = (0.5)*(V_fromCtoC2(i)^2)*(S)*(rho0)*(CLHT_fromCtoC2(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail           = CLHT_fromCtoC2(i);
                CL_wb             = CL_fromCtoC2(i);
                CL_new_fromCtoC2  = CL_wb - CL_tail;
                tol               = 1e-3; 
                n                 = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromCtoC2 = alpha_fullmodel(CL_new_fromCtoC2, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
                   CD_wb              = polyval(p_CD_wb1, alfa_new_fromCtoC2);
                   CMCL_new           = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromCtoC2), XAC, XCG, bCG, MAC);
                   CMCD_new           = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromCtoC2), XAC, XCG, bCG, MAC);
                   CMCT_new           = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new           = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new           = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromCtoC2));
                   CL_tail            = CLHT_new;
                   CL_new_fromCtoC2   = CL_fromCtoC2(i) - CL_tail;
                   CL_wb              = CL_wb + (CL_new_fromCtoC2 - CL_wb) * 1e-1;
                   n                  = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromCtoC2_fullvehicle(i)   = CL_new_fromCtoC2; 
                CLHT_fromCtoC2(i) = CL_tail;
                alfa_fromCtoC2_fullvehicle(i) = alfa_new_fromCtoC2;
%                 alfa_fromCtoC2(i) = alfa_new_fromCtoC2;
                WBL_fromCtoC2(i)  = q_fromCtoC2(i) * S * CL_fromCtoC2(i) * 1e-1;
                WBL_fromCtoC2_new(i)  = q_fromCtoC2(i) * S * CL_fromCtoC2_fullvehicle(i) * 1e-1;
                LHT_fromCtoC2(i)  = (0.5)*(V_fromCtoC2(i)^2)*(S)*(rho0)*(CLHT_fromCtoC2(i))*(1e-1);
                CMCL_fromCtoC2(i) = CMCL_new;
                CMCD_fromCtoC2(i) = CMCD_new;
                CMCG_fromCtoC2(i) = CMCG_new;
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoC2.value = CL_fromCtoC2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoC2.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoC2_new.value = CL_fromCtoC2_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoC2_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoC2.value = alfa_fromCtoC2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoC2.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoC2_rad.value = deg2rad(alfa_fromCtoC2);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoC2_rad.Attributes.unit = "rad";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoC2_new.value = alfa_fromCtoC2_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoC2_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoC2_new_rad.value = deg2rad(alfa_fromCtoC2_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoC2_new_rad.Attributes.unit = "rad";       
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtoC2.value = CD_fromCtoC2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtoC2.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromCtoC2.value = q_fromCtoC2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromCtoC2.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtoC2.value = WBL_fromCtoC2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtoC2.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoC2_new.value = WBL_fromCtoC2_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoC2_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromCtoC2.value = CMCL_fromCtoC2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromCtoC2.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromCtoC2.value = CMCD_fromCtoC2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromCtoC2.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromCtoC2.value = CMCT_fromCtoC2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromCtoC2.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromCtoC2.value = CMCG_fromCtoC2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromCtoC2.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromCtoC2.value = CLHT_fromCtoC2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromCtoC2.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromCtoC2.value = LHT_fromCtoC2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromCtoC2.Attributes.unit = "daN";  
            % =============================================================
            % FROM C2 TO D
            CL_fromC2toD   = zeros(length(V_fromC2toD), 1);
            alfa_fromC2toD = zeros(length(V_fromC2toD), 1);
            CD_fromC2toD   = zeros(length(V_fromC2toD), 1);
            q_fromC2toD    = zeros(length(V_fromC2toD), 1);
            WBL_fromC2toD  = zeros(length(V_fromC2toD), 1);
            WBL_fromC2toD_new  = zeros(length(V_fromC2toD), 1);
            CMCL_fromC2toD = zeros(length(V_fromC2toD), 1);
            CMCD_fromC2toD = zeros(length(V_fromC2toD), 1);
            CMCT_fromC2toD = zeros(length(V_fromC2toD), 1);
            CMCG_fromC2toD = zeros(length(V_fromC2toD), 1); 
            CLHT_fromC2toD = zeros(length(V_fromC2toD), 1);
            LHT_fromC2toD  = zeros(length(V_fromC2toD), 1);
            % =================================================================
            CL_fromC2toD_fullvehicle   = zeros(length(V_fromC2toD), 1);
            alfa_fromC2toD_fullvehicle = zeros(length(V_fromC2toD), 1);
            % =================================================================
            for i = 1:length(V_fromC2toD)
                CL_fromC2toD(i)   = CLmax_func(rho0, V_fromC2toD(i), WS, n_fromC2toD(i));
                if CL_fromC2toD(i) < CL_max_fullmodel
                    CL_fromC2toD(i) = CL_fromC2toD(i);
                elseif CL_fromC2toD(i) > CL_max_fullmodel
                    CL_fromC2toD(i) = CL_max_fullmodel;
                end
                alfa_fromC2toD(i) = alpha_fullmodel(CL_fromC2toD(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
                CD_fromC2toD(i)   = polyval(p_CD_wb1, alfa_fromC2toD(i));
                q_fromC2toD(i)    = 0.5*rho0*(V_fromC2toD(i))^2;
                WBL_fromC2toD(i)  = q_fromC2toD(i)*S*CL_fromC2toD(i)*1e-1;
                CMCL_fromC2toD(i) = CLWB_contrib(obj1, CL_fromC2toD(i), deg2rad(alfa_fromC2toD(i)), XAC, XCG, bCG, MAC);
                CMCD_fromC2toD(i) = CDWB_contrib(obj1, CL_fromC2toD(i), deg2rad(alfa_fromC2toD(i)), XAC, XCG, bCG, MAC);
                CMCT_fromC2toD(i) = CT_contr(obj1, CD_fromC2toD(i), Thrust_axes, MAC);
                CMCG_fromC2toD(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromC2toD(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromC2toD(i) = CL_Tail(obj1, CMCL_fromC2toD(i), CMCD_fromC2toD(i), ...
                                                 CMCT_fromC2toD(i), CMCG_fromC2toD(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromC2toD(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromC2toD(i) = (0.5)*(V_fromC2toD(i)^2)*(S)*(rho0)*(CLHT_fromC2toD(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail           = CLHT_fromC2toD(i);
                CL_wb             = CL_fromC2toD(i);
                CL_new_fromC2toD  = CL_wb - CL_tail;
                tol               = 1e-3; 
                n                 = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromC2toD = alpha_fullmodel(CL_new_fromC2toD, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
                   CD_wb              = polyval(p_CD_wb1, alfa_new_fromC2toD);
                   CMCL_new           = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromC2toD), XAC, XCG, bCG, MAC);
                   CMCD_new           = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromC2toD), XAC, XCG, bCG, MAC);
                   CMCT_new           = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new           = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new           = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromC2toD));
                   CL_tail            = CLHT_new;
                   CL_new_fromC2toD   = CL_fromC2toD(i) - CL_tail;
                   CL_wb              = CL_wb + (CL_new_fromC2toD - CL_wb) * 1e-1;
                   n                  = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromC2toD_fullvehicle(i)   = CL_new_fromC2toD; 
                CLHT_fromC2toD(i) = CL_tail;
                alfa_fromC2toD_fullvehicle(i) = alfa_new_fromC2toD;
%                 alfa_fromC2toD(i) = alfa_new_fromC2toD;
                WBL_fromC2toD(i)  = q_fromC2toD(i) * S * CL_fromC2toD(i) * 1e-1;
                WBL_fromC2toD_new(i)  = q_fromC2toD(i) * S * CL_fromC2toD_fullvehicle(i) * 1e-1;
                LHT_fromC2toD(i)  = (0.5)*(V_fromC2toD(i)^2)*(S)*(rho0)*(CLHT_fromC2toD(i))*(1e-1);
                CMCL_fromC2toD(i) = CMCL_new;
                CMCD_fromC2toD(i) = CMCD_new;
                CMCG_fromC2toD(i) = CMCG_new;
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC2toD.value = CL_fromC2toD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC2toD.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC2toD_new.value = CL_fromC2toD_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC2toD_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC2toD.value = alfa_fromC2toD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC2toD.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC2toD_rad.value = deg2rad(alfa_fromC2toD);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC2toD_rad.Attributes.unit = "rad";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC2toD_new.value = alfa_fromC2toD_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC2toD_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC2toD_new_rad.value = deg2rad(alfa_fromC2toD_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC2toD_new_rad.Attributes.unit = "rad";       
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromC2toD.value = CD_fromC2toD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromC2toD.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromC2toD.value = q_fromC2toD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromC2toD.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromC2toD.value = WBL_fromC2toD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromC2toD.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromC2toD_new.value = WBL_fromC2toD_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromC2toD_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromC2toD.value = CMCL_fromC2toD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromC2toD.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromC2toD.value = CMCD_fromC2toD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromC2toD.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromC2toD.value = CMCT_fromC2toD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromC2toD.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromC2toD.value = CMCG_fromC2toD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromC2toD.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromC2toD.value = CLHT_fromC2toD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromC2toD.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromC2toD.value = LHT_fromC2toD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromC2toD.Attributes.unit = "daN";
            % =============================================================
            % FROM D TO 0
            CL_fromDto0   = zeros(length(V_fromDto0), 1);
            alfa_fromDto0 = zeros(length(V_fromDto0), 1);
            CD_fromDto0   = zeros(length(V_fromDto0), 1);
            q_fromDto0    = zeros(length(V_fromDto0), 1);
            WBL_fromDto0  = zeros(length(V_fromDto0), 1);
            WBL_fromDto0_new  = zeros(length(V_fromDto0), 1);
            CMCL_fromDto0 = zeros(length(V_fromDto0), 1);
            CMCD_fromDto0 = zeros(length(V_fromDto0), 1);
            CMCT_fromDto0 = zeros(length(V_fromDto0), 1);
            CMCG_fromDto0 = zeros(length(V_fromDto0), 1); 
            CLHT_fromDto0 = zeros(length(V_fromDto0), 1);
            LHT_fromDto0  = zeros(length(V_fromDto0), 1);
            % =================================================================
            CL_fromDto0_fullvehicle   = zeros(length(V_fromDto0), 1);
            alfa_fromDto0_fullvehicle = zeros(length(V_fromDto0), 1);
            % =================================================================
            for i = 1:length(V_fromDto0)
                CL_fromDto0(i)   = CLmax_func(rho0, V_fromDto0(i), WS, n_fromDto0(i));
                if CL_fromDto0(i) < CL_max_fullmodel
                    CL_fromDto0(i) = CL_fromDto0(i);
                elseif CL_fromDto0(i) > CL_max_fullmodel
                    CL_fromDto0(i) = CL_max_fullmodel;
                end
                alfa_fromDto0(i) = alpha_fullmodel(CL_fromDto0(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
                CD_fromDto0(i)   = polyval(p_CD_wb1, alfa_fromDto0(i));
                q_fromDto0(i)    = 0.5*rho0*(V_fromDto0(i))^2;
                WBL_fromDto0(i)  = q_fromDto0(i)*S*CL_fromDto0(i)*1e-1;
                CMCL_fromDto0(i) = CLWB_contrib(obj1, CL_fromDto0(i), deg2rad(alfa_fromDto0(i)), XAC, XCG, bCG, MAC);
                CMCD_fromDto0(i) = CDWB_contrib(obj1, CL_fromDto0(i), deg2rad(alfa_fromDto0(i)), XAC, XCG, bCG, MAC);
                CMCT_fromDto0(i) = CT_contr(obj1, CD_fromDto0(i), Thrust_axes, MAC);
                CMCG_fromDto0(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromDto0(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromDto0(i) = CL_Tail(obj1, CMCL_fromDto0(i), CMCD_fromDto0(i), ...
                                                 CMCT_fromDto0(i), CMCG_fromDto0(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromDto0(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromDto0(i) = (0.5)*(V_fromDto0(i)^2)*(S)*(rho0)*(CLHT_fromDto0(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail           = CLHT_fromDto0(i);
                CL_wb             = CL_fromDto0(i);
                CL_new_fromDto0  = CL_wb - CL_tail;
                tol               = 1e-3; 
                n                 = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromDto0  = alpha_fullmodel(CL_new_fromDto0, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
                   CD_wb              = polyval(p_CD_wb1, alfa_new_fromDto0);
                   CMCL_new           = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromDto0), XAC, XCG, bCG, MAC);
                   CMCD_new           = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromDto0), XAC, XCG, bCG, MAC);
                   CMCT_new           = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new           = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new           = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromDto0));
                   CL_tail            = CLHT_new;
                   CL_new_fromDto0    = CL_fromDto0(i) - CL_tail;
                   CL_wb              = CL_wb + (CL_new_fromDto0 - CL_wb) * 1e-1;
                   n                  = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromDto0_fullvehicle(i)   = CL_new_fromDto0; 
                CLHT_fromDto0(i) = CL_tail;
                alfa_fromDto0_fullvehicle(i) = alfa_new_fromDto0;
%                 alfa_fromDto0(i) = alfa_new_fromDto0;
                WBL_fromDto0(i)  = q_fromDto0(i) * S * CL_fromDto0(i) * 1e-1;
                WBL_fromDto0_new(i)  = q_fromDto0(i) * S * CL_fromDto0_fullvehicle(i) * 1e-1;
                LHT_fromDto0(i)  = (0.5)*(V_fromDto0(i)^2)*(S)*(rho0)*(CLHT_fromDto0(i))*(1e-1);
                CMCL_fromDto0(i) = CMCL_new;
                CMCD_fromDto0(i) = CMCD_new;
                CMCG_fromDto0(i) = CMCG_new;
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0.value = CL_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.value = CL_fromDto0_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0.value = alfa_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_rad.value = deg2rad(alfa_fromDto0);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_rad.Attributes.unit = "rad";     
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_new.value = alfa_fromDto0_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_new_rad.value = deg2rad(alfa_fromDto0_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_new_rad.Attributes.unit = "rad";     
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDto0.value = CD_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromDto0.value = q_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromDto0.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDto0.value = WBL_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDto0.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.value = WBL_fromDto0_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromDto0.value = CMCL_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromDto0.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromDto0.value = CMCD_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromDto0.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromDto0.value = CMCT_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromDto0.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromDto0.value = CMCG_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromDto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromDto0.value = CLHT_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromDto0.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromDto0.value = LHT_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromDto0.Attributes.unit = "daN";  
            % =============================================================
            % =============================================================
        % TAIL AIRLOADS AT UNIT LOAD - N = 1
        
        % AIRSPEED
        VS = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.value;
        VD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
        V_unit_load_factor = zeros(numb, 1);
        V_unit_load_factor(1) = VS;
        for i = 2:length(V_unit_load_factor)
            V_unit_load_factor(i) = V_unit_load_factor(i-1) + (VD - VS)*(1/length(V_unit_load_factor));
        end  
        
        % DYNAMIC PRESSURE
        q_unit_load_factor = zeros(length(V_unit_load_factor), 1);
        for i = 1:length(V_unit_load_factor)
            q_unit_load_factor(i) = 0.5*rho0*(V_unit_load_factor(i)^2);
            % q_unit_load_factor(i) = 0.5*rho0*(V_unit_load_factor(i)^2);
        end
        
        % WING BODY LIFT COEFFICIENT
        Aircraft_weight = ones(length(V_unit_load_factor), 1)*Mass*g;
        CLWB_unit_load_factor = zeros(length(V_unit_load_factor), 1);
        for i = 1:length(V_unit_load_factor)
            Aero_force = q_unit_load_factor(i)*S;
            CLWB_unit_load_factor(i) = ((Aircraft_weight(i))*(1/Aero_force));
        end
        
        % ALFA WING BODY UNIT LOAD
        alfanew_unit_load_factor = zeros(length(V_unit_load_factor), 1);
        for i = 1:length(V_unit_load_factor)
            alfanew_unit_load_factor(i) = alpha_fullmodel(CLWB_unit_load_factor(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
        end
        
        % DRAG LIFT COEFFICIENT
        CD_unit_load_factor = zeros(length(V_unit_load_factor), 1);
        for i = 1:length(V_unit_load_factor)
            CD_unit_load_factor(i) = polyval(p_CD_wb1, alfa_new_fromDto0);
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value                   = V_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.Attributes.unit         = "m/s";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.value                   = q_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.Attributes.unit         = "Pa";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value                = CLWB_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.Attributes.unit      = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.value               = alfanew_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.Attributes.unit     = "deg"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value           = deg2rad(alfanew_unit_load_factor);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.Attributes.unit = "rad"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value                  = CD_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.Attributes.unit        = "Non dimensional"; 
        
%         % Interpolation coefficient from actual aerodynamic data
%         CL_supp    = str2num(Aircraft.Certification.Aerodynamic_data.CL.value);
%         alpha_supp = str2num(Aircraft.Certification.Aerodynamic_data.alpha.value);
%         x          = 0.027*ones(length(CL_supp), 1);
%         p = polyfit(CL_supp + x, alpha_supp, 2);
%         CLalfa_deg = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;
%         
%         % ALPHA CALCULATION 
%         alpha_unit_load_factor = zeros(length(V_unit_load_factor), 1);
%         for i = 1:length(alpha_unit_load_factor)
%             alpha_unit_load_factor(i) = alpha_calc(obj1, ...
%                                                    CLWB_unit_load_factor(i), CL0, CL_star, CLalfa_deg, p);
%         end
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.value = alpha_unit_load_factor;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.Attributes.unit = "deg";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value = deg2rad(alpha_unit_load_factor);
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.Attributes.unit = "rad"; 
        
        % WING BODY LIFT 
        WBLift_unit_load_factor = zeros(length(alfanew_unit_load_factor), 1);
        for i = 1:length(alfanew_unit_load_factor)   
            WBLift_unit_load_factor(i) = q_unit_load_factor(i) * S * (CLWB_unit_load_factor(i))*(1E-1);  
            % WBLift_unit_load_factor(i) = (0.5) * (V_unit_load_factor(i)^2) * rho0 * S * (CLWB_unit_load_factor(i))*(1E-1);    
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.value           = WBLift_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.Attributes.unit = "daN";    
        
        % -----------------------------------------------------------------
        alpha_unit_load_factor = alfanew_unit_load_factor;
        % EVALUATION OF TAIL LOADS CONTRIBUTION     
        CMCL_unit_load_factor             = zeros(length(alpha_unit_load_factor), 1);
        CMCD_unit_load_factor             = zeros(length(alpha_unit_load_factor), 1);
        CMCT_unit_load_factor             = zeros(length(alpha_unit_load_factor), 1);
        CMCG_unit_load_factor             = zeros(length(alpha_unit_load_factor), 1);
        CLHT_unit_load_factor             = zeros(length(alpha_unit_load_factor), 1);
        LHT_unit_load_factor              = zeros(length(alpha_unit_load_factor), 1);
        LW_unit_load_factor               = zeros(length(alpha_unit_load_factor), 1);
        CL_unit_load_factor_fullvehicle   = zeros(length(alpha_unit_load_factor), 1);
        alfa_unit_load_factor_fullvehicle = zeros(length(alpha_unit_load_factor), 1);
        LW_unit_load_factor_new           = zeros(length(alpha_unit_load_factor), 1);
        for i = 1:length(alpha_unit_load_factor)
            CMCL_unit_load_factor(i) = CLWB_contrib(obj1, CLWB_unit_load_factor(i), ...
                           deg2rad(alpha_unit_load_factor(i)), XAC, XCG, bCG, MAC);
            CMCD_unit_load_factor(i) = CDWB_contrib(obj1, CLWB_unit_load_factor(i), ...
                           deg2rad(alpha_unit_load_factor(i)), XAC, XCG, bCG, MAC);
            CMCT_unit_load_factor(i) = CT_contr(obj1, CD_unit_load_factor(i), Thrust_axes, MAC);      
            CMCG_unit_load_factor(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CMCG_unit_load_factor(i));   
            CLHT_unit_load_factor(i) = CL_Tail(obj1, CMCL_unit_load_factor(i), CMCD_unit_load_factor(i), ...
                                               CMCT_unit_load_factor(i), CMCG_unit_load_factor(i), ...
                                               l_ht, MAC, XAC, XCG, deg2rad(alpha_unit_load_factor(i)));
            LHT_unit_load_factor(i) = q_unit_load_factor(i) * S * (CLHT_unit_load_factor(i))*(1e-1);  
            
            % DRAFT VERSION OF ITERATION 
            CL_tail                 = CLHT_unit_load_factor(i);
            CL_wb                   = CLWB_unit_load_factor(i);
            CL_new_unit_load_factor = CL_wb - CL_tail;
            if CL_new_unit_load_factor > CL_max
                CL_new_unit_load_factor = CL_max;
            end
            tol                     = 1e-3; 
            n                       = 1;
            while abs(CL_wb - CL_tail) > tol
               alfa_new_unit_load_factor = alpha_fullmodel(CL_new_unit_load_factor, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
               CD_wb                     = polyval(p_CD_wb1, alfa_new_unit_load_factor);
               CMCL_new                  = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_unit_load_factor), XAC, XCG, bCG, MAC);
               CMCD_new                  = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_unit_load_factor), XAC, XCG, bCG, MAC);
               CMCT_new                  = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
               CMCG_new                  = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
               CLHT_new                  = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_new_unit_load_factor));
               CL_tail                   = CLHT_new;
               CL_new_unit_load_factor   = CLWB_unit_load_factor(i) - CL_tail;
               CL_wb                     = CL_wb + (CL_new_unit_load_factor - CL_wb) * 1e-1;
               n                         = n + 1;
               if n == 15
                   break
               end
            end
            CL_unit_load_factor_fullvehicle(i)   = CL_new_unit_load_factor; 
            CLHT_unit_load_factor(i)             = CL_tail;
            alfa_unit_load_factor_fullvehicle(i) = alfa_new_unit_load_factor;
%             alfa_fromA1toC(i) = alfa_new_fromA1toC;
            LW_unit_load_factor(i)     = q_unit_load_factor(i) * S * CLWB_unit_load_factor(i) * 1e-1;
            LW_unit_load_factor_new(i) = q_unit_load_factor(i) * S * CL_unit_load_factor_fullvehicle(i) * 1e-1;
            LHT_unit_load_factor(i)    = q_unit_load_factor(i) * S * (CLHT_unit_load_factor(i))*(1e-1);
            CMCL_unit_load_factor(i)   = CMCL_new;
            CMCD_unit_load_factor(i)   = CMCD_new;
            CMCG_unit_load_factor(i)   = CMCG_new;                         
            % LHT_unit_load_factor(i) = (0.5)*(V_unit_load_factor(i)^2) * S * rho0 * (CLHT_unit_load_factor(i))*(1e-1);
            % LHT_unit_load_factor(i) = (0.5)*(V_unit_load_factor(i)^2) * S_ht * rho0 * (CLHT_unit_load_factor(i))*(1e-1);
            % LW_unit_load_factor(i)  = WBLift_unit_load_factor(i) - LHT_unit_load_factor(i);                                                                                     
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_unit_load_factor.value = CMCL_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_unit_load_factor.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_unit_load_factor.value = CMCD_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_unit_load_factor.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_unit_load_factor.value = CMCT_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_unit_load_factor.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_unit_load_factor.value = CMCG_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_unit_load_factor.Attributes.unit = "Non dimensional";     
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_unit_load_factor.value = CLHT_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_unit_load_factor.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_unit_load_factor.value = LHT_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_unit_load_factor.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_unit_load_factor.value = LW_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_unit_load_factor.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_unit_load_factor_new.value = LW_unit_load_factor_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_unit_load_factor_new.Attributes.unit = "daN";
            % =================================================================================================================
            % TAIL LOADS DIAGRAM - POSITIVE SIDE   
            HT_balancing_loads = figure(5); 
            hold on; grid on; grid minor; 

    %         ylim([LHT_from0toS(1) LHT_fromDto0(1)])
    %         xlim([0.0 V_fromDto0(1)])

            plot(V_from0toS,  LHT_from0toS,    '-r', 'LineWidth', 1)
            plot(V_fromStoA1, LHT_fromStoA1,   '-r', 'LineWidth', 1)
            plot(V_fromA1toC1, LHT_fromA1toC1, '-r', 'LineWidth', 1)
            plot(V_fromC1toC, LHT_fromC1toC,   '-r', 'LineWidth', 1)
            plot(V_fromCtoC2, LHT_fromCtoC2,   '-r', 'LineWidth', 1)
            plot(V_fromC2toD,  LHT_fromC2toD,  '-r', 'LineWidth', 1)
            plot(V_fromDto0,  LHT_fromDto0,    '-r', 'LineWidth', 1)
            
            % UNIT LOAD DISTR.
            plot(V_unit_load_factor, LHT_unit_load_factor,  '--k', 'LineWidth', 2)
            % ---------------------------------------------------------------------
            plot(V_from0toS(1),    LHT_from0toS(1),      'k.', 'MarkerSize', 10)
            plot(V_from0toS(end),  LHT_from0toS(end),    'k.', 'MarkerSize', 10)
            plot(V_fromStoA1(end), LHT_fromStoA1(end),   'k.', 'MarkerSize', 10)
            plot(V_fromA1toC1(end), LHT_fromA1toC1(end), 'k.', 'MarkerSize', 10)
            plot(V_fromC1toC(end), LHT_fromC1toC(end),   'k.', 'MarkerSize', 10)
            plot(V_fromCtoC2(end), LHT_fromCtoC2(end),   'k.', 'MarkerSize', 10)
            plot(V_fromC2toD(end),  LHT_fromC2toD(end),  'k.', 'MarkerSize', 10)
            plot(V_fromDto0(end),  LHT_fromDto0(end),    'k.', 'MarkerSize', 10)
            % ---------------------------------------------------------------------
            text(V_fromStoA1(1),   LHT_fromStoA1(1),     '  S',  'FontSize', 6)
            text(V_fromStoA1(end), LHT_fromStoA1(end),   '  A1', 'FontSize', 6)
            text(V_fromA1toC1(end), LHT_fromA1toC1(end), '  C1', 'FontSize', 6)
            text(V_fromC1toC(end), LHT_fromC1toC(end),   '  C',  'FontSize', 6)
            text(V_fromCtoC2(end), LHT_fromCtoC2(end),     '  C2', 'FontSize', 6)
            text(V_fromC2toD(end), LHT_fromC2toD(end),     '  D',  'FontSize', 6)
            text(40.75, -18, 'n = 1')
            % ---------------------------------------------------------------------
            xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
            ylabel("Horizontal tail lift - $L_{ht}$ (daN)", "Interpreter", "latex")
            title("Horizontal empennage airloads per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct      
%              
%             % MAIN WING LOADS DIAGRAM 
%             CL_from0toS_new   = CL_from0toS   + CLHT_from0toS;        
%             CL_fromStoA1_new  = CL_fromStoA1  + CLHT_fromStoA1;        
%             CL_fromA1toC1_new = CL_fromA1toC1 + CLHT_fromA1toC1;        
%             CL_fromC1toC_new  = CL_fromC1toC  + CLHT_fromC1toC;      
%             CL_fromCtoC2_new  = CL_fromCtoC2  + CLHT_fromCtoC2;
%             CL_fromC2toD_new  = CL_fromC2toD  + CLHT_fromC2toD;
%             CL_fromDto0_new   = CL_fromDto0   + CLHT_fromDto0;
%             % MAIN WING LOADS DIAGRAM        
%             CL_from0toS_new   = CL_from0toS   - CLHT_from0toS;        
%             CL_fromStoA1_new  = CL_fromStoA1  - CLHT_fromStoA1;        
%             CL_fromA1toC1_new = CL_fromA1toC1 - CLHT_fromA1toC1;        
%             CL_fromC1toC_new  = CL_fromC1toC  - CLHT_fromC1toC;      
%             CL_fromCtoC2_new  = CL_fromCtoC2  - CLHT_fromCtoC2;
%             CL_fromC2toD_new  = CL_fromC2toD  - CLHT_fromC2toD;
%             CL_fromDto0_new   = CL_fromDto0   - CLHT_fromDto0;
%             LW_from0toS_new   = zeros(length(CL_from0toS_new), 1);
%             LW_fromStoA1_new  = zeros(length(CL_from0toS_new), 1);
%             LW_fromA1toC1_new = zeros(length(CL_from0toS_new), 1);
%             LW_fromC1toC_new  = zeros(length(CL_from0toS_new), 1);
%             LW_fromCtoC2_new  = zeros(length(CL_from0toS_new), 1);
%             LW_fromC2toD_new  = zeros(length(CL_from0toS_new), 1);
%             LW_fromDto0_new   = zeros(length(CL_from0toS_new), 1);
%              for i = 1:length(CL_from0toS_new)
%                 LW_from0toS_new(i)   = (0.5)*(V_from0toS(i)^2)* rho0 * S *(CL_from0toS_new(i))*(1e-1);  
%                 LW_fromStoA1_new(i)  = (0.5)*(V_fromStoA1(i)^2)* rho0 * S *(CL_fromStoA1_new(i))*(1e-1);
%                 LW_fromA1toC1_new(i) = (0.5)*(V_fromA1toC1(i)^2)* rho0 * S *(CL_fromA1toC1_new(i))*(1e-1);
%                 LW_fromC1toC_new(i)  = (0.5)*(V_fromC1toC(i)^2)* rho0 * S *(CL_fromC1toC_new(i))*(1e-1);
%                 LW_fromCtoC2_new(i)  = (0.5)*(V_fromCtoC2(i)^2)* rho0 * S *(CL_fromCtoC2_new(i))*(1e-1);
%                 LW_fromC2toD_new(i)  = (0.5)*(V_fromC2toD(i)^2)* rho0 * S *(CL_fromC2toD_new(i))*(1e-1);
%                 LW_fromDto0_new(i)   = (0.5)*(V_fromDto0(i)^2)* rho0 * S *(CL_fromDto0_new(i))*(1e-1);
% %                 LW_from0toS_new(i)   = (0.5)*(V_from0toS(i)^2)* rho0 * S *(CL_from0toS(i))*(1e-1);  
% %                 LW_fromStoA1_new(i)  = (0.5)*(V_fromStoA1(i)^2)* rho0 * S *(CL_fromStoA1(i))*(1e-1);
% %                 LW_fromA1toC1_new(i) = (0.5)*(V_fromA1toC1(i)^2)* rho0 * S *(CL_fromA1toC1(i))*(1e-1);
% %                 LW_fromC1toC_new(i)  = (0.5)*(V_fromC1toC(i)^2)* rho0 * S *(CL_fromC1toC(i))*(1e-1);
% %                 LW_fromCtoC2_new(i)  = (0.5)*(V_fromCtoC2(i)^2)* rho0 * S *(CL_fromCtoC2(i))*(1e-1);
% %                 LW_fromC2toD_new(i)  = (0.5)*(V_fromC2toD(i)^2)* rho0 * S *(CL_fromC2toD(i))*(1e-1);
% %                 LW_fromDto0_new(i)   = (0.5)*(V_fromDto0(i)^2)* rho0 * S *(CL_fromDto0(i))*(1e-1);
%             end       
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS_new.value = CL_from0toS_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1_new.value = CL_fromStoA1_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC1_new.value = CL_fromA1toC1_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC1_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC1toC_new.value = CL_fromC1toC_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC1toC_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoC2_new.value = CL_fromCtoC2_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoC2_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC2toD_new.value = CL_fromC2toD_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC2toD_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.value = CL_fromDto0_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.Attributes.unit = "Non dimensional";
%             % ------------------------------------------------------------------------------------------------------------------------
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.value = LW_from0toS_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromStoA1_new.value = LW_fromStoA1_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromStoA1_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC1_new.value = LW_fromA1toC1_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC1_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromC1toC_new.value = LW_fromC1toC_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromC1toC_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoC2_new.value = LW_fromCtoC2_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoC2_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromC2toD_new.value = LW_fromC2toD_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromC2toD_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.value = LW_fromDto0_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.Attributes.unit = "daN";  
            
            % TAIL LOADS DIAGRAM - POSITIVE SIDE   
            Wing_balancing_loads = figure(6); 
            hold on; grid on; grid minor; 

    %         ylim([LHT_from0toS(1) LHT_fromDto0(1)])
    %         xlim([0.0 V_fromDto0(1)])
%             % --------------------------------------------------------------------
%             plot(V_from0toS,   LW_from0toS_new,   '-r', 'LineWidth', 1)
%             plot(V_fromStoA1,  LW_fromStoA1_new,  '-r', 'LineWidth', 1)
%             plot(V_fromA1toC1, LW_fromA1toC1_new, '-r', 'LineWidth', 1)
%             plot(V_fromC1toC,  LW_fromC1toC_new,  '-r', 'LineWidth', 1)
%             plot(V_fromCtoC2,  LW_fromCtoC2_new,  '-r', 'LineWidth', 1)
%             plot(V_fromC2toD,  LW_fromC2toD_new,  '-r', 'LineWidth', 1)
%             plot(V_fromDto0,   LW_fromDto0_new,   '-r', 'LineWidth', 1)
%             % UNIT LOAD DISTR.
%             plot(V_unit_load_factor, LW_unit_load_factor,  '--k', 'LineWidth', 2)
%             % ---------------------------------------------------------------------
%             plot(V_from0toS(1),     LW_from0toS_new(1),     'k.', 'MarkerSize', 10)
%             plot(V_from0toS(end),   LW_from0toS_new(end),   'k.', 'MarkerSize', 10)
%             plot(V_fromStoA1(end),  LW_fromStoA1_new(end),  'k.', 'MarkerSize', 10)
%             plot(V_fromA1toC1(end), LW_fromA1toC1_new(end), 'k.', 'MarkerSize', 10)
%             plot(V_fromC1toC(end),  LW_fromC1toC_new(end),  'k.', 'MarkerSize', 10)
%             plot(V_fromCtoC2(end),  LW_fromCtoC2_new(end),  'k.', 'MarkerSize', 10)
%             plot(V_fromC2toD(end),  LW_fromC2toD_new(end),  'k.', 'MarkerSize', 10)
%             plot(V_fromDto0(end),   LW_fromDto0_new(end),   'k.', 'MarkerSize', 10)
%             % ---------------------------------------------------------------------
%             text(V_fromStoA1(1),   LW_fromStoA1_new(1),     '  S',  'FontSize', 6)
%             text(V_fromStoA1(end), LW_fromStoA1_new(end),   '  A1', 'FontSize', 6)
%             text(V_fromA1toC1(end), LW_fromA1toC1_new(end), '  C1', 'FontSize', 6)
%             text(V_fromC1toC(end),  LW_fromC1toC_new(end),  '  C',  'FontSize', 6)
%             text(V_fromCtoC2(end),  LW_fromCtoC2_new(end),  '  C2', 'FontSize', 6)
%             text(V_fromC2toD(end),  LW_fromC2toD_new(end),  '  D',  'FontSize', 6)
%             % text(40.75, -18, 'n = 1')
%             % ---------------------------------------------------------------------
            % --------------------------------------------------------------------
            plot(V_from0toS,   WBL_from0toS,   '-r', 'LineWidth', 1)
            plot(V_fromStoA1,  WBL_fromStoA1,  '-r', 'LineWidth', 1)
            plot(V_fromA1toC1, WBL_fromA1toC1, '-r', 'LineWidth', 1)
            plot(V_fromC1toC,  WBL_fromC1toC,  '-r', 'LineWidth', 1)
            plot(V_fromCtoC2,  WBL_fromCtoC2,  '-r', 'LineWidth', 1)
            plot(V_fromC2toD,  WBL_fromC2toD,  '-r', 'LineWidth', 1)
            plot(V_fromDto0,   WBL_fromDto0,   '-r', 'LineWidth', 1)
            % UNIT LOAD DISTR.
            plot(V_unit_load_factor, LW_unit_load_factor,  '--k', 'LineWidth', 2)
            % ---------------------------------------------------------------------
            plot(V_from0toS(1),     WBL_from0toS(1),     'k.', 'MarkerSize', 10)
            plot(V_from0toS(end),   WBL_from0toS(end),   'k.', 'MarkerSize', 10)
            plot(V_fromStoA1(end),  WBL_fromStoA1(end),  'k.', 'MarkerSize', 10)
            plot(V_fromA1toC1(end), WBL_fromA1toC1(end), 'k.', 'MarkerSize', 10)
            plot(V_fromC1toC(end),  WBL_fromC1toC(end),  'k.', 'MarkerSize', 10)
            plot(V_fromCtoC2(end),  WBL_fromCtoC2(end),  'k.', 'MarkerSize', 10)
            plot(V_fromC2toD(end),  WBL_fromC2toD(end),  'k.', 'MarkerSize', 10)
            plot(V_fromDto0(end),   WBL_fromDto0(end),   'k.', 'MarkerSize', 10)
            % ---------------------------------------------------------------------
            text(V_fromStoA1(1),    WBL_fromStoA1(1),       '  S',  'FontSize', 6)
            text(V_fromStoA1(end),  WBL_fromStoA1(end),     '  A1', 'FontSize', 6)
            text(V_fromA1toC1(end), WBL_fromA1toC1(end),    '  C1', 'FontSize', 6)
            text(V_fromC1toC(end),  WBL_fromC1toC(end),  '  C',  'FontSize', 6)
            text(V_fromCtoC2(end),  WBL_fromCtoC2(end),  '  C2', 'FontSize', 6)
            text(V_fromC2toD(end),  WBL_fromC2toD(end),  '  D',  'FontSize', 6)
            % text(40.75, -18, 'n = 1')
            % ---------------------------------------------------------------------
            xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
            ylabel("Wing lift - $L_{w}$ (daN)", "Interpreter", "latex")
            title("Wing airloads per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct 
            
            % STORE INSIDE THE AIRCRAFT STRUCTURE VARIABLE

            % POINT S
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.value = q_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS.value = WBL_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS_new.value = WBL_from0toS_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTS.value = LHT_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTS.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.value = alfa_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_rad.value = deg2rad(alfa_from0toS(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_rad.Attributes.unit = "rad";
            alfaS = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_new.value = alfa_from0toS_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_new_rad.value = deg2rad(alfa_from0toS_fullvehicle(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_new_rad.Attributes.unit = "rad";
            alfaS_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_new.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value = CL_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S_new.value = CL_from0toS_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CLHT_S.value = CLHT_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CLHT_S.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CD_S.value = CD_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CD_S.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CM_S.value = CMCG_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CM_S.Attributes.unit = "Non dimensional";
            % POINT A1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.value = q_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1.value = WBL_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1_new.value = WBL_fromStoA1_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.value = LHT_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.value = alfa_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_rad.value = deg2rad(alfa_fromStoA1(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_rad.Attributes.unit = "rad";
            alfaA1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_new.value = alfa_fromStoA1_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_new_rad.value = deg2rad(alfa_fromStoA1_fullvehicle(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_new_rad.Attributes.unit = "rad";
            alfaA1_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_new.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.value = CL_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1_new.value = CL_fromStoA1_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CLHT_A1.value = CLHT_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CLHT_A1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CD_A1.value = CD_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CD_A1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CM_A1.value = CMCG_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CM_A1.Attributes.unit = "Non dimensional";
            % POINT C1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.qC1.value = q_fromA1toC1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.qC1.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.LC1.value = WBL_fromA1toC1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.LC1.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.LC1_new.value = WBL_fromA1toC1_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.LC1_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.LHTC1.value = LHT_fromA1toC1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.LHTC1.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.alfaC1.value = alfa_fromA1toC1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.alfaC1.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.alfaC1_rad.value = deg2rad(alfa_fromA1toC1(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.alfaC1_rad.Attributes.unit = "rad"; 
            alfaC1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.alfaC1.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.alfaC1_new.value = alfa_fromA1toC1_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.alfaC1_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.alfaC1_new_rad.value = deg2rad(alfa_fromA1toC1_fullvehicle(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.alfaC1_new_rad.Attributes.unit = "rad"; 
            alfaC1_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.alfaC1_new.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CL_C1.value = CL_fromA1toC1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CL_C1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CL_C1_new.value = CL_fromA1toC1_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CL_C1_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CLHT_C1.value = CLHT_fromA1toC1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CLHT_C1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CD_C1.value = CD_fromA1toC1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CD_C1.Attributes.unit = "Non dimensional";     
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CM_C1.value = CMCG_fromA1toC1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CM_C1.Attributes.unit = "Non dimensional";    
            % POINT C
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value = q_fromC1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.value = WBL_fromC1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.value = WBL_fromC1toC_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.value = LHT_fromC1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value = alfa_fromC1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_rad.value = deg2rad(alfa_fromC1toC(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_rad.Attributes.unit = "rad";
            alfaC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_new.value = alfa_fromC1toC_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_new_rad.value = deg2rad(alfa_fromC1toC_fullvehicle(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_new_rad.Attributes.unit = "rad";
            alfaC_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_new.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value =  CL_fromC1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C_new.value = CL_fromC1toC_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CLHT_C.value = CLHT_fromC1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CLHT_C.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CD_C.value = CD_fromC1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CD_C.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CM_C.value = CMCG_fromC1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CM_C.Attributes.unit = "Non dimensional";              
            % POINT C2
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.qC2.value = q_fromCtoC2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.qC2.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.LC2.value = WBL_fromCtoC2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.LC2.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.LC2_new.value = WBL_fromCtoC2_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.LC2_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.LHTC2.value = LHT_fromCtoC2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.LHTC2.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.alfaC2.value = alfa_fromCtoC2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.alfaC2.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.alfaC2_rad.value = deg2rad(alfa_fromCtoC2(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.alfaC2_rad.Attributes.unit = "rad"; 
            alfaC2 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.alfaC2.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.alfaC2_new.value = alfa_fromCtoC2_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.alfaC2_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.alfaC2_new_rad.value = deg2rad(alfa_fromCtoC2_fullvehicle(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.alfaC2_new_rad.Attributes.unit = "rad"; 
            alfaC2_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.alfaC2_new.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CL_C2.value = CL_fromCtoC2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CL_C2.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CL_C2_new.value = CL_fromCtoC2_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CL_C2_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CLHT_C2.value = CLHT_fromCtoC2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CLHT_C2.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CD_C2.value = CD_fromCtoC2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CD_C2.Attributes.unit = "Non dimensional";  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CM_C2.value = CMCG_fromCtoC2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CM_C2.Attributes.unit = "Non dimensional";
            % POINT D
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value = q_fromC2toD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD.value = WBL_fromC2toD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.value = WBL_fromC2toD_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value = LHT_fromC2toD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value = alfa_fromC2toD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_rad.value = deg2rad(alfa_fromC2toD(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_rad.Attributes.unit = "rad";  
            alfaD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_new.value = alfa_fromC2toD_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_new_rad.value = deg2rad(alfa_fromC2toD_fullvehicle(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_new_rad.Attributes.unit = "rad";  
            alfaD_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_new.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value = CL_fromC2toD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D_new.value = CL_fromC2toD_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CLHT_D.value = CLHT_fromC2toD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CLHT_D.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_D.value = CD_fromC2toD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_D.Attributes.unit = "Non dimensional";  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CM_D.value = CMCG_fromC2toD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CM_D.Attributes.unit = "Non dimensional";    
        % =================================================================    
        elseif max(n_gust_cruise_plus) < nmax
        % =================================================================    
            % CL CALCULATIONS - POSITIVE LOAD FACTORS
            % ------------------------------------------------------------------------- 
            % Calculation of the CL for the wing - body configuration. It must be
            % noticed that the calculations relative to lift coefficient in this
            % particular range of airspeed V and load factor n are going to give a
            % constant as a results, which is the maximum lift coefficient following
            % the lift curve of the aircrat. The previously defined CLMAX is a little
            % bit greater than the one obtained from V - n diagram flight condition to
            % take into account the extra lift produced by the horizontal empennage.            
            % =============================================================
            % FROM A1 TO C
            CL_fromA1toC   = zeros(length(V_fromA1toC), 1);
            alfa_fromA1toC = zeros(length(V_fromA1toC), 1);
            CD_fromA1toC   = zeros(length(V_fromA1toC), 1);
            q_fromA1toC    = zeros(length(V_fromA1toC), 1);
            WBL_fromA1toC  = zeros(length(V_fromA1toC), 1);
            WBL_fromA1toC_new  = zeros(length(V_fromA1toC), 1);
            CMCL_fromA1toC = zeros(length(V_fromA1toC), 1);
            CMCD_fromA1toC = zeros(length(V_fromA1toC), 1);
            CMCT_fromA1toC = zeros(length(V_fromA1toC), 1);
            CMCG_fromA1toC = zeros(length(V_fromA1toC), 1); 
            CLHT_fromA1toC = zeros(length(V_fromA1toC), 1);
            LHT_fromA1toC  = zeros(length(V_fromA1toC), 1);
            % =================================================================
            CL_fromA1toC_fullvehicle   = zeros(length(V_fromA1toC), 1);
            alfa_fromA1toC_fullvehicle = zeros(length(V_fromA1toC), 1);
            % =================================================================
            for i = 1:length(V_fromA1toC)           
                CL_fromA1toC(i)   = CLmax_func(rho0, V_fromA1toC(i), WS, n_fromA1toC(i));
                if CL_fromA1toC(i) < CL_max_fullmodel
                    CL_fromA1toC(i) = CL_fromA1toC(i);
                elseif CL_fromA1toC(i) > CL_max_fullmodel
                    CL_fromA1toC(i) = CL_max_fullmodel;
                end
                alfa_fromA1toC(i) = alpha_fullmodel(CL_fromA1toC(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
                CD_fromA1toC(i)   = polyval(p_CD_wb1, alfa_fromA1toC(i));
                q_fromA1toC(i)    = 0.5*rho0*(V_fromA1toC(i))^2;
                WBL_fromA1toC(i)  = q_fromA1toC(i)*S*CL_fromA1toC(i)*1e-1;
                CMCL_fromA1toC(i) = CLWB_contrib(obj1, CL_fromA1toC(i), deg2rad(alfa_fromA1toC(i)), XAC, XCG, bCG, MAC);
                CMCD_fromA1toC(i) = CDWB_contrib(obj1, CL_fromA1toC(i), deg2rad(alfa_fromA1toC(i)), XAC, XCG, bCG, MAC);
                CMCT_fromA1toC(i) = CT_contr(obj1, CD_fromA1toC(i), Thrust_axes, MAC);
                CMCG_fromA1toC(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromA1toC(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromA1toC(i) = CL_Tail(obj1, CMCL_fromA1toC(i), CMCD_fromA1toC(i), ...
                                                 CMCT_fromA1toC(i), CMCG_fromA1toC(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromA1toC(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromA1toC(i) = (0.5)*(V_fromA1toC(i)^2)*(S)*(rho0)*(CLHT_fromA1toC(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail           = CLHT_fromA1toC(i);
                CL_wb             = CL_fromA1toC(i);
                CL_new_fromA1toC  = CL_wb - CL_tail;
                tol               = 1e-3; 
                n                 = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromA1toC = alpha_fullmodel(CL_new_fromA1toC, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
                   CD_wb              = polyval(p_CD_wb1, alfa_new_fromA1toC);
                   CMCL_new           = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromA1toC), XAC, XCG, bCG, MAC);
                   CMCD_new           = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromA1toC), XAC, XCG, bCG, MAC);
                   CMCT_new           = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new           = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new           = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromA1toC));
                   CL_tail            = CLHT_new;
                   CL_new_fromA1toC   = CL_fromA1toC(i) - CL_tail;
                   CL_wb              = CL_wb + (CL_new_fromA1toC - CL_wb) * 1e-1;
                   n                  = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromA1toC_fullvehicle(i)   = CL_new_fromA1toC; 
                CLHT_fromA1toC(i) = CL_tail;
                alfa_fromA1toC_fullvehicle(i) = alfa_new_fromA1toC;
%                 alfa_fromA1toC(i) = alfa_new_fromA1toC;
                WBL_fromA1toC(i)  = q_fromA1toC(i) * S * CL_fromA1toC(i) * 1e-1;
                WBL_fromA1toC_new(i)  = q_fromA1toC(i) * S * CL_fromA1toC_fullvehicle(i) * 1e-1;
                LHT_fromA1toC(i)  = (0.5)*(V_fromA1toC(i)^2)*(S)*(rho0)*(CLHT_fromA1toC(i))*(1e-1);
                CMCL_fromA1toC(i) = CMCL_new;
                CMCD_fromA1toC(i) = CMCD_new;
                CMCG_fromA1toC(i) = CMCG_new;
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC.value = CL_fromA1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC_new.value = CL_fromA1toC_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC.value = alfa_fromA1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC_rad.value = deg2rad(alfa_fromA1toC);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC_rad.Attributes.unit = "rad";    
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC_new.value = alfa_fromA1toC_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC_new_rad.value = deg2rad(alfa_fromA1toC_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC_new_rad.Attributes.unit = "rad";     
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromA1toC.value = CD_fromA1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromA1toC.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromA1toC.value = q_fromA1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromA1toC.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromA1toC.value = WBL_fromA1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromA1toC.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC_new.value = WBL_fromA1toC_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromA1toC.value = CMCL_fromA1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromA1toC.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromA1toC.value = CMCD_fromA1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromA1toC.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromA1toC.value = CMCT_fromA1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromA1toC.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromA1toC.value = CMCG_fromA1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromA1toC.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromA1toC.value = CLHT_fromA1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromA1toC.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromA1toC.value = LHT_fromA1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromA1toC.Attributes.unit = "daN";
            % =============================================================
            % FROM C TO D
            CL_fromCtoD   = zeros(length(V_fromCtoD), 1);
            alfa_fromCtoD = zeros(length(V_fromCtoD), 1);
            CD_fromCtoD   = zeros(length(V_fromCtoD), 1);
            q_fromCtoD    = zeros(length(V_fromCtoD), 1);
            WBL_fromCtoD  = zeros(length(V_fromCtoD), 1);
            WBL_fromCtoD_new  = zeros(length(V_fromCtoD), 1);
            CMCL_fromCtoD = zeros(length(V_fromCtoD), 1);
            CMCD_fromCtoD = zeros(length(V_fromCtoD), 1);
            CMCT_fromCtoD = zeros(length(V_fromCtoD), 1);
            CMCG_fromCtoD = zeros(length(V_fromCtoD), 1); 
            CLHT_fromCtoD = zeros(length(V_fromCtoD), 1);
            LHT_fromCtoD  = zeros(length(V_fromCtoD), 1);
            % =================================================================
            CL_fromCtoD_fullvehicle   = zeros(length(V_fromCtoD), 1);
            alfa_fromCtoD_fullvehicle = zeros(length(V_fromCtoD), 1);
            % =================================================================
            for i = 1:length(V_fromCtoD)
                CL_fromCtoD(i)   = CLmax_func(rho0, V_fromCtoD(i), WS, n_fromCtoD(i));
                if CL_fromCtoD(i) < CL_max_fullmodel
                    CL_fromCtoD(i) = CL_fromCtoD(i);
                elseif CL_fromCtoD(i) > CL_max_fullmodel
                    CL_fromCtoD(i) = CL_max_fullmodel;
                end
                alfa_fromCtoD(i) = alpha_fullmodel(CL_fromCtoD(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
                CD_fromCtoD(i)   = polyval(p_CD_wb1, alfa_fromCtoD(i));
                q_fromCtoD(i)    = 0.5*rho0*(V_fromCtoD(i))^2;
                WBL_fromCtoD(i)  = q_fromCtoD(i)*S*CL_fromCtoD(i)*1e-1;
                CMCL_fromCtoD(i) = CLWB_contrib(obj1, CL_fromCtoD(i), deg2rad(alfa_fromCtoD(i)), XAC, XCG, bCG, MAC);
                CMCD_fromCtoD(i) = CDWB_contrib(obj1, CL_fromCtoD(i), deg2rad(alfa_fromCtoD(i)), XAC, XCG, bCG, MAC);
                CMCT_fromCtoD(i) = CT_contr(obj1, CD_fromCtoD(i), Thrust_axes, MAC);
                CMCG_fromCtoD(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromCtoD(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromCtoD(i) = CL_Tail(obj1, CMCL_fromCtoD(i), CMCD_fromCtoD(i), ...
                                                 CMCT_fromCtoD(i), CMCG_fromCtoD(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromCtoD(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromCtoD(i) = (0.5)*(V_fromCtoD(i)^2)*(S)*(rho0)*(CLHT_fromCtoD(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail           = CLHT_fromCtoD(i);
                CL_wb             = CL_fromCtoD(i);
                CL_new_fromCtoD  = CL_wb - CL_tail;
                tol               = 1e-3; 
                n                 = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromCtoD  = alpha_fullmodel(CL_new_fromCtoD, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
                   CD_wb              = polyval(p_CD_wb1, alfa_new_fromCtoD);
                   CMCL_new           = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromCtoD), XAC, XCG, bCG, MAC);
                   CMCD_new           = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromCtoD), XAC, XCG, bCG, MAC);
                   CMCT_new           = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new           = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new           = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromCtoD));
                   CL_tail            = CLHT_new;
                   CL_new_fromCtoD    = CL_fromCtoD(i) - CL_tail;
                   CL_wb              = CL_wb + (CL_new_fromCtoD - CL_wb) * 1e-1;
                   n                  = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromCtoD_fullvehicle(i)   = CL_new_fromCtoD; 
                CLHT_fromCtoD(i) = CL_tail;
                alfa_fromCtoD_fullvehicle(i) = alfa_new_fromCtoD;
%                 alfa_fromCtoD(i) = alfa_new_fromCtoD;
                WBL_fromCtoD(i)  = q_fromCtoD(i) * S * CL_fromCtoD(i) * 1e-1;
                WBL_fromCtoD_new(i)  = q_fromCtoD(i) * S * CL_fromCtoD_fullvehicle(i) * 1e-1;
                LHT_fromCtoD(i)  = (0.5)*(V_fromCtoD(i)^2)*(S)*(rho0)*(CLHT_fromCtoD(i))*(1e-1);
                CMCL_fromCtoD(i) = CMCL_new;
                CMCD_fromCtoD(i) = CMCD_new;
                CMCG_fromCtoD(i) = CMCG_new;
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoD.value = CL_fromCtoD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoD.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoD_new.value = CL_fromCtoD_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoD_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoD.value = alfa_fromCtoD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoD.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoD_rad.value = deg2rad(alfa_fromCtoD);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoD_rad.Attributes.unit = "rad";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoD_new.value = alfa_fromCtoD_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoD_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoD_new_rad.value = deg2rad(alfa_fromCtoD_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoD_new_rad.Attributes.unit = "rad";       
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtoD.value = CD_fromCtoD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtoD.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromCtoD.value = q_fromCtoD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromCtoD.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtoD.value = WBL_fromCtoD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtoD.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoD_new.value = WBL_fromCtoD_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoD_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromCtoD.value = CMCL_fromCtoD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromCtoD.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromCtoD.value = CMCD_fromCtoD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromCtoD.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromCtoD.value = CMCT_fromCtoD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromCtoD.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromCtoD.value = CMCG_fromCtoD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromCtoD.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromCtoD.value = CLHT_fromCtoD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromCtoD.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromCtoD.value = LHT_fromCtoD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromCtoD.Attributes.unit = "daN";
            % =============================================================
            % FROM D TO 0
            
            CL_fromDto0   = zeros(length(V_fromDto0), 1);
            alfa_fromDto0 = zeros(length(V_fromDto0), 1);
            CD_fromDto0   = zeros(length(V_fromDto0), 1);
            q_fromDto0    = zeros(length(V_fromDto0), 1);
            WBL_fromDto0  = zeros(length(V_fromDto0), 1);
            WBL_fromDto0_new  = zeros(length(V_fromDto0), 1);
            CMCL_fromDto0 = zeros(length(V_fromDto0), 1);
            CMCD_fromDto0 = zeros(length(V_fromDto0), 1);
            CMCT_fromDto0 = zeros(length(V_fromDto0), 1);
            CMCG_fromDto0 = zeros(length(V_fromDto0), 1); 
            CLHT_fromDto0 = zeros(length(V_fromDto0), 1);
            LHT_fromDto0  = zeros(length(V_fromDto0), 1);
            % =================================================================
            CL_fromDto0_fullvehicle   = zeros(length(V_fromDto0), 1);
            alfa_fromDto0_fullvehicle = zeros(length(V_fromDto0), 1);
            % =================================================================
            for i = 1:length(V_fromDto0) 
                CL_fromDto0(i)   = CLmax_func(rho0, V_fromDto0(i), WS, n_fromDto0(i));
                if CL_fromDto0(i) < CL_max_fullmodel
                    CL_fromDto0(i) = CL_fromDto0(i);
                elseif CL_fromDto0(i) > CL_max_fullmodel
                    CL_fromDto0(i) = CL_max_fullmodel;
                end
                alfa_fromDto0(i) = alpha_fullmodel(CL_fromDto0(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
                CD_fromDto0(i)   = polyval(p_CD_wb1, alfa_fromDto0(i));
                q_fromDto0(i)    = 0.5*rho0*(V_fromDto0(i))^2;
                WBL_fromDto0(i)  = q_fromDto0(i)*S*CL_fromDto0(i)*1e-1;
                CMCL_fromDto0(i) = CLWB_contrib(obj1, CL_fromDto0(i), deg2rad(alfa_fromDto0(i)), XAC, XCG, bCG, MAC);
                CMCD_fromDto0(i) = CDWB_contrib(obj1, CL_fromDto0(i), deg2rad(alfa_fromDto0(i)), XAC, XCG, bCG, MAC);
                CMCT_fromDto0(i) = CT_contr(obj1, CD_fromDto0(i), Thrust_axes, MAC);
                CMCG_fromDto0(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromDto0(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromDto0(i) = CL_Tail(obj1, CMCL_fromDto0(i), CMCD_fromDto0(i), ...
                                                 CMCT_fromDto0(i), CMCG_fromDto0(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromDto0(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromDto0(i) = (0.5)*(V_fromDto0(i)^2)*(S)*(rho0)*(CLHT_fromDto0(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail           = CLHT_fromDto0(i);
                CL_wb             = CL_fromDto0(i);
                CL_new_fromDto0  = CL_wb - CL_tail;
                tol               = 1e-3; 
                n                 = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromDto0  = alpha_fullmodel(CL_new_fromDto0, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
                   CD_wb              = polyval(p_CD_wb1, alfa_new_fromDto0);
                   CMCL_new           = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromDto0), XAC, XCG, bCG, MAC);
                   CMCD_new           = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromDto0), XAC, XCG, bCG, MAC);
                   CMCT_new           = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new           = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new           = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromDto0));
                   CL_tail            = CLHT_new;
                   CL_new_fromDto0    = CL_fromDto0(i) - CL_tail;
                   CL_wb              = CL_wb + (CL_new_fromDto0 - CL_wb) * 1e-1;
                   n                  = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromDto0_fullvehicle(i)   = CL_new_fromDto0; 
                CLHT_fromDto0(i) = CL_tail;
                alfa_fromDto0_fullvehicle(i) = alfa_new_fromDto0;
%                 alfa_fromDto0(i) = alfa_new_fromDto0;
                WBL_fromDto0(i)  = q_fromDto0(i) * S * CL_fromDto0(i) * 1e-1;
                WBL_fromDto0_new(i)  = q_fromDto0(i) * S * CL_fromDto0_fullvehicle(i) * 1e-1;
                LHT_fromDto0(i)  = (0.5)*(V_fromDto0(i)^2)*(S)*(rho0)*(CLHT_fromDto0(i))*(1e-1);
                CMCL_fromDto0(i) = CMCL_new;
                CMCD_fromDto0(i) = CMCD_new;
                CMCG_fromDto0(i) = CMCG_new;
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0.value = CL_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.value = CL_fromDto0_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0.value = alfa_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_rad.value = deg2rad(alfa_fromDto0);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_rad.Attributes.unit = "rad";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_new.value = alfa_fromDto0_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_new_rad.value = deg2rad(alfa_fromDto0_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_new_rad.Attributes.unit = "rad";         
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDto0.value = CD_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromDto0.value = q_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromDto0.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDto0.value = WBL_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDto0.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.value = WBL_fromDto0_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromDto0.value = CMCL_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromDto0.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromDto0.value = CMCD_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromDto0.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromDto0.value = CMCT_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromDto0.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromDto0.value = CMCG_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromDto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromDto0.value = CLHT_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromDto0.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromDto0.value = LHT_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromDto0.Attributes.unit = "daN";
            % =============================================================
            % =============================================================
        % TAIL AIRLOADS AT UNIT LOAD - N = 1
        
        % AIRSPEED
        VS = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.value;
        VD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
        V_unit_load_factor = zeros(numb, 1);
        V_unit_load_factor(1) = VS;
        for i = 2:length(V_unit_load_factor)
            V_unit_load_factor(i) = V_unit_load_factor(i-1) + (VD - VS)*(1/length(V_unit_load_factor));
        end  
        
        % DYNAMIC PRESSURE
        q_unit_load_factor = zeros(length(V_unit_load_factor), 1);
        for i = 1:length(V_unit_load_factor)
            q_unit_load_factor(i) = 0.5*rho0*(V_unit_load_factor(i)^2);
            % q_unit_load_factor(i) = 0.5*rho0*(V_unit_load_factor(i)^2);
        end
        
        % WING BODY LIFT COEFFICIENT
        Aircraft_weight = ones(length(V_unit_load_factor), 1)*Mass*g;
        CLWB_unit_load_factor = zeros(length(V_unit_load_factor), 1);
        for i = 1:length(V_unit_load_factor)
            Aero_force = q_unit_load_factor(i)*S;
            CLWB_unit_load_factor(i) = ((Aircraft_weight(i))*(1/Aero_force));
        end
        
        % ALFA WING BODY UNIT LOAD
        alfanew_unit_load_factor = zeros(length(V_unit_load_factor), 1);
        for i = 1:length(V_unit_load_factor)
            alfanew_unit_load_factor(i) = alpha_fullmodel(CLWB_unit_load_factor(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
        end
        
        % DRAG LIFT COEFFICIENT
        CD_unit_load_factor = zeros(length(V_unit_load_factor), 1);
        for i = 1:length(V_unit_load_factor)
            CD_unit_load_factor(i) = polyval(p_CD_wb1, alfa_new_fromDto0);
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value                   = V_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.Attributes.unit         = "m/s";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.value                   = q_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.Attributes.unit         = "Pa";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value                = CLWB_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.Attributes.unit      = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.value               = alfanew_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.Attributes.unit     = "deg"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value           = deg2rad(alfanew_unit_load_factor);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.Attributes.unit = "rad"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value                  = CD_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.Attributes.unit        = "Non dimensional"; 
        
%         % Interpolation coefficient from actual aerodynamic data
%         CL_supp    = str2num(Aircraft.Certification.Aerodynamic_data.CL.value);
%         alpha_supp = str2num(Aircraft.Certification.Aerodynamic_data.alpha.value);
%         x          = 0.027*ones(length(CL_supp), 1);
%         p = polyfit(CL_supp + x, alpha_supp, 2);
%         CLalfa_deg = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;
%         
%         % ALPHA CALCULATION 
%         alpha_unit_load_factor = zeros(length(V_unit_load_factor), 1);
%         for i = 1:length(alpha_unit_load_factor)
%             alpha_unit_load_factor(i) = alpha_calc(obj1, ...
%                                                    CLWB_unit_load_factor(i), CL0, CL_star, CLalfa_deg, p);
%         end
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.value = alpha_unit_load_factor;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.Attributes.unit = "deg";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value = deg2rad(alpha_unit_load_factor);
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.Attributes.unit = "rad"; 
        
        % WING BODY LIFT 
        WBLift_unit_load_factor = zeros(length(alfanew_unit_load_factor), 1);
        for i = 1:length(alfanew_unit_load_factor)   
            WBLift_unit_load_factor(i) = q_unit_load_factor(i) * S * (CLWB_unit_load_factor(i))*(1E-1);  
            % WBLift_unit_load_factor(i) = (0.5) * (V_unit_load_factor(i)^2) * rho0 * S * (CLWB_unit_load_factor(i))*(1E-1);    
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.value           = WBLift_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.Attributes.unit = "daN";    
        
        % -----------------------------------------------------------------
        alpha_unit_load_factor = alfanew_unit_load_factor;
        % EVALUATION OF TAIL LOADS CONTRIBUTION     
        CMCL_unit_load_factor             = zeros(length(alpha_unit_load_factor), 1);
        CMCD_unit_load_factor             = zeros(length(alpha_unit_load_factor), 1);
        CMCT_unit_load_factor             = zeros(length(alpha_unit_load_factor), 1);
        CMCG_unit_load_factor             = zeros(length(alpha_unit_load_factor), 1);
        CLHT_unit_load_factor             = zeros(length(alpha_unit_load_factor), 1);
        LHT_unit_load_factor              = zeros(length(alpha_unit_load_factor), 1);
        LW_unit_load_factor               = zeros(length(alpha_unit_load_factor), 1);
        CL_unit_load_factor_fullvehicle   = zeros(length(alpha_unit_load_factor), 1);
        alfa_unit_load_factor_fullvehicle = zeros(length(alpha_unit_load_factor), 1);
        LW_unit_load_factor_new           = zeros(length(alpha_unit_load_factor), 1);
        for i = 1:length(alpha_unit_load_factor)
            CMCL_unit_load_factor(i) = CLWB_contrib(obj1, CLWB_unit_load_factor(i), ...
                           deg2rad(alpha_unit_load_factor(i)), XAC, XCG, bCG, MAC);
            CMCD_unit_load_factor(i) = CDWB_contrib(obj1, CLWB_unit_load_factor(i), ...
                           deg2rad(alpha_unit_load_factor(i)), XAC, XCG, bCG, MAC);
            CMCT_unit_load_factor(i) = CT_contr(obj1, CD_unit_load_factor(i), Thrust_axes, MAC);      
            CMCG_unit_load_factor(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CMCG_unit_load_factor(i));   
            CLHT_unit_load_factor(i) = CL_Tail(obj1, CMCL_unit_load_factor(i), CMCD_unit_load_factor(i), ...
                                               CMCT_unit_load_factor(i), CMCG_unit_load_factor(i), ...
                                               l_ht, MAC, XAC, XCG, deg2rad(alpha_unit_load_factor(i)));
            LHT_unit_load_factor(i) = q_unit_load_factor(i) * S * (CLHT_unit_load_factor(i))*(1e-1);  
            
            % DRAFT VERSION OF ITERATION 
            CL_tail                 = CLHT_unit_load_factor(i);
            CL_wb                   = CLWB_unit_load_factor(i);
            CL_new_unit_load_factor = CL_wb - CL_tail;
            if CL_new_unit_load_factor > CL_max
                CL_new_unit_load_factor = CL_max;
            end
            tol                     = 1e-3; 
            n                       = 1;
            while abs(CL_wb - CL_tail) > tol
               alfa_new_unit_load_factor = alpha_fullmodel(CL_new_unit_load_factor, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
               CD_wb                     = polyval(p_CD_wb1, alfa_new_unit_load_factor);
               CMCL_new                  = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_unit_load_factor), XAC, XCG, bCG, MAC);
               CMCD_new                  = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_unit_load_factor), XAC, XCG, bCG, MAC);
               CMCT_new                  = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
               CMCG_new                  = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
               CLHT_new                  = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_new_unit_load_factor));
               CL_tail                   = CLHT_new;
               CL_new_unit_load_factor   = CLWB_unit_load_factor(i) - CL_tail;
               CL_wb                     = CL_wb + (CL_new_unit_load_factor - CL_wb) * 1e-1;
               n                         = n + 1;
               if n == 15
                   break
               end
            end
            CL_unit_load_factor_fullvehicle(i)   = CL_new_unit_load_factor; 
            CLHT_unit_load_factor(i)             = CL_tail;
            alfa_unit_load_factor_fullvehicle(i) = alfa_new_unit_load_factor;
%             alfa_fromA1toC(i) = alfa_new_fromA1toC;
            LW_unit_load_factor(i)     = q_unit_load_factor(i) * S * CLWB_unit_load_factor(i) * 1e-1;
            LW_unit_load_factor_new(i) = q_unit_load_factor(i) * S * CL_unit_load_factor_fullvehicle(i) * 1e-1;
            LHT_unit_load_factor(i)    = q_unit_load_factor(i) * S * (CLHT_unit_load_factor(i))*(1e-1);
            CMCL_unit_load_factor(i)   = CMCL_new;
            CMCD_unit_load_factor(i)   = CMCD_new;
            CMCG_unit_load_factor(i)   = CMCG_new;                         
            % LHT_unit_load_factor(i) = (0.5)*(V_unit_load_factor(i)^2) * S * rho0 * (CLHT_unit_load_factor(i))*(1e-1);
            % LHT_unit_load_factor(i) = (0.5)*(V_unit_load_factor(i)^2) * S_ht * rho0 * (CLHT_unit_load_factor(i))*(1e-1);
            % LW_unit_load_factor(i)  = WBLift_unit_load_factor(i) - LHT_unit_load_factor(i);                                                                                     
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_unit_load_factor.value = CMCL_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_unit_load_factor.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_unit_load_factor.value = CMCD_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_unit_load_factor.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_unit_load_factor.value = CMCT_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_unit_load_factor.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_unit_load_factor.value = CMCG_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_unit_load_factor.Attributes.unit = "Non dimensional";     
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_unit_load_factor.value = CLHT_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_unit_load_factor.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_unit_load_factor.value = LHT_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_unit_load_factor.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_unit_load_factor.value = LW_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_unit_load_factor.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_unit_load_factor_new.value = LW_unit_load_factor_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_unit_load_factor_new.Attributes.unit = "daN";          
            
            % =================================================================================================================
            % TAIL LOADS DIAGRAM - POSITIVE SIDE   
            HT_balancing_loads = figure(5); 
            hold on; grid on; grid minor; 

    %         ylim([LHT_from0toS(1) LHT_fromDto0(1)])
    %         xlim([0.0 V_fromDto0(1)])

            plot(V_from0toS,  LHT_from0toS,  '-r', 'LineWidth', 1)
            plot(V_fromStoA1, LHT_fromStoA1, '-r', 'LineWidth', 1)
            plot(V_fromA1toC, LHT_fromA1toC, '-r', 'LineWidth', 1)
            plot(V_fromCtoD,  LHT_fromCtoD,  '-r', 'LineWidth', 1)
            plot(V_fromDto0,  LHT_fromDto0,  '-r', 'LineWidth', 1)
            
            % UNIT LOAD DISTR.
            plot(V_unit_load_factor, LHT_unit_load_factor,  '--k', 'LineWidth', 2)
            % ---------------------------------------------------------------------
            plot(V_from0toS(1),    LHT_from0toS(1),    'k.', 'MarkerSize', 10)
            plot(V_from0toS(end),  LHT_from0toS(end),  'k.', 'MarkerSize', 10)
            plot(V_fromStoA1(end), LHT_fromStoA1(end), 'k.', 'MarkerSize', 10)
            plot(V_fromA1toC(end), LHT_fromA1toC(end), 'k.', 'MarkerSize', 10)
            plot(V_fromCtoD(end),  LHT_fromCtoD(end),  'k.', 'MarkerSize', 10)
            plot(V_fromDto0(end),  LHT_fromDto0(end),  'k.', 'MarkerSize', 10)
            % ---------------------------------------------------------------------
            text(V_fromStoA1(1),   LHT_fromStoA1(1),     '  S',  'FontSize', 6)
            text(V_fromStoA1(end), LHT_fromStoA1(end),   '  A1', 'FontSize', 6)
            text(V_fromA1toC(end), LHT_fromA1toC(end), '  C', 'FontSize', 6)
            text(V_fromCtoD(end),  LHT_fromCtoD(end),   '  D',  'FontSize', 6)
            text(40.75, -18, 'n = 1')
            % ---------------------------------------------------------------------
            xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
            ylabel("Horizontal tail lift - $L_{ht}$ (daN)", "Interpreter", "latex")
            title("Horizontal empennage airloads per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct   
            
%             % MAIN WING LOADS DIAGRAM        
%             CL_from0toS_new  = CL_from0toS  + CLHT_from0toS;        
%             CL_fromStoA1_new = CL_fromStoA1 + CLHT_fromStoA1;        
%             CL_fromA1toC_new = CL_fromA1toC + CLHT_fromA1toC;        
%             CL_fromCtoD_new  = CL_fromCtoD  + CLHT_fromCtoD;      
%             CL_fromDto0_new  = CL_fromDto0  + CLHT_fromDto0;
%             % MAIN WING LOADS DIAGRAM        
%             CL_from0toS_new  = CL_from0toS  - CLHT_from0toS;        
%             CL_fromStoA1_new = CL_fromStoA1 - CLHT_fromStoA1;        
%             CL_fromA1toC_new = CL_fromA1toC - CLHT_fromA1toC;        
%             CL_fromCtoD_new  = CL_fromCtoD  - CLHT_fromCtoD;      
%             CL_fromDto0_new  = CL_fromDto0  - CLHT_fromDto0;
%             LW_from0toS_new  = zeros(length(CL_from0toS_new), 1);
%             LW_fromStoA1_new = zeros(length(CL_from0toS_new), 1);
%             LW_fromA1toC_new = zeros(length(CL_from0toS_new), 1);
%             LW_fromCtoD_new  = zeros(length(CL_from0toS_new), 1);
%             LW_fromDto0_new  = zeros(length(CL_from0toS_new), 1);
%              for i = 1:length(CL_from0toS_new)
%                 LW_from0toS_new(i)  = (0.5)*(V_from0toS(i)^2)* rho0 * S *(CL_from0toS_new(i))*(1e-1);  
%                 LW_fromStoA1_new(i) = (0.5)*(V_fromStoA1(i)^2)* rho0 * S *(CL_fromStoA1_new(i))*(1e-1);
%                 LW_fromA1toC_new(i) = (0.5)*(V_fromA1toC(i)^2)* rho0 * S *(CL_fromA1toC_new(i))*(1e-1);
%                 LW_fromCtoD_new(i) = (0.5)*(V_fromCtoD(i)^2)* rho0 * S *(CL_fromCtoD_new(i))*(1e-1);
%                 LW_fromDto0_new(i)  = (0.5)*(V_fromDto0(i)^2)* rho0 * S *(CL_fromDto0_new(i))*(1e-1);
% 
% %                 LW_from0toS_new(i)  = (0.5)*(V_from0toS(i)^2)* rho0 * S *(CL_from0toS(i))*(1e-1);  
% %                 LW_fromStoA1_new(i) = (0.5)*(V_fromStoA1(i)^2)* rho0 * S *(CL_fromStoA1(i))*(1e-1);
% %                 LW_fromA1toC_new(i) = (0.5)*(V_fromA1toC(i)^2)* rho0 * S *(CL_fromA1toC(i))*(1e-1);
% %                 LW_fromCtoD_new(i) = (0.5)*(V_fromCtoD(i)^2)* rho0 * S *(CL_fromCtoD(i))*(1e-1);
% %                 LW_fromDto0_new(i)  = (0.5)*(V_fromDto0(i)^2)* rho0 * S *(CL_fromDto0(i))*(1e-1);
% 
%             end       
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS_new.value = CL_from0toS_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1_new.value = CL_fromStoA1_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC_new.value = CL_fromA1toC_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoD_new.value = CL_fromCtoD_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoD_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.value = CL_fromDto0_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.Attributes.unit = "Non dimensional";
%             % ------------------------------------------------------------------------------------------------------------------------
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.value = LW_from0toS_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromStoA1_new.value = LW_fromStoA1_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromStoA1_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC_new.value = LW_fromA1toC_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoD_new.value = LW_fromCtoD_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoD_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.value = LW_fromDto0_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.Attributes.unit = "daN";    
            
            % TAIL LOADS DIAGRAM - POSITIVE SIDE   
            Wing_balancing_loads = figure(6); 
            hold on; grid on; grid minor; 

    %         ylim([LHT_from0toS(1) LHT_fromDto0(1)])
    %         xlim([0.0 V_fromDto0(1)])
% 
%             plot(V_from0toS,  LW_from0toS_new,  '-r', 'LineWidth', 1)
%             plot(V_fromStoA1, LW_fromStoA1_new, '-r', 'LineWidth', 1)
%             plot(V_fromA1toC, LW_fromA1toC_new, '-r', 'LineWidth', 1)
%             plot(V_fromCtoD,  LW_fromCtoD_new, '-r', 'LineWidth', 1)
%             plot(V_fromDto0,  LW_fromDto0_new,  '-r', 'LineWidth', 1)
%             % UNIT LOAD DISTR.
%             plot(V_unit_load_factor, LW_unit_load_factor,  '--k', 'LineWidth', 2)
%             % ---------------------------------------------------------------------
%             plot(V_from0toS(1),    LW_from0toS_new(1),    'k.', 'MarkerSize', 10)
%             plot(V_from0toS(end),  LW_from0toS_new(end),  'k.', 'MarkerSize', 10)
%             plot(V_fromStoA1(end), LW_fromStoA1_new(end), 'k.', 'MarkerSize', 10)
%             plot(V_fromA1toC(end), LW_fromA1toC_new(end), 'k.', 'MarkerSize', 10)
%             plot(V_fromCtoD(end),  LW_fromCtoD_new(end),  'k.', 'MarkerSize', 10)
%             plot(V_fromDto0(end),  LW_fromDto0_new(end),  'k.', 'MarkerSize', 10)
%             % ---------------------------------------------------------------------
%             text(V_fromStoA1(1),   LW_fromStoA1_new(1),   '  S',  'FontSize', 6)
%             text(V_fromStoA1(end), LW_fromStoA1_new(end), '  A1', 'FontSize', 6)
%             text(V_fromA1toC(end), LW_fromA1toC_new(end), '  C',  'FontSize', 6)
%             text(V_fromCtoD(end),  LW_fromCtoD_new(end),  '  D', 'FontSize', 6)
%             % text(40.75, -18, 'n = 1')
            % ---------------------------------------------------------------------
            % --------------------------------------------------------------------
            plot(V_from0toS,   WBL_from0toS,  '-r', 'LineWidth', 1)
            plot(V_fromStoA1,  WBL_fromStoA1, '-r', 'LineWidth', 1)
            plot(V_fromA1toC,  WBL_fromA1toC, '-r', 'LineWidth', 1)
            plot(V_fromCtoD,   WBL_fromCtoD,  '-r', 'LineWidth', 1)
            plot(V_fromDto0,   WBL_fromDto0,  '-r', 'LineWidth', 1)
            % UNIT LOAD DISTR.
            plot(V_unit_load_factor, LW_unit_load_factor,  '--k', 'LineWidth', 2)
            % ---------------------------------------------------------------------
            plot(V_from0toS(1),     WBL_from0toS(1),     'k.', 'MarkerSize', 10)
            plot(V_from0toS(end),   WBL_from0toS(end),   'k.', 'MarkerSize', 10)
            plot(V_fromStoA1(end),  WBL_fromStoA1(end),  'k.', 'MarkerSize', 10)
            plot(V_fromA1toC(end),  WBL_fromA1toC(end), 'k.', 'MarkerSize', 10)
            plot(V_fromCtoD(end),   WBL_fromCtoD(end),  'k.', 'MarkerSize', 10)
            plot(V_fromDto0(end),   WBL_fromDto0(end),  'k.', 'MarkerSize', 10)
            % ---------------------------------------------------------------------
            text(V_fromStoA1(1),    WBL_fromStoA1(1),       '  S',  'FontSize', 6)
            text(V_fromStoA1(end),  WBL_fromStoA1(end),     '  A1', 'FontSize', 6)
            text(V_fromA1toC(end),  WBL_fromA1toC(end),    '  C1', 'FontSize', 6)
            text(V_fromCtoD(end),   WBL_fromCtoD(end),  '  C',  'FontSize', 6)
            text(V_fromDto0(end),   WBL_fromDto0(end),  '  C2', 'FontSize', 6)
            % text(40.75, -18, 'n = 1')
            % ---------------------------------------------------------------------
            xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
            ylabel("Wing lift - $L_{w}$ (daN)", "Interpreter", "latex")
            title("Wing lift per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct  
            
            % STORE INSIDE THE AIRCRAFT STRUCTURE VARIABLE

            % POINT S
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.value = q_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS.value = WBL_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS_new.value = WBL_from0toS_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTS.value = LHT_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTS.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.value = alfa_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_rad.value = deg2rad(alfa_from0toS(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_rad.Attributes.unit = "rad"; 
            alfaS = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_new.value = alfa_from0toS_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_new_rad.value = deg2rad(alfa_from0toS_fullvehicle(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_new_rad.Attributes.unit = "rad"; 
            alfaS_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_new.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value = CL_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S_new.value = CL_from0toS_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CLHT_S.value = CLHT_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CLHT_S.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CD_S.value = CD_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CD_S.Attributes.unit = "Non dimensional";  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CM_S.value = CMCG_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CM_S.Attributes.unit = "Non dimensional";     
            % POINT A1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.value = q_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1.value = WBL_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1_new.value = WBL_fromStoA1_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.value = LHT_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.value = alfa_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_rad.value = deg2rad(alfa_fromStoA1(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_rad.Attributes.unit = "rad";
            alfaA1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaA1_new.value = alfa_fromStoA1_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaA1_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaA1_new_rad.value = deg2rad(alfa_fromStoA1_fullvehicle(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaA1_new_rad.Attributes.unit = "rad"; 
            alfaA1_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_new.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.value = CL_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1_new.value = CL_fromStoA1_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CLHT_A1.value = CLHT_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CLHT_A1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CD_A1.value = CD_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CD_A1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CM_A1.value = CMCG_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CM_A1.Attributes.unit = "Non dimensional"; 
            % POINT C
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value = q_fromA1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.value = WBL_fromA1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.value = WBL_fromA1toC_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.value = LHT_fromA1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value = alfa_fromA1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_rad.value = deg2rad(alfa_fromA1toC(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_rad.Attributes.unit = "rad"; 
            alfaC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_new.value = alfa_fromA1toC_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_new_rad.value = deg2rad(alfa_fromA1toC_fullvehicle(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_new_rad.Attributes.unit = "rad"; 
            alfaC_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_new.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value = CL_fromA1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C_new.value = CL_fromA1toC_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CLHT_C.value = CLHT_fromA1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CLHT_C.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CD_C.value = CD_fromA1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CD_C.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CM_C.value = CMCG_fromA1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CM_C.Attributes.unit = "Non dimensional";        
            % POINT D
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value = q_fromCtoD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD.value = WBL_fromCtoD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.value = WBL_fromCtoD_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value = LHT_fromCtoD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value = alfa_fromCtoD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_rad.value = deg2rad(alfa_fromCtoD(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_rad.Attributes.unit = "rad";  
            alfaD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_new.value = alfa_fromCtoD_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_new_rad.value = deg2rad(alfa_fromCtoD_fullvehicle(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_new_rad.Attributes.unit = "rad";  
            alfaD_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_new.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value = CL_fromCtoD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D_new.value = CL_fromCtoD_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CLHT_D.value = CLHT_fromCtoD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CLHT_D.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_D.value = CD_fromCtoD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_D.Attributes.unit = "Non dimensional";      
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CM_D.value = CMCG_fromCtoD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CM_D.Attributes.unit = "Non dimensional";       
            
        end
     
    % CASE 2: VA lower than the intercept    
    case 'Case 2'
        % CL, ALFA, CD CALCULATIONS - POSITIVE LOAD FACTORS
        % ----------------------------------------------------------------- 
        % Calculation of the CL for the wing - body configuration. It must be
        % noticed that the calculations relative to lift coefficient in this
        % particular range of airspeed V and load factor n are going to give a
        % constant as a results, which is the maximum lift coefficient following
        % the lift curve of the aircrat. The previously defined CLMAX is a little
        % bit greater than the one obtained from V - n diagram flight condition to
        % take into account the extra lift produced by the horizontal empennage.
        % WING BODY LIFT CALCULATION
        % ----------------------------------------------------------------- 
        % In this section of the code the following formula is applied to evaluate
        % the wing-body lift with the formula 
        %        
        %   L = 0.5*rho*(V^2)*Sw*CL 
        %
        % using the previously calculated lift coefficients, coming from the 
        % various curves of the final envelope diagram.
        % =================================================================
        % FROM 0 TO S
        CL_from0toS      = zeros(length(V_from0toS), 1);
        alfa_from0toS    = zeros(length(V_from0toS), 1);
        CD_from0toS      = zeros(length(V_from0toS), 1);
        q_from0toS       = zeros(length(V_from0toS), 1);
        WBL_from0toS     = zeros(length(V_from0toS), 1);
        WBL_from0toS_new = zeros(length(V_from0toS), 1);
        CMCL_from0toS    = zeros(length(V_from0toS), 1);
        CMCD_from0toS    = zeros(length(V_from0toS), 1);
        CMCT_from0toS    = zeros(length(V_from0toS), 1);
        CMCG_from0toS    = zeros(length(V_from0toS), 1); 
        CLHT_from0toS    = zeros(length(V_from0toS), 1);
        LHT_from0toS     = zeros(length(V_from0toS), 1);
        % =================================================================
        CL_from0toS_fullvehicle   = zeros(length(V_from0toS), 1);
        alfa_from0toS_fullvehicle = zeros(length(V_from0toS), 1);
        % =================================================================
        for i = 1:length(V_from0toS)
            CL_from0toS(i)   = CLmax_func(rho0, V_from0toS(i), WS, n_from0toS(i));
            if CL_from0toS(i) < CL_max_fullmodel
                CL_from0toS(i) = CL_from0toS(i);
            elseif CL_from0toS(i) > CL_max_fullmodel
                CL_from0toS(i) = CL_max_fullmodel;
            end
            alfa_from0toS(i) = alpha_fullmodel(CL_from0toS(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
            CD_from0toS(i)   = polyval(p_CD_wb1, alfa_from0toS(i));
            q_from0toS(i)    = 0.5*rho0*(V_from0toS(i))^2;
            WBL_from0toS(i)  = q_from0toS(i)*S*CL_from0toS(i)*1e-1;
            CMCL_from0toS(i) = CLWB_contrib(obj1, CL_from0toS(i), deg2rad(alfa_from0toS(i)), XAC, XCG, bCG, MAC);
            CMCD_from0toS(i) = CDWB_contrib(obj1, CL_from0toS(i), deg2rad(alfa_from0toS(i)), XAC, XCG, bCG, MAC);
            CMCT_from0toS(i) = CT_contr(obj1, CD_from0toS(i), Thrust_axes, MAC);
            CMCG_from0toS(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_from0toS(i));
            % HORIZONTAL TAIL LIFT COEFFICIENT
            CLHT_from0toS(i) = CL_Tail(obj1, CMCL_from0toS(i), CMCD_from0toS(i), ...
                                             CMCT_from0toS(i), CMCG_from0toS(i), ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_from0toS(i)));
            % HORIZONTAL TAIL LIFT
            LHT_from0toS(i) = (0.5)*(V_from0toS(i)^2)*(S)*(rho0)*(CLHT_from0toS(i))*(1e-1);
            
            % DRAFT VERSION OF ITERATION 
            CL_tail         = CLHT_from0toS(i);
            CL_wb           = CL_from0toS(i);
            CL_new_from0toS = CL_wb - CL_tail;
            tol             = 1e-3; 
            n               = 1;
            while abs(CL_wb - CL_tail) > tol
               alfa_new_from0toS = alpha_fullmodel(CL_new_from0toS, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
               CD_wb             = polyval(p_CD_wb1, alfa_new_from0toS);
               CMCL_new          = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_from0toS), XAC, XCG, bCG, MAC);
               CMCD_new          = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_from0toS), XAC, XCG, bCG, MAC);
               CMCT_new          = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
               CMCG_new          = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
               CLHT_new          = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_new_from0toS));
               CL_tail           = CLHT_new;
               CL_new_from0toS   = CL_from0toS(i) - CL_tail;
               CL_wb             = CL_wb + (CL_new_from0toS - CL_wb) * 1e-1;
               n                 = n + 1;
               if n == 15
                   break
               end
            end
            CL_from0toS_fullvehicle(i)   = CL_new_from0toS; 
            CLHT_from0toS(i) = CL_tail;
%             alfa_from0toS(i) = alfa_new_from0toS;
            alfa_from0toS_fullvehicle(i) = alfa_new_from0toS;
            WBL_from0toS(i)     = q_from0toS(i) * S * CL_from0toS(i) * 1e-1;
            WBL_from0toS_new(i) = q_from0toS(i) * S * CL_from0toS_fullvehicle(i) * 1e-1;
            LHT_from0toS(i)     = (0.5)*(V_from0toS(i)^2)*(S)*(rho0)*(CLHT_from0toS(i))*(1e-1);
            CMCL_from0toS(i)    = CMCL_new;
            CMCD_from0toS(i)    = CMCD_new;
            CMCG_from0toS(i)    = CMCG_new;
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS.value = CL_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS_new.value = CL_from0toS_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS.value = alfa_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS_rad.value = deg2rad(alfa_from0toS);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS_rad.Attributes.unit = "rad";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS_new.value = alfa_from0toS_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS_new_rad.value = deg2rad(alfa_from0toS_fullvehicle);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS_new_rad.Attributes.unit = "rad";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toS.value = CD_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toS.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toS.value = q_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toS.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toS.value = WBL_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toS.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.value = WBL_from0toS_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_from0toS.value = CMCL_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_from0toS.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_from0toS.value = CMCD_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_from0toS.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_from0toS.value = CMCT_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_from0toS.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_from0toS.value = CMCG_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_from0toS.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_from0toS.value = CLHT_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_from0toS.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_from0toS.value = LHT_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_from0toS.Attributes.unit = "daN";        
        % =================================================================
        % FROM S TO A1
        CL_fromStoA1      = zeros(length(V_fromStoA1), 1);
        alfa_fromStoA1    = zeros(length(V_fromStoA1), 1);
        CD_fromStoA1      = zeros(length(V_fromStoA1), 1);
        q_fromStoA1       = zeros(length(V_fromStoA1), 1);
        WBL_fromStoA1     = zeros(length(V_fromStoA1), 1); 
        WBL_fromStoA1_new = zeros(length(V_fromStoA1), 1);
        CMCL_fromStoA1    = zeros(length(V_fromStoA1), 1);
        CMCD_fromStoA1    = zeros(length(V_fromStoA1), 1);
        CMCT_fromStoA1    = zeros(length(V_fromStoA1), 1);
        CMCG_fromStoA1    = zeros(length(V_fromStoA1), 1);    
        CLHT_fromStoA1    = zeros(length(V_fromStoA1), 1); 
        LHT_fromStoA1     = zeros(length(V_fromStoA1), 1);   
        % =================================================================
        CL_fromStoA1_fullvehicle   = zeros(length(V_fromStoA1), 1);
        alfa_fromStoA1_fullvehicle = zeros(length(V_fromStoA1), 1);
        % ================================================================= 
        for i = 1:length(V_fromStoA1)
            CL_fromStoA1(i)   = CLmax_func(rho0, V_fromStoA1(i), WS, n_fromStoA1(i));
            if CL_fromStoA1(i) < CL_max_fullmodel
                CL_fromStoA1(i) = CL_fromStoA1(i);
            elseif CL_fromStoA1(i) > CL_max_fullmodel
                CL_fromStoA1(i) = CL_max_fullmodel;
            end
            alfa_fromStoA1(i) = alpha_fullmodel(CL_fromStoA1(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
            CD_fromStoA1(i)   = polyval(p_CD_wb1, alfa_fromStoA1(i));
            q_fromStoA1(i)    = 0.5*rho0*(V_fromStoA1(i))^2;
            WBL_fromStoA1(i)  = q_fromStoA1(i)*S*CL_fromStoA1(i)*1e-1;   
            CMCL_fromStoA1(i) = CLWB_contrib(obj1, CL_fromStoA1(i), deg2rad(alfa_fromStoA1(i)), XAC, XCG, bCG, MAC);
            CMCD_fromStoA1(i) = CDWB_contrib(obj1, CL_fromStoA1(i), deg2rad(alfa_fromStoA1(i)), XAC, XCG, bCG, MAC);
            CMCT_fromStoA1(i) = CT_contr(obj1, CD_fromStoA1(i), Thrust_axes, MAC);
            CMCG_fromStoA1(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromStoA1(i));
            % HORIZONTAL TAIL LIFT COEFFICIENT
            CLHT_fromStoA1(i) = CL_Tail(obj1, CMCL_fromStoA1(i), CMCD_fromStoA1(i), ...
                                              CMCT_fromStoA1(i), CMCG_fromStoA1(i), ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_fromStoA1(i)));
            % HORIZONTAL TAIL LIFT
            LHT_fromStoA1(i) = (0.5)*(V_fromStoA1(i)^2)*(S)*(rho0)*(CLHT_fromStoA1(i))*(1e-1);
            
            % DRAFT VERSION OF ITERATION 
            CL_tail          = CLHT_fromStoA1(i);
            CL_wb            = CL_fromStoA1(i);
            CL_new_fromStoA1 = CL_wb - CL_tail;
            tol              = 1e-3; 
            n                = 1;
            while abs(CL_wb - CL_tail) > tol
               alfa_new_fromStoA1 = alpha_fullmodel(CL_new_fromStoA1, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
               CD_wb              = polyval(p_CD_wb1, alfa_new_fromStoA1);
               CMCL_new           = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromStoA1), XAC, XCG, bCG, MAC);
               CMCD_new           = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromStoA1), XAC, XCG, bCG, MAC);
               CMCT_new           = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
               CMCG_new           = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
               CLHT_new           = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromStoA1));
               CL_tail           = CLHT_new;
               CL_new_fromStoA1  = CL_fromStoA1(i) - CL_tail;
               CL_wb             = CL_wb + (CL_new_fromStoA1 - CL_wb) * 1e-1;
               n                 = n + 1;
               if n == 15
                   break
               end
            end
            CL_fromStoA1_fullvehicle(i)   = CL_new_fromStoA1; 
            CLHT_fromStoA1(i) = CL_tail;
%             alfa_fromStoA1(i) = alfa_new_fromStoA1;
            alfa_fromStoA1_fullvehicle(i) = alfa_new_fromStoA1;
            WBL_fromStoA1(i)  = q_fromStoA1(i) * S * CL_fromStoA1(i) * 1e-1;
            WBL_fromStoA1_new = q_fromStoA1(i) * S * CL_fromStoA1_fullvehicle(i) * 1e-1;
            LHT_fromStoA1(i)  = (0.5)*(V_fromStoA1(i)^2)*(S)*(rho0)*(CLHT_fromStoA1(i))*(1e-1);
            CMCL_fromStoA1(i) = CMCL_new;
            CMCD_fromStoA1(i) = CMCD_new;
            CMCG_fromStoA1(i) = CMCG_new;
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1.value = CL_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1_new.value = CL_fromStoA1_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1.value = alfa_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1_rad.value = deg2rad(alfa_fromStoA1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1_rad.Attributes.unit = "rad";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1_new.value = alfa_fromStoA1_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1_new_rad.value = deg2rad(alfa_fromStoA1_fullvehicle);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1_new_rad.Attributes.unit = "rad";          
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromStoA1.value = CD_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromStoA1.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromStoA1.value = q_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromStoA1.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromStoA1.value = WBL_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromStoA1.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromStoA1_new.value = WBL_fromStoA1_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromStoA1_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromStoA1.value = CMCL_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromStoA1.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromStoA1.value = CMCD_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromStoA1.Attributes.unit = "Non dimensional";      
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromStoA1.value = CMCT_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromStoA1.Attributes.unit = "Non dimensional";       
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromStoA1.value = CMCG_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromStoA1.Attributes.unit = "Non dimensional";         
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromStoA1.value = CLHT_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromStoA1.Attributes.unit = "Non dimensional";          
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromStoA1.value = LHT_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromStoA1.Attributes.unit = "daN";   
        % =================================================================  
        % POINT A 
        CL_A   = CLmax_func(rho0, VA, WS, nA);
            if CL_A < CL_max_fullmodel
                CL_A = CL_A;
            elseif CL_A > CL_max_fullmodel
                CL_A = CL_max_fullmodel;
            end
        alfa_A = alpha_fullmodel(CL_A, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
        CD_A   = polyval(p_CD_wb1, alfa_A);
        q_A    = 0.5*rho0*(VA)^2;
        WBL_A  = q_A*S*CL_A*1e-1; 
        CMCL_A = CLWB_contrib(obj1, CL_A, deg2rad(alfa_A), XAC, XCG, bCG, MAC);
        CMCD_A = CDWB_contrib(obj1, CL_A, deg2rad(alfa_A), XAC, XCG, bCG, MAC);
        CMCT_A = CT_contr(obj1, CD_A, Thrust_axes, MAC);
        CMCG_A = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_A);
        % HORIZONTAL TAIL LIFT COEFFICIENT
        CLHT_A = CL_Tail(obj1, CMCL_A, CMCD_A, CMCT_A, CMCG_A, l_ht, MAC, XAC, XCG, deg2rad(alfa_A)); 
        % HORIZONTAL TAIL LIFT
        LHT_A = (0.5)*(VA^2)*(S)*(rho0)*(CLHT_A)*(1e-1);
        
%         % NEW LIFT COEFFICIENT
%         CL_A_new = CL_A - CLHT_A;
%         % NEW LIFT COEFFICIENT
%         LW_A_new = q_A*S*CL_A_new*1e-1;
        
        % DRAFT VERSION OF ITERATION 
        CL_tail  = CLHT_A;
        CL_wb    = CL_A;
        CL_new_A = CL_wb - CL_tail;
        tol      = 1e-3; 
        n        = 1;
        while abs(CL_wb - CL_tail) > tol
           alfa_new_A = alpha_fullmodel(CL_new_A, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
           CD_wb      = polyval(p_CD_wb1, alfa_new_A);
           CMCL_new   = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_A), XAC, XCG, bCG, MAC);
           CMCD_new   = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_A), XAC, XCG, bCG, MAC);
           CMCT_new   = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
           CMCG_new   = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
           CLHT_new   = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                l_ht, MAC, XAC, XCG, deg2rad(alfa_new_A));
           CL_tail    = CLHT_new;
           CL_new_A   = CL_A - CL_tail;
           CL_wb      = CL_wb + (CL_new_A - CL_wb) * 1e-1;
           n          = n + 1;
           if n == 100
               break
           end
        end
        CL_A_fullvehicle   = CL_new_A; 
        CLHT_A = CL_tail;
        alfa_A_fullvehicle = alfa_new_A;
%         alfa_A = alfa_new_A;
        WBL_A     = q_A * S * CL_A * 1e-1;
        WBL_A_new = q_A * S * CL_A_fullvehicle * 1e-1;
        LHT_A     = (0.5) * ((VA)^2) * (S) * (rho0) * (CLHT_A) * (1e-1);
        CMCL_A    = CMCL_new;
        CMCD_A    = CMCD_new;
        CMCG_A    = CMCG_new;        
        
        % STORE ALL THE DATA INSIDE THE STRUCT VARIABLE
        % POINT A
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value = q_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.Attributes.unit = "Pa";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA.value = WBL_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA_new.value = WBL_A_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.value = LHT_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.value = alfa_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA_rad.value = deg2rad(alfa_A);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA_new.value = alfa_A_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA_new_rad.value = deg2rad(alfa_A_fullvehicle);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA_new_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A_new.value = CL_A_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value = CL_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CLHT_A.value = CLHT_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CLHT_A.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CD_A.value = CD_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CD_A.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CM_A.value = CMCG_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CM_A.Attributes.unit = "Non dimensional";          
        % =================================================================  
        % FROM A1 TO C
        CL_fromA1toC      = zeros(length(V_fromA1toC), 1);
        alfa_fromA1toC    = zeros(length(V_fromA1toC), 1);
        CD_fromA1toC      = zeros(length(V_fromA1toC), 1);
        q_fromA1toC       = zeros(length(V_fromA1toC), 1);
        WBL_fromA1toC     = zeros(length(V_fromA1toC), 1);   
        WBL_fromA1toC_new = zeros(length(V_fromA1toC), 1); 
        CMCL_fromA1toC    = zeros(length(V_fromA1toC), 1);
        CMCD_fromA1toC    = zeros(length(V_fromA1toC), 1);
        CMCT_fromA1toC    = zeros(length(V_fromA1toC), 1);
        CMCG_fromA1toC    = zeros(length(V_fromA1toC), 1);    
        CLHT_fromA1toC    = zeros(length(V_fromA1toC), 1); 
        LHT_fromA1toC     = zeros(length(V_fromA1toC), 1);
        % =================================================================
        CL_fromA1toC_fullvehicle   = zeros(length(V_fromA1toC), 1);
        alfa_fromA1toC_fullvehicle = zeros(length(V_fromA1toC), 1);
        % ================================================================= 
        for i = 1:length(V_fromA1toC)
            CL_fromA1toC(i)   = CLmax_func(rho0, V_fromA1toC(i), WS, n_fromA1toC(i));
            if CL_fromA1toC(i) < CL_max_fullmodel
                CL_fromA1toC(i) = CL_fromA1toC(i);
            elseif CL_fromA1toC(i) > CL_max_fullmodel
                CL_fromA1toC(i) = CL_max_fullmodel;
            end
            alfa_fromA1toC(i) = alpha_fullmodel(CL_fromA1toC(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
            CD_fromA1toC(i)   = polyval(p_CD_wb1, alfa_fromA1toC(i));
            q_fromA1toC(i)    = 0.5*rho0*(V_fromA1toC(i))^2;
            WBL_fromA1toC(i)  = q_fromA1toC(i)*S*CL_fromA1toC(i)*1e-1; 
            CMCL_fromA1toC(i) = CLWB_contrib(obj1, CL_fromA1toC(i), deg2rad(alfa_fromA1toC(i)), XAC, XCG, bCG, MAC);
            CMCD_fromA1toC(i) = CDWB_contrib(obj1, CL_fromA1toC(i), deg2rad(alfa_fromA1toC(i)), XAC, XCG, bCG, MAC);
            CMCT_fromA1toC(i) = CT_contr(obj1, CD_fromA1toC(i), Thrust_axes, MAC);
            CMCG_fromA1toC(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromA1toC(i));
            % HORIZONTAL TAIL LIFT COEFFICIENT
            CLHT_fromA1toC(i) = CL_Tail(obj1, CMCL_fromA1toC(i), CMCD_fromA1toC(i), ...
                                              CMCT_fromA1toC(i), CMCG_fromA1toC(i), ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_fromA1toC(i))); 
            % HORIZONTAL TAIL LIFT
            LHT_fromA1toC(i) = (0.5)*(V_fromA1toC(i)^2)*(S)*(rho0)*(CLHT_fromA1toC(i))*(1e-1);    
            
            % DRAFT VERSION OF ITERATION 
            CL_tail          = CLHT_fromA1toC(i);
            CL_wb            = CL_fromA1toC(i);
            CL_new_fromA1toC = CL_wb - CL_tail;
            tol              = 1e-3; 
            n                = 1;
            while abs(CL_wb - CL_tail) > tol
               alfa_new_fromA1toC = alpha_fullmodel(CL_new_fromA1toC, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
               CD_wb              = polyval(p_CD_wb1, alfa_new_fromA1toC);
               CMCL_new           = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromA1toC), XAC, XCG, bCG, MAC);
               CMCD_new           = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromA1toC), XAC, XCG, bCG, MAC);
               CMCT_new           = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
               CMCG_new           = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
               CLHT_new           = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromA1toC));
               CL_tail           = CLHT_new;
               CL_new_fromA1toC  = CL_fromA1toC(i) - CL_tail;
               CL_wb             = CL_wb + (CL_new_fromA1toC - CL_wb) * 1e-1;
               n                 = n + 1;
               if n == 15
                   break
               end
            end
            CL_fromA1toC_fullvehicle(i) = CL_new_fromA1toC; 
            CLHT_fromA1toC(i) = CL_tail;
            alfa_fromA1toC_fullvehicle(i) = alfa_new_fromA1toC;
%             alfa_fromA1toC(i) = alfa_new_fromA1toC;
            WBL_fromA1toC(i)  = q_fromA1toC(i) * S * CL_fromA1toC(i) * 1e-1;
            WBL_fromA1toC_new = q_fromA1toC(i) * S * CL_fromA1toC_fullvehicle(i) * 1e-1;
            LHT_fromA1toC(i)  = (0.5)*(V_fromA1toC(i)^2)*(S)*(rho0)*(CLHT_fromA1toC(i))*(1e-1);
            CMCL_fromA1toC(i) = CMCL_new;
            CMCD_fromA1toC(i) = CMCD_new;
            CMCG_fromA1toC(i) = CMCG_new;
            
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC.value = CL_fromA1toC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC_new.value = CL_fromA1toC_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC.value = alfa_fromA1toC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC.Attributes.unit = "deg"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC_rad.value = deg2rad(alfa_fromA1toC);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC_new.value = alfa_fromA1toC_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC_new.Attributes.unit = "deg"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC_new_rad.value = deg2rad(alfa_fromA1toC_fullvehicle);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC_new_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromA1toC.value = CD_fromA1toC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromA1toC.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromA1toC.value = q_fromA1toC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromA1toC.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromA1toC.value = WBL_fromA1toC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromA1toC.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC_new.value = WBL_fromA1toC_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC_new.Attributes.unit = "daN"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromA1toC.value = CMCL_fromA1toC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromA1toC.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromA1toC.value = CMCD_fromA1toC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromA1toC.Attributes.unit = "Non dimensional";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromA1toC.value = CMCT_fromA1toC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromA1toC.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromA1toC.value = CMCG_fromA1toC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromA1toC.Attributes.unit = "Non dimensional";    
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromA1toC.value = CLHT_fromA1toC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromA1toC.Attributes.unit = "Non dimensional";            
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromA1toC.value = LHT_fromA1toC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromA1toC.Attributes.unit = "daN";    
        % =================================================================  
        % FROM C TO A2
        CL_fromCtoA2      = zeros(length(V_fromCtoA2), 1);
        alfa_fromCtoA2    = zeros(length(V_fromCtoA2), 1);
        CD_fromCtoA2      = zeros(length(V_fromCtoA2), 1);
        q_fromCtoA2       = zeros(length(V_fromCtoA2), 1);
        WBL_fromCtoA2     = zeros(length(V_fromCtoA2), 1);  
        WBL_fromCtoA2_new = zeros(length(V_fromCtoA2), 1);
        CMCL_fromCtoA2    = zeros(length(V_fromCtoA2), 1);
        CMCD_fromCtoA2    = zeros(length(V_fromCtoA2), 1);
        CMCT_fromCtoA2    = zeros(length(V_fromCtoA2), 1);
        CMCG_fromCtoA2    = zeros(length(V_fromCtoA2), 1);    
        CLHT_fromCtoA2    = zeros(length(V_fromCtoA2), 1); 
        LHT_fromCtoA2     = zeros(length(V_fromCtoA2), 1);
        % =================================================================
        CL_fromCtoA2_fullvehicle   = zeros(length(V_fromCtoA2), 1);
        alfa_fromCtoA2_fullvehicle = zeros(length(V_fromA1toC), 1);
        % ================================================================= 
        for i = 1:length(V_fromCtoA2)
            CL_fromCtoA2(i)   = CLmax_func(rho0, V_fromCtoA2(i), WS, n_fromCtoA2(i));
            if CL_fromCtoA2(i) < CL_max_fullmodel
                CL_fromCtoA2(i) = CL_fromCtoA2(i);
            elseif CL_fromCtoA2(i) > CL_max_fullmodel
                CL_fromCtoA2(i) = CL_max_fullmodel;
            end
            alfa_fromCtoA2(i) = alpha_fullmodel(CL_fromCtoA2(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
            CD_fromCtoA2(i)   = polyval(p_CD_wb1, alfa_fromCtoA2(i));
            q_fromCtoA2(i)    = 0.5*rho0*(V_fromCtoA2(i))^2;
            WBL_fromCtoA2(i)  = q_fromCtoA2(i)*S*CL_fromCtoA2(i)*1e-1;    
            CMCL_fromCtoA2(i) = CLWB_contrib(obj1, CL_fromCtoA2(i), deg2rad(alfa_fromCtoA2(i)), XAC, XCG, bCG, MAC);    
            CMCD_fromCtoA2(i) = CDWB_contrib(obj1, CL_fromCtoA2(i), deg2rad(alfa_fromCtoA2(i)), XAC, XCG, bCG, MAC);
            CMCT_fromCtoA2(i) = CT_contr(obj1, CD_fromCtoA2(i), Thrust_axes, MAC);   
            CMCG_fromCtoA2(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromCtoA2(i));      
            % HORIZONTAL TAIL LIFT COEFFICIENT
            CLHT_fromCtoA2(i) = CL_Tail(obj1, CMCL_fromCtoA2(i), CMCD_fromCtoA2(i), ...
                                              CMCT_fromCtoA2(i), CMCG_fromCtoA2(i), ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_fromCtoA2(i)));   
            % HORIZONTAL TAIL LIFT
            LHT_fromCtoA2(i) = (0.5)*(V_fromCtoA2(i)^2)*(S)*(rho0)*(CLHT_fromCtoA2(i))*(1e-1); 
            
            % DRAFT VERSION OF ITERATION 
            CL_tail          = CLHT_fromCtoA2(i);
            CL_wb            = CL_fromCtoA2(i);
            CL_new_fromCtoA2 = CL_wb - CL_tail;
            tol              = 1e-3; 
            n                = 1;
            while abs(CL_wb - CL_tail) > tol
               alfa_new_fromCtoA2 = alpha_fullmodel(CL_new_fromCtoA2, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
               CD_wb              = polyval(p_CD_wb1, alfa_new_fromCtoA2);
               CMCL_new           = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromCtoA2), XAC, XCG, bCG, MAC);
               CMCD_new           = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromCtoA2), XAC, XCG, bCG, MAC);
               CMCT_new           = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
               CMCG_new           = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
               CLHT_new           = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromCtoA2));
               CL_tail           = CLHT_new;
               CL_new_fromCtoA2  = CL_fromCtoA2(i) - CL_tail;
               CL_wb             = CL_wb + (CL_new_fromCtoA2 - CL_wb) * 1e-1;
               n                 = n + 1;
               if n == 15
                   break
               end
            end
            CL_fromCtoA2_fullvehicle(i)   = CL_new_fromCtoA2; 
            CLHT_fromCtoA2(i) = CL_tail;
            alfa_fromCtoA2_fullvehicle(i) = alfa_new_fromCtoA2;
%             alfa_fromCtoA2(i) = alfa_new_fromCtoA2;
            WBL_fromCtoA2(i)  = q_fromCtoA2(i) * S * CL_fromCtoA2(i) * 1e-1;
            WBL_fromCtoA2_new = q_fromCtoA2(i) * S * CL_fromCtoA2_fullvehicle(i) * 1e-1;
            LHT_fromCtoA2(i)  = (0.5)*(V_fromCtoA2(i)^2)*(S)*(rho0)*(CLHT_fromCtoA2(i))*(1e-1);
            CMCL_fromCtoA2(i) = CMCL_new;
            CMCD_fromCtoA2(i) = CMCD_new;
            CMCG_fromCtoA2(i) = CMCG_new;            
            
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoA2.value = CL_fromCtoA2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoA2.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoA2_new.value = CL_fromCtoA2_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoA2_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoA2.value = alfa_fromCtoA2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoA2.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoA2_rad.value = deg2rad(alfa_fromCtoA2);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoA2_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoA2_new.value = alfa_fromCtoA2_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoA2_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoA2_new_rad.value = deg2rad(alfa_fromCtoA2_fullvehicle);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoA2_new_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtoA2.value = CD_fromCtoA2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtoA2.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromCtoA2.value = q_fromCtoA2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromCtoA2.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromCtoA2.value = WBL_fromCtoA2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromCtoA2.Attributes.unit = "daN";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoA2_new.value = WBL_fromCtoA2_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoA2_new.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromCtoA2.value = CMCL_fromCtoA2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromCtoA2.Attributes.unit = "Non dimensional";       
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromCtoA2.value = CMCD_fromCtoA2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromCtoA2.Attributes.unit = "Non dimensional";          
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromCtoA2.value = CMCT_fromCtoA2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromCtoA2.Attributes.unit = "Non dimensional";           
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromCtoA2.value = CMCG_fromCtoA2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromCtoA2.Attributes.unit = "Non dimensional";           
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromCtoA2.value = CLHT_fromCtoA2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromCtoA2.Attributes.unit = "Non dimensional";             
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromCtoA2.value = LHT_fromCtoA2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromCtoA2.Attributes.unit = "daN";  
        % =================================================================
        % FROM A2 TO D
        CL_fromA2toD      = zeros(length(V_fromA2toD), 1);
        alfa_fromA2toD    = zeros(length(V_fromA2toD), 1);
        CD_fromA2toD      = zeros(length(V_fromA2toD), 1);
        q_fromA2toD       = zeros(length(V_fromA2toD), 1);
        WBL_fromA2toD     = zeros(length(V_fromA2toD), 1);   
        WBL_fromA2toD_new = zeros(length(V_fromA2toD), 1);
        CMCL_fromA2toD    = zeros(length(V_fromA2toD), 1);
        CMCD_fromA2toD    = zeros(length(V_fromA2toD), 1);
        CMCT_fromA2toD    = zeros(length(V_fromA2toD), 1);
        CMCG_fromA2toD    = zeros(length(V_fromA2toD), 1);      
        CLHT_fromA2toD    = zeros(length(V_fromA2toD), 1);    
        LHT_fromA2toD     = zeros(length(V_fromA2toD), 1); 
        % =================================================================
        CL_fromA2toD_fullvehicle   = zeros(length(V_fromA2toD), 1);
        alfa_fromA2toD_fullvehicle = zeros(length(V_fromA1toC), 1);
        % =================================================================   
        for i = 1:length(V_fromA2toD)
            CL_fromA2toD(i)   = CLmax_func(rho0, V_fromA2toD(i), WS, n_fromA2toD(i));
            if CL_fromA2toD(i) < CL_max_fullmodel
                CL_fromA2toD(i) = CL_fromA2toD(i);
            elseif CL_fromA2toD(i) > CL_max_fullmodel
                CL_fromA2toD(i) = CL_max_fullmodel;
            end
            alfa_fromA2toD(i) = alpha_fullmodel(CL_fromA2toD(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
            CD_fromA2toD(i)   = polyval(p_CD_wb1, alfa_fromA2toD(i));
            q_fromA2toD(i)    = 0.5*rho0*(V_fromA2toD(i))^2;
            WBL_fromA2toD(i)  = q_fromA2toD(i)*S*CL_fromA2toD(i)*1e-1;   
            CMCL_fromA2toD(i) = CLWB_contrib(obj1, CL_fromA2toD(i), deg2rad(alfa_fromA2toD(i)), XAC, XCG, bCG, MAC);    
            CMCD_fromA2toD(i) = CDWB_contrib(obj1, CL_fromA2toD(i), deg2rad(alfa_fromA2toD(i)), XAC, XCG, bCG, MAC);
            CMCT_fromA2toD(i) = CT_contr(obj1, CD_fromA2toD(i), Thrust_axes, MAC);   
            CMCG_fromA2toD(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromA2toD(i));      
            % HORIZONTAL TAIL LIFT COEFFICIENT
            CLHT_fromA2toD(i) = CL_Tail(obj1, CMCL_fromA2toD(i), CMCD_fromA2toD(i), ...
                                              CMCT_fromA2toD(i), CMCG_fromA2toD(i), ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_fromA2toD(i)));   
            % HORIZONTAL TAIL LIFT
            LHT_fromA2toD(i) = (0.5)*(V_fromA2toD(i)^2)*(S)*(rho0)*(CLHT_fromA2toD(i))*(1e-1);
            
            % DRAFT VERSION OF ITERATION 
            CL_tail          = CLHT_fromA2toD(i);
            CL_wb            = CL_fromA2toD(i);
            CL_new_fromA2toD = CL_wb - CL_tail;
            tol              = 1e-3; 
            n                = 1;
            while abs(CL_wb - CL_tail) > tol
               alfa_new_fromA2toD = alpha_fullmodel(CL_new_fromA2toD, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
               CD_wb              = polyval(p_CD_wb1, alfa_new_fromA2toD);
               CMCL_new           = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromA2toD), XAC, XCG, bCG, MAC);
               CMCD_new           = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromA2toD), XAC, XCG, bCG, MAC);
               CMCT_new           = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
               CMCG_new           = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
               CLHT_new           = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromA2toD));
               CL_tail           = CLHT_new;
               CL_new_fromA2toD  = CL_fromA2toD(i) - CL_tail;
               CL_wb             = CL_wb + (CL_new_fromA2toD - CL_wb) * 1e-1;
               n                 = n + 1;
               if n == 15
                   break
               end
            end
            CL_fromA2toD_fullvehicle(i)   = CL_new_fromA2toD; 
            CLHT_fromA2toD(i) = CL_tail;
            alfa_fromA2toD_fullvehicle(i) = alfa_new_fromA2toD;
%             alfa_fromA2toD(i) = alfa_new_fromA2toD;
            WBL_fromA2toD(i)  = q_fromA2toD(i) * S * CL_fromA2toD(i) * 1e-1;
            WBL_fromA2toD_new = q_fromA2toD(i) * S * CL_fromA2toD_fullvehicle(i) * 1e-1;
            LHT_fromA2toD(i)  = (0.5)*(V_fromA2toD(i)^2)*(S)*(rho0)*(CLHT_fromA2toD(i))*(1e-1);
            CMCL_fromA2toD(i) = CMCL_new;
            CMCD_fromA2toD(i) = CMCD_new;
            CMCG_fromA2toD(i) = CMCG_new;             
            
        end        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA2toD.value = CL_fromA2toD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA2toD.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA2toD_new.value = CL_fromA2toD_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA2toD_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA2toD.value = alfa_fromA2toD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA2toD.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA2toD_rad.value = deg2rad(alfa_fromA2toD);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA2toD_rad.Attributes.unit = "rad";    
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA2toD_new.value = alfa_fromA2toD_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA2toD_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA2toD_new_rad.value = deg2rad(alfa_fromA2toD_fullvehicle);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA2toD_new_rad.Attributes.unit = "rad";      
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromA2toD.value = CD_fromA2toD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromA2toD.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromA2toD.value = q_fromA2toD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromA2toD.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromA2toD.value = WBL_fromA2toD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromA2toD.Attributes.unit = "daN"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA2toD_new.value = WBL_fromA2toD_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA2toD_new.Attributes.unit = "daN"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromA2toD.value = CMCL_fromA2toD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromA2toD.Attributes.unit = "Non dimensional";    
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromA2toD.value = CMCD_fromA2toD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromA2toD.Attributes.unit = "Non dimensional";    
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromA2toD.value = CMCT_fromA2toD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromA2toD.Attributes.unit = "Non dimensional";       
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromA2toD.value = CMCG_fromA2toD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromA2toD.Attributes.unit = "Non dimensional";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromA2toD.value = CLHT_fromA2toD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromA2toD.Attributes.unit = "Non dimensional";              
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromA2toD.value = LHT_fromA2toD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromA2toD.Attributes.unit = "daN";
        % =================================================================
        % FROM D TO 0
        CL_fromDto0      = zeros(length(V_fromDto0), 1);
        alfa_fromDto0    = zeros(length(V_fromDto0), 1);
        CD_fromDto0      = zeros(length(V_fromDto0), 1);
        q_fromDto0       = zeros(length(V_fromDto0), 1);
        WBL_fromDto0     = zeros(length(V_fromDto0), 1);  
        WBL_fromDto0_new = zeros(length(V_fromDto0), 1);   
        CMCL_fromDto0    = zeros(length(V_fromDto0), 1);  
        CMCD_fromDto0    = zeros(length(V_fromDto0), 1);
        CMCT_fromDto0    = zeros(length(V_fromDto0), 1); 
        CMCG_fromDto0    = zeros(length(V_fromDto0), 1);        
        CLHT_fromDto0    = zeros(length(V_fromDto0), 1);      
        LHT_fromDto0     = zeros(length(V_fromDto0), 1);  
        % =================================================================
        CL_fromDto0_fullvehicle    = zeros(length(V_fromDto0), 1);
        alfa_fromDto0_fullvehicle  = zeros(length(V_fromA1toC), 1);
        % =================================================================  
        for i = 1:length(V_fromDto0)
            CL_fromDto0(i)   = CLmax_func(rho0, V_fromDto0(i), WS, n_fromDto0(i));
            if CL_fromDto0(i) < CL_max_fullmodel
                CL_fromDto0(i) = CL_fromDto0(i);
            elseif CL_fromDto0(i) > CL_max_fullmodel
                CL_fromDto0(i) = CL_max_fullmodel;
            end
            alfa_fromDto0(i) = alpha_fullmodel(CL_fromDto0(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
            CD_fromDto0(i)   = polyval(p_CD_wb1, alfa_fromDto0(i));
            q_fromDto0(i)    = 0.5*rho0*(V_fromDto0(i))^2;
            WBL_fromDto0(i)  = q_fromDto0(i)*S*CL_fromDto0(i)*1e-1;     
            CMCL_fromDto0(i) = CLWB_contrib(obj1, CL_fromDto0(i), deg2rad(alfa_fromDto0(i)), XAC, XCG, bCG, MAC);      
            CMCD_fromDto0(i) = CDWB_contrib(obj1, CL_fromDto0(i), deg2rad(alfa_fromDto0(i)), XAC, XCG, bCG, MAC);
            CMCT_fromDto0(i) = CT_contr(obj1, CD_fromDto0(i), Thrust_axes, MAC);   
            CMCG_fromDto0(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromDto0(i));       
            % HORIZONTAL TAIL LIFT COEFFICIENT
            CLHT_fromDto0(i) = CL_Tail(obj1, CMCL_fromDto0(i), CMCD_fromDto0(i), ...
                                             CMCT_fromDto0(i), CMCG_fromDto0(i), ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_fromDto0(i)));    
            % HORIZONTAL TAIL LIFT
            LHT_fromDto0(i) = (0.5)*(V_fromDto0(i)^2)*(S)*(rho0)*(CLHT_fromDto0(i))*(1e-1);      
            
            % DRAFT VERSION OF ITERATION 
            CL_tail          = CLHT_fromDto0(i);
            CL_wb            = CL_fromDto0(i);
            CL_new_fromDto0 = CL_wb - CL_tail;
            tol              = 1e-3; 
            n                = 1;
            while abs(CL_wb - CL_tail) > tol
               alfa_new_fromDto0  = alpha_fullmodel(CL_new_fromDto0, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
               CD_wb              = polyval(p_CD_wb1, alfa_new_fromDto0);
               CMCL_new           = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromDto0), XAC, XCG, bCG, MAC);
               CMCD_new           = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromDto0), XAC, XCG, bCG, MAC);
               CMCT_new           = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
               CMCG_new           = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
               CLHT_new           = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromDto0));
               CL_tail           = CLHT_new;
               CL_new_fromDto0   = CL_fromDto0(i) - CL_tail;
               CL_wb             = CL_wb + (CL_new_fromDto0 - CL_wb) * 1e-1;
               n                 = n + 1;
               if n == 15
                   break
               end
            end
            CL_fromDto0_fullvehicle(i)   = CL_new_fromDto0; 
            CLHT_fromDto0(i) = CL_tail;
            alfa_fromDto0_fullvehicle(i) = alfa_new_fromDto0;
%             alfa_fromDto0(i) = alfa_new_fromDto0;
            WBL_fromDto0(i)  = q_fromDto0(i) * S * CL_fromDto0(i) * 1e-1;
            WBL_fromDto0_new = q_fromDto0(i) * S * CL_fromDto0_fullvehicle(i) * 1e-1;
            LHT_fromDto0(i)  = (0.5)*(V_fromDto0(i)^2)*(S)*(rho0)*(CLHT_fromDto0(i))*(1e-1);
            CMCL_fromDto0(i) = CMCL_new;
            CMCD_fromDto0(i) = CMCD_new;
            CMCG_fromDto0(i) = CMCG_new;              
            
        end           
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0.value = CL_fromDto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.value = CL_fromDto0_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0.value = alfa_fromDto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_rad.value = deg2rad(alfa_fromDto0);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_rad.Attributes.unit = "rad";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_new.value = alfa_fromDto0_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_new_rad.value = deg2rad(alfa_fromDto0_fullvehicle);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_new_rad.Attributes.unit = "rad";          
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDto0.value = CD_fromDto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDto0.Attributes.unit = "Non dimensional";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromDto0.value = q_fromDto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromDto0.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromDto0.value = WBL_fromDto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromDto0.Attributes.unit = "daN";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.value = WBL_fromDto0_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.Attributes.unit = "daN"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromDto0.value = CMCL_fromDto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromDto0.Attributes.unit = "Non dimensional";    
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromDto0.value = CMCD_fromDto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromDto0.Attributes.unit = "Non dimensional";     
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromDto0.value = CMCT_fromDto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromDto0.Attributes.unit = "Non dimensional";      
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromDto0.value = CMCG_fromDto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromDto0.Attributes.unit = "Non dimensional";      
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromDto0.value = CLHT_fromDto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromDto0.Attributes.unit = "Non dimensional";               
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromDto0.value = LHT_fromDto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromDto0.Attributes.unit = "daN";     
        % =================================================================
        % TAIL AIRLOADS AT UNIT LOAD - N = 1
        
        % AIRSPEED
        VS = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.value;
        VD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
        V_unit_load_factor = zeros(numb, 1);
        V_unit_load_factor(1) = VS;
        for i = 2:length(V_unit_load_factor)
            V_unit_load_factor(i) = V_unit_load_factor(i-1) + (VD - VS)*(1/length(V_unit_load_factor));
        end  
        
        % DYNAMIC PRESSURE
        q_unit_load_factor = zeros(length(V_unit_load_factor), 1);
        for i = 1:length(V_unit_load_factor)
            q_unit_load_factor(i) = 0.5*rho0*(V_unit_load_factor(i)^2);
            % q_unit_load_factor(i) = 0.5*rho0*(V_unit_load_factor(i)^2);
        end
        
        % WING BODY LIFT COEFFICIENT
        Aircraft_weight = ones(length(V_unit_load_factor), 1)*Mass*g;
        CLWB_unit_load_factor = zeros(length(V_unit_load_factor), 1);
        for i = 1:length(V_unit_load_factor)
            Aero_force = q_unit_load_factor(i)*S;
            CLWB_unit_load_factor(i) = ((Aircraft_weight(i))*(1/Aero_force));
        end
        
        % ALFA WING BODY UNIT LOAD
        alfanew_unit_load_factor = zeros(length(V_unit_load_factor), 1);
        for i = 1:length(V_unit_load_factor)
            alfanew_unit_load_factor(i) = alpha_fullmodel(CLWB_unit_load_factor(i), a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
        end
        
        % DRAG LIFT COEFFICIENT
        CD_unit_load_factor = zeros(length(V_unit_load_factor), 1);
        for i = 1:length(V_unit_load_factor)
            CD_unit_load_factor(i) = polyval(p_CD_wb1, alfa_new_fromDto0);
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value                   = V_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.Attributes.unit         = "m/s";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.value                   = q_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.Attributes.unit         = "Pa";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value                = CLWB_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.Attributes.unit      = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.value               = alfanew_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.Attributes.unit     = "deg"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value           = deg2rad(alfanew_unit_load_factor);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.Attributes.unit = "rad"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value                  = CD_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.Attributes.unit        = "Non dimensional"; 
        
%         % Interpolation coefficient from actual aerodynamic data
%         CL_supp    = str2num(Aircraft.Certification.Aerodynamic_data.CL.value);
%         alpha_supp = str2num(Aircraft.Certification.Aerodynamic_data.alpha.value);
%         x          = 0.027*ones(length(CL_supp), 1);
%         p = polyfit(CL_supp + x, alpha_supp, 2);
%         CLalfa_deg = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;
%         
%         % ALPHA CALCULATION 
%         alpha_unit_load_factor = zeros(length(V_unit_load_factor), 1);
%         for i = 1:length(alpha_unit_load_factor)
%             alpha_unit_load_factor(i) = alpha_calc(obj1, ...
%                                                    CLWB_unit_load_factor(i), CL0, CL_star, CLalfa_deg, p);
%         end
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.value = alpha_unit_load_factor;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.Attributes.unit = "deg";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value = deg2rad(alpha_unit_load_factor);
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.Attributes.unit = "rad"; 
        
        % WING BODY LIFT 
        WBLift_unit_load_factor = zeros(length(alfanew_unit_load_factor), 1);
        for i = 1:length(alfanew_unit_load_factor)   
            WBLift_unit_load_factor(i) = q_unit_load_factor(i) * S * (CLWB_unit_load_factor(i))*(1E-1);  
            % WBLift_unit_load_factor(i) = (0.5) * (V_unit_load_factor(i)^2) * rho0 * S * (CLWB_unit_load_factor(i))*(1E-1);    
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.value           = WBLift_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.Attributes.unit = "daN";    
        
        % -----------------------------------------------------------------
        alpha_unit_load_factor = alfanew_unit_load_factor;
        % EVALUATION OF TAIL LOADS CONTRIBUTION     
        CMCL_unit_load_factor             = zeros(length(alpha_unit_load_factor), 1);
        CMCD_unit_load_factor             = zeros(length(alpha_unit_load_factor), 1);
        CMCT_unit_load_factor             = zeros(length(alpha_unit_load_factor), 1);
        CMCG_unit_load_factor             = zeros(length(alpha_unit_load_factor), 1);
        CLHT_unit_load_factor             = zeros(length(alpha_unit_load_factor), 1);
        LHT_unit_load_factor              = zeros(length(alpha_unit_load_factor), 1);
        LW_unit_load_factor               = zeros(length(alpha_unit_load_factor), 1);
        CL_unit_load_factor_fullvehicle   = zeros(length(alpha_unit_load_factor), 1);
        alfa_unit_load_factor_fullvehicle = zeros(length(alpha_unit_load_factor), 1);
        LW_unit_load_factor_new           = zeros(length(alpha_unit_load_factor), 1);
        for i = 1:length(alpha_unit_load_factor)
            CMCL_unit_load_factor(i) = CLWB_contrib(obj1, CLWB_unit_load_factor(i), ...
                           deg2rad(alpha_unit_load_factor(i)), XAC, XCG, bCG, MAC);
            CMCD_unit_load_factor(i) = CDWB_contrib(obj1, CLWB_unit_load_factor(i), ...
                           deg2rad(alpha_unit_load_factor(i)), XAC, XCG, bCG, MAC);
            CMCT_unit_load_factor(i) = CT_contr(obj1, CD_unit_load_factor(i), Thrust_axes, MAC);      
            CMCG_unit_load_factor(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CMCG_unit_load_factor(i));   
            CLHT_unit_load_factor(i) = CL_Tail(obj1, CMCL_unit_load_factor(i), CMCD_unit_load_factor(i), ...
                                               CMCT_unit_load_factor(i), CMCG_unit_load_factor(i), ...
                                               l_ht, MAC, XAC, XCG, deg2rad(alpha_unit_load_factor(i)));
            LHT_unit_load_factor(i) = q_unit_load_factor(i) * S * (CLHT_unit_load_factor(i))*(1e-1);  
            
            % DRAFT VERSION OF ITERATION 
            CL_tail                 = CLHT_unit_load_factor(i);
            CL_wb                   = CLWB_unit_load_factor(i);
            CL_new_unit_load_factor = CL_wb - CL_tail;
            if CL_new_unit_load_factor > CL_max
                CL_new_unit_load_factor = CL_max;
            end
            tol                     = 1e-3; 
            n                       = 1;
            while abs(CL_wb - CL_tail) > tol
               alfa_new_unit_load_factor = alpha_fullmodel(CL_new_unit_load_factor, a2, b2, c2, CL_max, CL_star, CL0, CLalfa);
               CD_wb                     = polyval(p_CD_wb1, alfa_new_unit_load_factor);
               CMCL_new                  = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_unit_load_factor), XAC, XCG, bCG, MAC);
               CMCD_new                  = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_unit_load_factor), XAC, XCG, bCG, MAC);
               CMCT_new                  = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
               CMCG_new                  = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
               CLHT_new                  = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_new_unit_load_factor));
               CL_tail                   = CLHT_new;
               CL_new_unit_load_factor   = CLWB_unit_load_factor(i) - CL_tail;
               CL_wb                     = CL_wb + (CL_new_unit_load_factor - CL_wb) * 1e-1;
               n                         = n + 1;
               if n == 15
                   break
               end
            end
            CL_unit_load_factor_fullvehicle(i)   = CL_new_unit_load_factor; 
            CLHT_unit_load_factor(i)             = CL_tail;
            alfa_unit_load_factor_fullvehicle(i) = alfa_new_unit_load_factor;
%             alfa_fromA1toC(i) = alfa_new_fromA1toC;
            LW_unit_load_factor(i)     = q_unit_load_factor(i) * S * CLWB_unit_load_factor(i) * 1e-1;
            LW_unit_load_factor_new(i) = q_unit_load_factor(i) * S * CL_unit_load_factor_fullvehicle(i) * 1e-1;
            LHT_unit_load_factor(i)    = q_unit_load_factor(i) * S * (CLHT_unit_load_factor(i))*(1e-1);
            CMCL_unit_load_factor(i)   = CMCL_new;
            CMCD_unit_load_factor(i)   = CMCD_new;
            CMCG_unit_load_factor(i)   = CMCG_new;                         
            % LHT_unit_load_factor(i) = (0.5)*(V_unit_load_factor(i)^2) * S * rho0 * (CLHT_unit_load_factor(i))*(1e-1);
            % LHT_unit_load_factor(i) = (0.5)*(V_unit_load_factor(i)^2) * S_ht * rho0 * (CLHT_unit_load_factor(i))*(1e-1);
            % LW_unit_load_factor(i)  = WBLift_unit_load_factor(i) - LHT_unit_load_factor(i);                                                                                     
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_unit_load_factor.value = CMCL_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_unit_load_factor.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_unit_load_factor.value = CMCD_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_unit_load_factor.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_unit_load_factor.value = CMCT_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_unit_load_factor.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_unit_load_factor.value = CMCG_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_unit_load_factor.Attributes.unit = "Non dimensional";     
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_unit_load_factor.value = CLHT_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_unit_load_factor.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_unit_load_factor.value = LHT_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_unit_load_factor.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_unit_load_factor.value = LW_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_unit_load_factor.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_unit_load_factor_new.value = LW_unit_load_factor_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_unit_load_factor_new.Attributes.unit = "daN";
        
        % TAIL LOADS DIAGRAM - POSITIVE SIDE   
        HT_balancing_loads = figure(5); 
        hold on; grid on; grid minor; 
        
%         ylim([LHT_from0toS(1) LHT_fromDto0(1)])
%         xlim([0.0 V_fromDto0(1)])
        
        % plot(V_from0toS,  LHT_from0toS,  '-r', 'LineWidth', 1)
        plot(V_fromStoA1, LHT_fromStoA1, '-r', 'LineWidth', 1)
        plot(V_fromA1toC, LHT_fromA1toC, '-r', 'LineWidth', 1)
        plot(V_fromCtoA2, LHT_fromCtoA2, '-r', 'LineWidth', 1)
        plot(V_fromA2toD, LHT_fromA2toD, '-r', 'LineWidth', 1)
        plot(V_fromDto0,  LHT_fromDto0,  '-r', 'LineWidth', 1)
        % UNIT LOAD DISTR.
        plot(V_unit_load_factor, LHT_unit_load_factor,  '--k', 'LineWidth', 2)
        % ---------------------------------------------------------------------
        plot(V_from0toS(1),    LHT_from0toS(1),    'k.', 'MarkerSize', 10)
        plot(V_from0toS(end),  LHT_from0toS(end),  'k.', 'MarkerSize', 10)
        plot(V_fromStoA1(end), LHT_fromStoA1(end), 'k.', 'MarkerSize', 10)
        plot(VA,               LHT_A,              'k.', 'MarkerSize', 10)
        plot(V_fromA1toC(end), LHT_fromA1toC(end), 'k.', 'MarkerSize', 10)
        plot(V_fromCtoA2(end), LHT_fromCtoA2(end), 'k.', 'MarkerSize', 10)
        plot(V_fromA2toD(end), LHT_fromA2toD(end), 'k.', 'MarkerSize', 10)
        plot(V_fromDto0(end),  LHT_fromDto0(end),  'k.', 'MarkerSize', 10)
        % ---------------------------------------------------------------------
        text(V_fromStoA1(1),   LHT_fromStoA1(1),   '  S',  'FontSize', 6)
        text(V_fromStoA1(end), LHT_fromStoA1(end), '  A1', 'FontSize', 6)
        text(VA,               LHT_A,              '  A',  'FontSize', 6)
        text(V_fromA1toC(end), LHT_fromA1toC(end), '  C',  'FontSize', 6)
        text(V_fromCtoA2(end), LHT_fromCtoA2(end), '  A2', 'FontSize', 6)
        text(V_fromA2toD(end), LHT_fromA2toD(end), '  D',  'FontSize', 6)
        text(40.75, -18, 'n = 1')
        % ---------------------------------------------------------------------
% 
%             % --------------------------------------------------------------------
%             plot(V_from0toS,  WBL_from0toS,  '-r', 'LineWidth', 1)
%             plot(V_fromStoA1, WBL_fromStoA1, '-r', 'LineWidth', 1)
%             plot(V_fromA1toC, WBL_fromA1toC, '-r', 'LineWidth', 1)
%             plot(V_fromCtoA2, WBL_fromCtoA2, '-r', 'LineWidth', 1)
%             plot(V_fromA2toD, WBL_fromA2toD, '-r', 'LineWidth', 1)
%             plot(V_fromDto0,  WBL_fromDto0,  '-r', 'LineWidth', 1)
%             % UNIT LOAD DISTR.
%             plot(V_unit_load_factor, LW_unit_load_factor,  '--k', 'LineWidth', 2)
%             % ---------------------------------------------------------------------
%             plot(V_from0toS(1),     WBL_from0toS(1),     'k.', 'MarkerSize', 10)
%             plot(V_from0toS(end),   WBL_from0toS(end),   'k.', 'MarkerSize', 10)
%             plot(V_fromStoA1(end),  WBL_fromStoA1(end),  'k.', 'MarkerSize', 10)
%             plot(V_fromA1toC(end),  WBL_fromA1toC(end),  'k.', 'MarkerSize', 10)
%             plot(V_fromCtoA2(end),  WBL_fromCtoA2(end),   'k.', 'MarkerSize', 10)
%             plot(V_fromA2toD(end),  WBL_fromA2toD(end),   'k.', 'MarkerSize', 10)
%             plot(V_fromDto0(end),   WBL_fromDto0(end),   'k.', 'MarkerSize', 10)
%             % ---------------------------------------------------------------------
%             text(V_fromStoA1(1),    WBL_fromStoA1(1),    '  S',  'FontSize', 6)
%             text(V_fromStoA1(end),  WBL_fromStoA1(end),  '  A1', 'FontSize', 6)
%             text(V_fromA1toC(end),  WBL_fromA1toC(end),  '  C', 'FontSize', 6)
%             text(V_fromCtoA2(end),  WBL_fromCtoA2(end),  '  A2',  'FontSize', 6)
%             text(V_fromA2toD(end),  WBL_fromA2toD(end),  '  D', 'FontSize', 6)
%             % text(40.75, -18, 'n = 1')
%             % ---------------------------------------------------------------------        
        xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
        ylabel("Horizontal tail lift - $L_{ht}$ (daN)", "Interpreter", "latex")
        title("Horizontal empennage airloads per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
        
%         % MAIN WING LOADS DIAGRAM        
%         CL_from0toS_new  = CL_from0toS  + CLHT_from0toS;        
%         CL_fromStoA1_new = CL_fromStoA1 + CLHT_fromStoA1;        
%         CL_fromA1toC_new = CL_fromA1toC + CLHT_fromA1toC;        
%         CL_fromCtoA2_new = CL_fromCtoA2 + CLHT_fromCtoA2;        
%         CL_fromA2toD_new = CL_fromA2toD + CLHT_fromA2toD;        
%         CL_fromDto0_new  = CL_fromDto0  + CLHT_fromDto0;
        
%         % MAIN WING LOADS DIAGRAM        
%         CL_from0toS_new  = CL_from0toS  - CLHT_from0toS;        
%         CL_fromStoA1_new = CL_fromStoA1 - CLHT_fromStoA1;        
%         CL_fromA1toC_new = CL_fromA1toC - CLHT_fromA1toC;        
%         CL_fromCtoA2_new = CL_fromCtoA2 - CLHT_fromCtoA2;        
%         CL_fromA2toD_new = CL_fromA2toD - CLHT_fromA2toD;        
%         CL_fromDto0_new  = CL_fromDto0  - CLHT_fromDto0;
%         LW_from0toS_new  = zeros(length(CL_from0toS_new), 1);
%         LW_fromStoA1_new = zeros(length(CL_from0toS_new), 1);
%         LW_fromA1toC_new = zeros(length(CL_from0toS_new), 1);
%         LW_fromCtoA2_new = zeros(length(CL_from0toS_new), 1);
%         LW_fromA2toD_new = zeros(length(CL_from0toS_new), 1);
%         LW_fromDto0_new  = zeros(length(CL_from0toS_new), 1);
%         for i = 1:length(CL_from0toS_new)
%             LW_from0toS_new(i)  = (0.5)*(V_from0toS(i)^2)* rho0 * S *(CL_from0toS_new(i))*(1e-1);  
%             LW_fromStoA1_new(i) = (0.5)*(V_fromStoA1(i)^2)* rho0 * S *(CL_fromStoA1_new(i))*(1e-1);
%             LW_fromA1toC_new(i) = (0.5)*(V_fromA1toC(i)^2)* rho0 * S *(CL_fromA1toC_new(i))*(1e-1);
%             LW_fromCtoA2_new(i) = (0.5)*(V_fromCtoA2(i)^2)* rho0 * S *(CL_fromCtoA2_new(i))*(1e-1);
%             LW_fromA2toD_new(i) = (0.5)*(V_fromA2toD(i)^2)* rho0 * S *(CL_fromA2toD_new(i))*(1e-1);
%             LW_fromDto0_new(i)  = (0.5)*(V_fromDto0(i)^2)* rho0 * S *(CL_fromDto0_new(i))*(1e-1);
% %             LW_from0toS_new(i)  = (0.5)*(V_from0toS(i)^2)* rho0 * S *(CL_from0toS(i))*(1e-1);  
% %             LW_fromStoA1_new(i) = (0.5)*(V_fromStoA1(i)^2)* rho0 * S *(CL_fromStoA1(i))*(1e-1);
% %             LW_fromA1toC_new(i) = (0.5)*(V_fromA1toC(i)^2)* rho0 * S *(CL_fromA1toC(i))*(1e-1);
% %             LW_fromCtoA2_new(i) = (0.5)*(V_fromCtoA2(i)^2)* rho0 * S *(CL_fromCtoA2(i))*(1e-1);
% %             LW_fromA2toD_new(i) = (0.5)*(V_fromA2toD(i)^2)* rho0 * S *(CL_fromA2toD(i))*(1e-1);
% %             LW_fromDto0_new(i)  = (0.5)*(V_fromDto0(i)^2)* rho0 * S *(CL_fromDto0(i))*(1e-1);
%         end       
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS_new.value = CL_from0toS_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS_new.Attributes.unit = "Non dimensional";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1_new.value = CL_fromStoA1_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1_new.Attributes.unit = "Non dimensional";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC_new.value = CL_fromA1toC_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC_new.Attributes.unit = "Non dimensional";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoA2_new.value = CL_fromCtoA2_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoA2_new.Attributes.unit = "Non dimensional";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA2toD_new.value = CL_fromA2toD_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA2toD_new.Attributes.unit = "Non dimensional";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.value = CL_fromDto0_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.Attributes.unit = "Non dimensional";
%         % ------------------------------------------------------------------------------------------------------------------------
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.value = LW_from0toS_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.Attributes.unit = "daN";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromStoA1_new.value = LW_fromStoA1_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromStoA1_new.Attributes.unit = "daN";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC_new.value = LW_fromA1toC_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC_new.Attributes.unit = "daN";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoA2_new.value = LW_fromCtoA2_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoA2_new.Attributes.unit = "daN";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA2toD_new.value = LW_fromA2toD_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA2toD_new.Attributes.unit = "daN";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.value = LW_fromDto0_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.Attributes.unit = "daN";
        
        % TAIL LOADS DIAGRAM - POSITIVE SIDE   
        Wing_balancing_loads = figure(6); 
        hold on; grid on; grid minor; 
        
%         ylim([LHT_from0toS(1) LHT_fromDto0(1)])
%         xlim([0.0 V_fromDto0(1)])
        
%         plot(V_from0toS,  LW_from0toS_new,  '-r', 'LineWidth', 1)
%         plot(V_fromStoA1, LW_fromStoA1_new, '-r', 'LineWidth', 1)
%         plot(V_fromA1toC, LW_fromA1toC_new, '-r', 'LineWidth', 1)
%         plot(V_fromCtoA2, LW_fromCtoA2_new, '-r', 'LineWidth', 1)
%         plot(V_fromA2toD, LW_fromA2toD_new, '-r', 'LineWidth', 1)
%         plot(V_fromDto0,  LW_fromDto0_new,  '-r', 'LineWidth', 1)
%         % UNIT LOAD DISTR.
%         plot(V_unit_load_factor, LW_unit_load_factor,  '--k', 'LineWidth', 2)
%         % ---------------------------------------------------------------------
%         plot(V_from0toS(1),    LW_from0toS_new(1),    'k.', 'MarkerSize', 10)
%         plot(V_from0toS(end),  LW_from0toS_new(end),  'k.', 'MarkerSize', 10)
% %         plot(VA,               WBL_A,                 'k.', 'MarkerSize', 10)
%         plot(V_fromStoA1(end), LW_fromStoA1_new(end), 'k.', 'MarkerSize', 10)
%         plot(V_fromA1toC(end), LW_fromA1toC_new(end), 'k.', 'MarkerSize', 10)
%         plot(V_fromCtoA2(end), LW_fromCtoA2_new(end), 'k.', 'MarkerSize', 10)
%         plot(V_fromA2toD(end), LW_fromA2toD_new(end), 'k.', 'MarkerSize', 10)
%         plot(V_fromDto0(end),  LW_fromDto0_new(end),  'k.', 'MarkerSize', 10)
%         % ---------------------------------------------------------------------
%         text(V_fromStoA1(1),   LW_fromStoA1_new(1),   '  S',  'FontSize',  6)
% %         text(VA,               WBL_A,                  '  A',  'FontSize', 6)
%         text(V_fromStoA1(end), LW_fromStoA1_new(end), '  A1', 'FontSize',  6)
%         text(V_fromA1toC(end), LW_fromA1toC_new(end), '  C',  'FontSize',  6)
%         text(V_fromCtoA2(end), LW_fromCtoA2_new(end), '  A2', 'FontSize',  6)
%         text(V_fromA2toD(end), LW_fromA2toD_new(end), '  D',  'FontSize',  6)
%         % text(40.75, -18, 'n = 1')
%         % ---------------------------------------------------------------------

            % --------------------------------------------------------------------
            plot(V_from0toS,  WBL_from0toS,  '-r', 'LineWidth', 1)
            plot(V_fromStoA1, WBL_fromStoA1, '-r', 'LineWidth', 1)
            plot(V_fromA1toC, WBL_fromA1toC, '-r', 'LineWidth', 1)
            plot(V_fromCtoA2, WBL_fromCtoA2, '-r', 'LineWidth', 1)
            plot(V_fromA2toD, WBL_fromA2toD, '-r', 'LineWidth', 1)
            plot(V_fromDto0,  WBL_fromDto0,  '-r', 'LineWidth', 1)
            % UNIT LOAD DISTR.
            plot(V_unit_load_factor, LW_unit_load_factor,  '--k', 'LineWidth', 2)
            % ---------------------------------------------------------------------
            plot(V_from0toS(1)   ,  WBL_from0toS(1)   ,  'k.', 'MarkerSize', 10)
            plot(V_from0toS(end) ,  WBL_from0toS(end) ,  'k.', 'MarkerSize', 10)
            plot(V_fromStoA1(end),  WBL_fromStoA1(end),  'k.', 'MarkerSize', 10)
            plot(V_fromA1toC(end),  WBL_fromA1toC(end),  'k.', 'MarkerSize', 10)
            plot(V_fromCtoA2(end),  WBL_fromCtoA2(end),  'k.', 'MarkerSize', 10)
            plot(V_fromA2toD(end),  WBL_fromA2toD(end),  'k.', 'MarkerSize', 10)
            plot(V_fromDto0(end) ,  WBL_fromDto0(end) ,  'k.', 'MarkerSize', 10)
            plot(VA              ,  WBL_A             ,  'k.', 'MarkerSize', 10)
            % ---------------------------------------------------------------------
            text(V_fromStoA1(1)  ,  WBL_fromStoA1(1),    '  S',  'FontSize', 6)
            text(V_fromA1toC(end),  WBL_fromA1toC(end),  '  C', 'FontSize', 6)
            text(V_fromA2toD(end),  WBL_fromA2toD(end),  '  D', 'FontSize', 6)
            text(VA              ,  WBL_A             ,  '  A', 'FontSize', 6)
            % text(V_fromStoA1(end),  WBL_fromStoA1(end),  '  A1', 'FontSize', 6)
            % text(V_fromCtoA2(end),  WBL_fromCtoA2(end),  '  A2',  'FontSize', 6)
            % text(40.75, -18, 'n = 1')
            % ---------------------------------------------------------------------            
        
        xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
        ylabel("Wing lift - $L_{w}$ (daN)", "Interpreter", "latex")
        title("Wing airloads per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct  
        
        % STORE INSIDE THE AIRCRAFT STRUCTURE VARIABLE

        % POINT S
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.value = q_from0toS(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.Attributes.unit = "Pa";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS.value = WBL_from0toS(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS_new.value = WBL_from0toS_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTS.value = LHT_from0toS(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTS.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.value = alfa_from0toS(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_rad.value = deg2rad(alfa_from0toS(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_rad.Attributes.unit = "rad";
        alfaS = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_new.value = alfa_from0toS_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_new_rad.value = deg2rad(alfa_from0toS_fullvehicle(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_new_rad.Attributes.unit = "rad";
        alfaS_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_new.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value = CL_from0toS(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S_new.value = CL_from0toS_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CLHT_S.value = CLHT_from0toS(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CLHT_S.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CD_S.value = CD_from0toS(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CD_S.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CM_S.value = CMCG_from0toS(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CM_S.Attributes.unit = "Non dimensional";       
        % POINT A1
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.value = q_fromStoA1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.Attributes.unit = "Pa";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1.value = WBL_fromStoA1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1_new.value = WBL_fromStoA1_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.value = LHT_fromStoA1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.value = alfa_fromStoA1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_rad.value = deg2rad(alfa_fromStoA1(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_rad.Attributes.unit = "rad";
        alfaA1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_new.value = alfa_fromStoA1_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_new_rad.value = deg2rad(alfa_fromStoA1_fullvehicle(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_new_rad.Attributes.unit = "rad";
        alfaA1_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_new.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.value = CL_fromStoA1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1_new.value = CL_fromStoA1_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CLHT_A1.value = CLHT_fromStoA1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CLHT_A1.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CD_A1.value = CD_fromStoA1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CD_A1.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CM_A1.value = CMCG_fromStoA1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CM_A1.Attributes.unit = "Non dimensional";
        % POINT C
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value = q_fromA1toC(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.Attributes.unit = "Pa"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.value = WBL_fromA1toC(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.value = WBL_fromA1toC_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.value = LHT_fromA1toC(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value = alfa_fromA1toC(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_rad.value = deg2rad(alfa_fromA1toC(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_rad.Attributes.unit = "rad"; 
        alfaC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_new.value = alfa_fromA1toC_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_new_rad.value = deg2rad(alfa_fromA1toC_fullvehicle(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_new_rad.Attributes.unit = "rad";
        alfaC_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_new.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value = CL_fromA1toC(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C_new.value = CL_fromA1toC_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CLHT_C.value = CLHT_fromA1toC(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CLHT_C.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CD_C.value = CD_fromA1toC(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CD_C.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CM_C.value = CMCG_fromA1toC(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CM_C.Attributes.unit = "Non dimensional";        
        % POINT A2
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.qA2.value = q_fromCtoA2(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.qA2.Attributes.unit = "Pa";     
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.LA2.value = WBL_fromCtoA2(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.LA2.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.LA2_new.value = WBL_fromCtoA2_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.LA2_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.LHTA2.value = LHT_fromCtoA2(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.LHTA2.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.alfaA2.value = alfa_fromCtoA2(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.alfaA2.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.alfaA2_rad.value = deg2rad(alfa_fromCtoA2(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.alfaA2_rad.Attributes.unit = "rad";
        alfaA2 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.alfaA2.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.alfaA2_new.value = alfa_fromCtoA2_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.alfaA2_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.alfaA2_new_rad.value = deg2rad(alfa_fromCtoA2_fullvehicle(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.alfaA2_new_rad.Attributes.unit = "rad";
        alfaA2_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.alfaA2_new.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CL_A2.value = CL_fromCtoA2(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CL_A2.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CL_A2_new.value = CL_fromCtoA2_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CL_A2_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CLHT_A2.value = CLHT_fromCtoA2(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CLHT_A2.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CD_A2.value = CD_fromCtoA2(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CD_A2.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CM_A2.value = CMCG_fromCtoA2(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CM_A2.Attributes.unit = "Non dimensional";
        % POINT D
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value = q_fromA2toD(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.Attributes.unit = "Pa"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD.value = WBL_fromA2toD(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.value = WBL_fromA2toD_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value = LHT_fromA2toD(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value = alfa_fromA2toD(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_rad.value = deg2rad(alfa_fromA2toD(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_rad.Attributes.unit = "rad"; 
        alfaD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_new.value = alfa_fromA2toD_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_new_rad.value = deg2rad(alfa_fromA2toD_fullvehicle(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_new_rad.Attributes.unit = "rad";
        alfaD_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_new.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value = CL_fromA2toD(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D_new.value = CL_fromA2toD_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CLHT_D.value = CLHT_fromA2toD(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CLHT_D.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_D.value = CD_fromA2toD(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_D.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CM_D.value = CMCG_fromA2toD(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CM_D.Attributes.unit = "Non dimensional";      
        
end

Inverted_flight_Case = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Inverted_flight.value; 
switch (Inverted_flight_Case)
    % CASE 1: Complex solutions of the intercept
    case 'Case 1'
        if abs(min(n_gust_cruise_neg)) > abs(nmin)
            % S_inv - G - G1 - F - G2 - E
            % =============================================================
            % FROM 0 TO S_INVERTED 
            CL_from0toSinv      = zeros(length(V_from0toSinv), 1);
            alfa_from0toSinv    = zeros(length(V_from0toSinv), 1);
            CD_from0toSinv      = zeros(length(V_from0toSinv), 1);
            q_from0toSinv       = zeros(length(V_from0toSinv), 1);
            WBL_from0toSinv     = zeros(length(V_from0toSinv), 1);
            WBL_from0toSinv_new = zeros(length(V_from0toSinv), 1);
            CMCL_from0toSinv    = zeros(length(V_from0toSinv), 1);
            CMCD_from0toSinv    = zeros(length(V_from0toSinv), 1);
            CMCT_from0toSinv    = zeros(length(V_from0toSinv), 1);
            CMCG_from0toSinv    = zeros(length(V_from0toSinv), 1); 
            CLHT_from0toSinv    = zeros(length(V_from0toSinv), 1);
            LHT_from0toSinv     = zeros(length(V_from0toSinv), 1);
            % =================================================================
            CL_from0toSinv_fullvehicle   = zeros(length(V_from0toSinv), 1);
            alfa_from0toSinv_fullvehicle = zeros(length(V_from0toSinv), 1);
            % ================================================================= 
            for i = 1:length(V_from0toSinv)
                CL_from0toSinv(i)   = CLmax_func(rho0, V_from0toSinv(i), WS, n_from0toSinv(i));
                alfa_from0toSinv(i) = alpha_fullmodel_inv(CL_from0toSinv(i), CL_max_inv, CL0, CLalfa, alfa_0l);
                CD_from0toSinv(i)   = polyval(p_CD_wb1, alfa_from0toSinv(i));
                q_from0toSinv(i)    = 0.5*rho0*(V_from0toSinv(i))^2;
                WBL_from0toSinv(i)  = q_from0toSinv(i)*S*CL_from0toSinv(i)*1e-1;
                CMCL_from0toSinv(i) = CLWB_contrib(obj1, CL_from0toSinv(i), deg2rad(alfa_from0toSinv(i)), XAC, XCG, bCG, MAC);
                CMCD_from0toSinv(i) = CDWB_contrib(obj1, CL_from0toSinv(i), deg2rad(alfa_from0toSinv(i)), XAC, XCG, bCG, MAC);
                CMCT_from0toSinv(i) = CT_contr(obj1, CD_from0toSinv(i), Thrust_axes, MAC);
                CMCG_from0toSinv(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_from0toSinv(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_from0toSinv(i) = CL_Tail(obj1, CMCL_from0toSinv(i), CMCD_from0toSinv(i), ...
                                                 CMCT_from0toSinv(i), CMCG_from0toSinv(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_from0toSinv(i)));
                % HORIZONTAL TAIL LIFT
                LHT_from0toSinv(i) = (0.5)*(V_from0toSinv(i)^2)*(S)*(rho0)*(CLHT_from0toSinv(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail            = CLHT_from0toSinv(i);
                CL_wb              = CL_from0toSinv(i);
                CL_new_from0toSinv = CL_wb - CL_tail;
                tol                = 1e-3; 
                n                  = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_from0toSinv = alpha_fullmodel_inv(CL_new_from0toSinv, CL_max_inv, CL0, CLalfa, alfa_0l);
                   CD_wb                = polyval(p_CD_wb1, alfa_new_from0toSinv);
                   CMCL_new             = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_from0toSinv), XAC, XCG, bCG, MAC);
                   CMCD_new             = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_from0toSinv), XAC, XCG, bCG, MAC);
                   CMCT_new             = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new             = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new             = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_from0toSinv));
                   CL_tail              = CLHT_new;
                   CL_new_from0toSinv   = CL_from0toSinv(i) - CL_tail;
                   CL_wb                = CL_wb + (CL_new_from0toSinv - CL_wb) * 1e-1;
                   n                    = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_from0toSinv_fullvehicle(i)   = CL_new_from0toSinv; 
                CLHT_from0toSinv(i) = CL_tail;
                alfa_from0toSinv_fullvehicle(i) = alfa_new_from0toSinv;
%                 alfa_from0toSinv(i) = alfa_new_from0toSinv;
                WBL_from0toSinv(i)     = q_from0toSinv(i) * S * CL_from0toSinv(i) * 1e-1;
                WBL_from0toSinv_new(i) = q_from0toSinv(i) * S * CL_from0toSinv_fullvehicle(i) * 1e-1;
                LHT_from0toSinv(i)     = (0.5)*(V_from0toSinv(i)^2)*(S)*(rho0)*(CLHT_from0toSinv(i))*(1e-1);
                CMCL_from0toSinv(i)    = CMCL_new;
                CMCD_from0toSinv(i)    = CMCD_new;
                CMCG_from0toSinv(i)    = CMCG_new;
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv.value = CL_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.value = CL_from0toSinv_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv.value = alfa_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_rad.value = deg2rad(alfa_from0toSinv);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_rad.Attributes.unit = "rad";    
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_new.value = alfa_from0toSinv_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_new_rad.value = deg2rad(alfa_from0toSinv_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_new_rad.Attributes.unit = "rad";     
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toSinv.value = CD_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toSinv.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toSinv.value = q_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toSinv.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toSinv.value = WBL_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toSinv.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.value = WBL_from0toSinv_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_from0toSinv.value = CMCL_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_from0toSinv.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_from0toSinv.value = CMCD_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_from0toSinv.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_from0toSinv.value = CMCT_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_from0toSinv.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_from0toSinv.value = CMCG_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_from0toSinv.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_from0toSinv.value = CLHT_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_from0toSinv.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_from0toSinv.value = LHT_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_from0toSinv.Attributes.unit = "daN";
            % =============================================================
            % FROM S_INVERTED TO G
            CL_fromSinvtoG      = zeros(length(V_fromSinvtoG), 1);
            alfa_fromSinvtoG    = zeros(length(V_fromSinvtoG), 1);
            CD_fromSinvtoG      = zeros(length(V_fromSinvtoG), 1);
            q_fromSinvtoG       = zeros(length(V_fromSinvtoG), 1);
            WBL_fromSinvtoG     = zeros(length(V_fromSinvtoG), 1);
            WBL_fromSinvtoG_new = zeros(length(V_fromSinvtoG), 1);
            CMCL_fromSinvtoG    = zeros(length(V_fromSinvtoG), 1);
            CMCD_fromSinvtoG    = zeros(length(V_fromSinvtoG), 1);
            CMCT_fromSinvtoG    = zeros(length(V_fromSinvtoG), 1);
            CMCG_fromSinvtoG    = zeros(length(V_fromSinvtoG), 1); 
            CLHT_fromSinvtoG    = zeros(length(V_fromSinvtoG), 1);
            LHT_fromSinvtoG     = zeros(length(V_fromSinvtoG), 1);
            % =================================================================
            CL_fromSinvtoG_fullvehicle   = zeros(length(V_fromSinvtoG), 1);
            alfa_fromSinvtoG_fullvehicle = zeros(length(V_fromSinvtoG), 1);
            % ================================================================= 
            for i = 1:length(V_fromSinvtoG)
                CL_fromSinvtoG(i)   = CLmax_func(rho0, V_fromSinvtoG(i), WS, n_fromSinvtoG(i));
                alfa_fromSinvtoG(i) = alpha_fullmodel_inv(CL_fromSinvtoG(i), CL_max_inv, CL0, CLalfa, alfa_0l);
                CD_fromSinvtoG(i)   =  polyval(p_CD_wb1, alfa_fromSinvtoG(i));
                q_fromSinvtoG(i)    = 0.5*rho0*(V_fromSinvtoG(i))^2;
                WBL_fromSinvtoG(i)  = q_fromSinvtoG(i)*S*CL_fromSinvtoG(i)*1e-1;
                CMCL_fromSinvtoG(i) = CLWB_contrib(obj1, CL_fromSinvtoG(i), deg2rad(alfa_fromSinvtoG(i)), XAC, XCG, bCG, MAC);
                CMCD_fromSinvtoG(i) = CDWB_contrib(obj1, CL_fromSinvtoG(i), deg2rad(alfa_fromSinvtoG(i)), XAC, XCG, bCG, MAC);
                CMCT_fromSinvtoG(i) = CT_contr(obj1, CD_fromSinvtoG(i), Thrust_axes, MAC);
                CMCG_fromSinvtoG(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromSinvtoG(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromSinvtoG(i) = CL_Tail(obj1, CMCL_fromSinvtoG(i), CMCD_fromSinvtoG(i), ...
                                                 CMCT_fromSinvtoG(i), CMCG_fromSinvtoG(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromSinvtoG(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromSinvtoG(i) = (0.5)*(V_fromSinvtoG(i)^2)*(S)*(rho0)*(CLHT_fromSinvtoG(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail            = CLHT_fromSinvtoG(i);
                CL_wb              = CL_fromSinvtoG(i);
                CL_new_fromSinvtoG = CL_wb - CL_tail;
                tol                = 1e-3; 
                n                  = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromSinvtoG = alpha_fullmodel_inv(CL_new_fromSinvtoG, CL_max_inv, CL0, CLalfa, alfa_0l);
                   CD_wb                = polyval(p_CD_wb1, alfa_new_fromSinvtoG);
                   CMCL_new             = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromSinvtoG), XAC, XCG, bCG, MAC);
                   CMCD_new             = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromSinvtoG), XAC, XCG, bCG, MAC);
                   CMCT_new             = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new             = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new             = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromSinvtoG));
                   CL_tail              = CLHT_new;
                   CL_new_fromSinvtoG   = CL_fromSinvtoG(i) - CL_tail;
                   CL_wb                = CL_wb + (CL_new_fromSinvtoG - CL_wb) * 1e-1;
                   n                    = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromSinvtoG_fullvehicle(i)   = CL_new_fromSinvtoG; 
                CLHT_fromSinvtoG(i) = CL_tail;
                alfa_fromSinvtoG_fullvehicle(i) = alfa_new_fromSinvtoG;
%                 alfa_fromSinvtoG(i) = alfa_new_fromSinvtoG;
                WBL_fromSinvtoG(i)     = q_fromSinvtoG(i) * S * CL_fromSinvtoG(i) * 1e-1;
                WBL_fromSinvtoG_new(i) = q_fromSinvtoG(i) * S * CL_fromSinvtoG_fullvehicle(i) * 1e-1;
                LHT_fromSinvtoG(i)  = (0.5)*(V_fromSinvtoG(i)^2)*(S)*(rho0)*(CLHT_fromSinvtoG(i))*(1e-1);
                CMCL_fromSinvtoG(i) = CMCL_new;
                CMCD_fromSinvtoG(i) = CMCD_new;
                CMCG_fromSinvtoG(i) = CMCG_new;
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG.value = CL_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG_new.value = CL_fromSinvtoG_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG.value = alfa_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG_rad.value = deg2rad(alfa_fromSinvtoG);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG_rad.Attributes.unit = "rad";       
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG_new.value = alfa_fromSinvtoG_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG_new_rad.value = deg2rad(alfa_fromSinvtoG_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG_new_rad.Attributes.unit = "rad";    
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromSinvtoG.value = CD_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromSinvtoG.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromSinvtoG.value = q_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromSinvtoG.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromSinvtoG.value = WBL_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromSinvtoG.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromSinvtoG_new.value = WBL_fromSinvtoG_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromSinvtoG_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromSinvtoG.value = CMCL_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromSinvtoG.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromSinvtoG.value = CMCD_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromSinvtoG.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromSinvtoG.value = CMCT_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromSinvtoG.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromSinvtoG.value = CMCG_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromSinvtoG.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromSinvtoG.value = CLHT_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromSinvtoG.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromSinvtoG.value = LHT_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromSinvtoG.Attributes.unit = "daN";
            % =============================================================
            % FROM G TO G1
            CL_fromGtoG1      = zeros(length(V_fromGtoG1), 1);
            alfa_fromGtoG1    = zeros(length(V_fromGtoG1), 1);
            CD_fromGtoG1      = zeros(length(V_fromGtoG1), 1);
            q_fromGtoG1       = zeros(length(V_fromGtoG1), 1);
            WBL_fromGtoG1     = zeros(length(V_fromGtoG1), 1);
            WBL_fromGtoG1_new = zeros(length(V_fromGtoG1), 1);
            CMCL_fromGtoG1    = zeros(length(V_fromGtoG1), 1);
            CMCD_fromGtoG1    = zeros(length(V_fromGtoG1), 1);
            CMCT_fromGtoG1    = zeros(length(V_fromGtoG1), 1);
            CMCG_fromGtoG1    = zeros(length(V_fromGtoG1), 1); 
            CLHT_fromGtoG1    = zeros(length(V_fromGtoG1), 1);
            LHT_fromGtoG1     = zeros(length(V_fromGtoG1), 1);
            % =================================================================
            CL_fromGtoG1_fullvehicle   = zeros(length(V_fromGtoG1), 1);
            alfa_fromGtoG1_fullvehicle = zeros(length(V_fromGtoG1), 1);
            % ================================================================= 
            for i = 1:length(V_fromGtoG1) 
                CL_fromGtoG1(i)   = CLmax_func(rho0, V_fromGtoG1(i), WS, n_fromGtoG1(i));
                alfa_fromGtoG1(i) = alpha_fullmodel_inv(CL_fromGtoG1(i), CL_max_inv, CL0, CLalfa, alfa_0l);
                CD_fromGtoG1(i)   = polyval(p_CD_wb1, alfa_fromGtoG1(i));
                q_fromGtoG1(i)    = 0.5*rho0*(V_fromGtoG1(i))^2;
                WBL_fromGtoG1(i)  = q_fromGtoG1(i)*S*CL_fromGtoG1(i)*1e-1;
                CMCL_fromGtoG1(i) = CLWB_contrib(obj1, CL_fromGtoG1(i), deg2rad(alfa_fromGtoG1(i)), XAC, XCG, bCG, MAC);
                CMCD_fromGtoG1(i) = CDWB_contrib(obj1, CL_fromGtoG1(i), deg2rad(alfa_fromGtoG1(i)), XAC, XCG, bCG, MAC);
                CMCT_fromGtoG1(i) = CT_contr(obj1, CD_fromGtoG1(i), Thrust_axes, MAC);
                CMCG_fromGtoG1(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromGtoG1(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromGtoG1(i) = CL_Tail(obj1, CMCL_fromGtoG1(i), CMCD_fromGtoG1(i), ...
                                                 CMCT_fromGtoG1(i), CMCG_fromGtoG1(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromGtoG1(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromGtoG1(i) = (0.5)*(V_fromGtoG1(i)^2)*(S)*(rho0)*(CLHT_fromGtoG1(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail          = CLHT_fromGtoG1(i);
                CL_wb            = CL_fromGtoG1(i);
                CL_new_fromGtoG1 = CL_wb - CL_tail;
                tol              = 1e-3; 
                n                = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromGtoG1   = alpha_fullmodel_inv(CL_new_fromGtoG1, CL_max_inv, CL0, CLalfa, alfa_0l);
                   CD_wb                = polyval(p_CD_wb1, alfa_new_fromGtoG1);
                   CMCL_new             = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromGtoG1), XAC, XCG, bCG, MAC);
                   CMCD_new             = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromGtoG1), XAC, XCG, bCG, MAC);
                   CMCT_new             = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new             = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new             = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromGtoG1));
                   CL_tail              = CLHT_new;
                   CL_new_fromGtoG1     = CL_fromGtoG1(i) - CL_tail;
                   CL_wb                = CL_wb + (CL_new_fromGtoG1 - CL_wb) * 1e-1;
                   n                    = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromGtoG1_fullvehicle(i)   = CL_new_fromGtoG1; 
                CLHT_fromGtoG1(i) = CL_tail;
                alfa_fromGtoG1_fullvehicle(i) = alfa_new_fromGtoG1;
%                 alfa_fromGtoG1(i) = alfa_new_fromGtoG1;
                WBL_fromGtoG1(i)     = q_fromGtoG1(i) * S * CL_fromGtoG1(i) * 1e-1;
                WBL_fromGtoG1_new(i) = q_fromGtoG1(i) * S * CL_fromGtoG1_fullvehicle(i) * 1e-1;
                LHT_fromGtoG1(i)     = (0.5)*(V_fromGtoG1(i)^2)*(S)*(rho0)*(CLHT_fromGtoG1(i))*(1e-1);
                CMCL_fromGtoG1(i)    = CMCL_new;
                CMCD_fromGtoG1(i)    = CMCD_new;
                CMCG_fromGtoG1(i)    = CMCG_new;
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoG1.value = CL_fromGtoG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoG1.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoG1_new.value = CL_fromGtoG1_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoG1_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoG1.value = alfa_fromGtoG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoG1.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoG1_rad.value = deg2rad(alfa_fromGtoG1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoG1_rad.Attributes.unit = "rad";    
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoG1_new.value = alfa_fromGtoG1_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoG1_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoG1_new_rad.value = deg2rad(alfa_fromGtoG1_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoG1_new_rad.Attributes.unit = "rad";     
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromGtoG1.value = CD_fromGtoG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromGtoG1.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromGtoG1.value = q_fromGtoG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromGtoG1.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromGtoG1.value = WBL_fromGtoG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromGtoG1.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromGtoG1_new.value = WBL_fromGtoG1_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromGtoG1_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromGtoG1.value = CMCL_fromGtoG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromGtoG1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromGtoG1.value = CMCD_fromGtoG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromGtoG1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromGtoG1.value = CMCT_fromGtoG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromGtoG1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromGtoG1.value = CMCG_fromGtoG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromGtoG1.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromGtoG1.value = CLHT_fromGtoG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromGtoG1.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromGtoG1.value = LHT_fromGtoG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromGtoG1.Attributes.unit = "daN";
            % =============================================================
            % FROM G1 TO F 
            CL_fromG1toF   = zeros(length(V_fromG1toF), 1);
            alfa_fromG1toF = zeros(length(V_fromG1toF), 1);
            CD_fromG1toF   = zeros(length(V_fromG1toF), 1);
            q_fromG1toF    = zeros(length(V_fromG1toF), 1);
            WBL_fromG1toF  = zeros(length(V_fromG1toF), 1);
            WBL_fromG1toF_new = zeros(length(V_fromG1toF), 1);
            CMCL_fromG1toF = zeros(length(V_fromG1toF), 1);
            CMCD_fromG1toF = zeros(length(V_fromG1toF), 1);
            CMCT_fromG1toF = zeros(length(V_fromG1toF), 1);
            CMCG_fromG1toF = zeros(length(V_fromG1toF), 1); 
            CLHT_fromG1toF = zeros(length(V_fromG1toF), 1);
            LHT_fromG1toF  = zeros(length(V_fromG1toF), 1);
            % =================================================================
            CL_fromG1toF_fullvehicle   = zeros(length(V_fromG1toF), 1);
            alfa_fromG1toF_fullvehicle = zeros(length(V_fromG1toF), 1);
            % ================================================================= 
            for i = 1:length(V_fromG1toF)
                CL_fromG1toF(i)   = CLmax_func(rho0, V_fromG1toF(i), WS, n_fromG1toF(i));
                alfa_fromG1toF(i) = alpha_fullmodel_inv(CL_fromG1toF(i), CL_max_inv, CL0, CLalfa, alfa_0l);
                CD_fromG1toF(i)   = polyval(p_CD_wb1, alfa_fromG1toF(i));
                q_fromG1toF(i)    = 0.5*rho0*(V_fromG1toF(i))^2;
                WBL_fromG1toF(i)  = q_fromG1toF(i)*S*CL_fromG1toF(i)*1e-1;
                CMCL_fromG1toF(i) = CLWB_contrib(obj1, CL_fromG1toF(i), deg2rad(alfa_fromG1toF(i)), XAC, XCG, bCG, MAC);
                CMCD_fromG1toF(i) = CDWB_contrib(obj1, CL_fromG1toF(i), deg2rad(alfa_fromG1toF(i)), XAC, XCG, bCG, MAC);
                CMCT_fromG1toF(i) = CT_contr(obj1, CD_fromG1toF(i), Thrust_axes, MAC);
                CMCG_fromG1toF(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromG1toF(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromG1toF(i) = CL_Tail(obj1, CMCL_fromG1toF(i), CMCD_fromG1toF(i), ...
                                                 CMCT_fromG1toF(i), CMCG_fromG1toF(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromG1toF(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromG1toF(i) = (0.5)*(V_fromG1toF(i)^2)*(S)*(rho0)*(CLHT_fromG1toF(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail          = CLHT_fromG1toF(i);
                CL_wb            = CL_fromG1toF(i);
                CL_new_fromG1toF = CL_wb - CL_tail;
                tol              = 1e-3; 
                n                = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromG1toF   = alpha_fullmodel_inv(CL_new_fromG1toF, CL_max_inv, CL0, CLalfa, alfa_0l);
                   CD_wb                = polyval(p_CD_wb1, alfa_new_fromG1toF);
                   CMCL_new             = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromG1toF), XAC, XCG, bCG, MAC);
                   CMCD_new             = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromG1toF), XAC, XCG, bCG, MAC);
                   CMCT_new             = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new             = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new             = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromG1toF));
                   CL_tail              = CLHT_new;
                   CL_new_fromG1toF     = CL_fromG1toF(i) - CL_tail;
                   CL_wb                = CL_wb + (CL_new_fromG1toF - CL_wb) * 1e-1;
                   n                    = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromG1toF_fullvehicle(i)   = CL_new_fromG1toF; 
                CLHT_fromG1toF(i) = CL_tail;
                alfa_fromG1toF_fullvehicle(i) = alfa_new_fromG1toF;
%                 alfa_fromG1toF(i) = alfa_new_fromG1toF;
                WBL_fromG1toF(i)  = q_fromG1toF(i) * S * CL_fromG1toF(i) * 1e-1;
                WBL_fromG1toF_new(i)  = q_fromG1toF(i) * S * CL_fromG1toF_fullvehicle(i) * 1e-1;
                LHT_fromG1toF(i)  = (0.5)*(V_fromG1toF(i)^2)*(S)*(rho0)*(CLHT_fromG1toF(i))*(1e-1);
                CMCL_fromG1toF(i) = CMCL_new;
                CMCD_fromG1toF(i) = CMCD_new;
                CMCG_fromG1toF(i) = CMCG_new;
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF.value = CL_fromG1toF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF_new.value = CL_fromG1toF_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF.value = alfa_fromG1toF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF_rad.value = deg2rad(alfa_fromG1toF);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF_rad.Attributes.unit = "rad";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF_new.value = alfa_fromG1toF_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF_new_rad.value = deg2rad(alfa_fromG1toF_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF_new_rad.Attributes.unit = "rad";      
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromG1toF.value = CD_fromG1toF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromG1toF.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromG1toF.value = q_fromG1toF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromG1toF.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromG1toF.value = WBL_fromG1toF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromG1toF.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG1toF_new.value = WBL_fromG1toF_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG1toF_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromG1toF.value = CMCL_fromG1toF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromG1toF.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromG1toF.value = CMCD_fromG1toF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromG1toF.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromG1toF.value = CMCT_fromG1toF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromG1toF.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromG1toF.value = CMCG_fromG1toF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromG1toF.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromG1toF.value = CLHT_fromG1toF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromG1toF.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromG1toF.value = LHT_fromG1toF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromG1toF.Attributes.unit = "daN";  
            % =============================================================
            % FROM F TO G2 
            CL_fromFtoG2   = zeros(length(V_fromFtoG2), 1);
            alfa_fromFtoG2 = zeros(length(V_fromFtoG2), 1);
            CD_fromFtoG2   = zeros(length(V_fromFtoG2), 1);
            q_fromFtoG2    = zeros(length(V_fromFtoG2), 1);
            WBL_fromFtoG2  = zeros(length(V_fromFtoG2), 1);
            WBL_fromFtoG2_new = zeros(length(V_fromFtoG2), 1);
            CMCL_fromFtoG2 = zeros(length(V_fromFtoG2), 1);
            CMCD_fromFtoG2 = zeros(length(V_fromFtoG2), 1);
            CMCT_fromFtoG2 = zeros(length(V_fromFtoG2), 1);
            CMCG_fromFtoG2 = zeros(length(V_fromFtoG2), 1); 
            CLHT_fromFtoG2 = zeros(length(V_fromFtoG2), 1);
            LHT_fromFtoG2  = zeros(length(V_fromFtoG2), 1);
            % =================================================================
            CL_fromFtoG2_fullvehicle   = zeros(length(V_fromFtoG2), 1);
            alfa_fromFtoG2_fullvehicle = zeros(length(V_fromFtoG2), 1);
            % ================================================================= 
            for i = 1:length(V_fromFtoG2)
                CL_fromFtoG2(i)   = CLmax_func(rho0, V_fromFtoG2(i), WS, n_fromFtoG2(i));
                alfa_fromFtoG2(i) = alpha_fullmodel_inv(CL_fromFtoG2(i), CL_max_inv, CL0, CLalfa, alfa_0l);
                CD_fromFtoG2(i)   = polyval(p_CD_wb1, alfa_fromFtoG2(i));
                q_fromFtoG2(i)    = 0.5*rho0*(V_fromFtoG2(i))^2;
                WBL_fromFtoG2(i)  = q_fromFtoG2(i)*S*CL_fromFtoG2(i)*1e-1;
                CMCL_fromFtoG2(i) = CLWB_contrib(obj1, CL_fromFtoG2(i), deg2rad(alfa_fromFtoG2(i)), XAC, XCG, bCG, MAC);
                CMCD_fromFtoG2(i) = CDWB_contrib(obj1, CL_fromFtoG2(i), deg2rad(alfa_fromFtoG2(i)), XAC, XCG, bCG, MAC);
                CMCT_fromFtoG2(i) = CT_contr(obj1, CD_fromFtoG2(i), Thrust_axes, MAC);
                CMCG_fromFtoG2(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromFtoG2(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromFtoG2(i) = CL_Tail(obj1, CMCL_fromFtoG2(i), CMCD_fromFtoG2(i), ...
                                                 CMCT_fromFtoG2(i), CMCG_fromFtoG2(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromFtoG2(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromFtoG2(i) = (0.5)*(V_fromFtoG2(i)^2)*(S)*(rho0)*(CLHT_fromFtoG2(i))*(1e-1);
 
                % DRAFT VERSION OF ITERATION 
                CL_tail          = CLHT_fromFtoG2(i);
                CL_wb            = CL_fromFtoG2(i);
                CL_new_fromFtoG2 = CL_wb - CL_tail;
                tol              = 1e-3; 
                n                = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromFtoG2   = alpha_fullmodel_inv(CL_new_fromFtoG2, CL_max_inv, CL0, CLalfa, alfa_0l);
                   CD_wb                = polyval(p_CD_wb1, alfa_new_fromFtoG2);
                   CMCL_new             = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromFtoG2), XAC, XCG, bCG, MAC);
                   CMCD_new             = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromFtoG2), XAC, XCG, bCG, MAC);
                   CMCT_new             = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new             = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new             = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromFtoG2));
                   CL_tail              = CLHT_new;
                   CL_new_fromFtoG2     = CL_fromFtoG2(i) - CL_tail;
                   CL_wb                = CL_wb + (CL_new_fromFtoG2 - CL_wb) * 1e-1;
                   n                    = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromFtoG2_fullvehicle(i)   = CL_new_fromFtoG2; 
                CLHT_fromFtoG2(i) = CL_tail;
                alfa_fromFtoG2_fullvehicle(i) = alfa_new_fromFtoG2;
%                 alfa_fromFtoG2(i) = alfa_new_fromFtoG2;
                WBL_fromFtoG2(i)  = q_fromFtoG2(i) * S * CL_fromFtoG2(i) * 1e-1;
                WBL_fromFtoG2_new(i)  = q_fromFtoG2(i) * S * CL_fromFtoG2_fullvehicle(i) * 1e-1;
                LHT_fromFtoG2(i)  = (0.5)*(V_fromFtoG2(i)^2)*(S)*(rho0)*(CLHT_fromFtoG2(i))*(1e-1);
                CMCL_fromFtoG2(i) = CMCL_new;
                CMCD_fromFtoG2(i) = CMCD_new;
                CMCG_fromFtoG2(i) = CMCG_new;               
                
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoG2.value = CL_fromFtoG2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoG2.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoG2_new.value = CL_fromFtoG2_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoG2_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoG2.value = alfa_fromFtoG2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoG2.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoG2_rad.value = deg2rad(alfa_fromFtoG2);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoG2_rad.Attributes.unit = "rad";     
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoG2_new.value = alfa_fromFtoG2_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoG2_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoG2_new_rad.value = deg2rad(alfa_fromFtoG2_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoG2_new_rad.Attributes.unit = "rad";    
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoG2.value = CD_fromFtoG2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoG2.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromFtoG2.value = q_fromFtoG2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromFtoG2.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoG2.value = WBL_fromFtoG2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoG2.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoG2_new.value = WBL_fromFtoG2_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoG2_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromFtoG2.value = CMCL_fromFtoG2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromFtoG2.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromFtoG2.value = CMCD_fromFtoG2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromFtoG2.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromFtoG2.value = CMCT_fromFtoG2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromFtoG2.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromFtoG2.value = CMCG_fromFtoG2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromFtoG2.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromFtoG2.value = CLHT_fromFtoG2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromFtoG2.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromFtoG2.value = LHT_fromFtoG2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromFtoG2.Attributes.unit = "daN";
            % =============================================================
            % FROM G2 TO E 
            CL_fromG2toE   = zeros(length(V_fromG2toE), 1);
            alfa_fromG2toE = zeros(length(V_fromG2toE), 1);
            CD_fromG2toE   = zeros(length(V_fromG2toE), 1);
            q_fromG2toE    = zeros(length(V_fromG2toE), 1);
            WBL_fromG2toE  = zeros(length(V_fromG2toE), 1);
            WBL_fromG2toE_new = zeros(length(V_fromG2toE), 1);
            CMCL_fromG2toE = zeros(length(V_fromG2toE), 1);
            CMCD_fromG2toE = zeros(length(V_fromG2toE), 1);
            CMCT_fromG2toE = zeros(length(V_fromG2toE), 1);
            CMCG_fromG2toE = zeros(length(V_fromG2toE), 1); 
            CLHT_fromG2toE = zeros(length(V_fromG2toE), 1);
            LHT_fromG2toE  = zeros(length(V_fromG2toE), 1);
            % =================================================================
            CL_fromG2toE_fullvehicle   = zeros(length(V_fromG2toE), 1);
            alfa_fromG2toE_fullvehicle = zeros(length(V_fromG2toE), 1);
            % ================================================================= 
            for i = 1:length(V_fromG2toE)
                CL_fromG2toE(i)   = CLmax_func(rho0, V_fromG2toE(i), WS, n_fromG2toE(i));
                alfa_fromG2toE(i) = alpha_fullmodel_inv(CL_fromG2toE(i), CL_max_inv, CL0, CLalfa, alfa_0l);
                CD_fromG2toE(i)   = polyval(p_CD_wb1, alfa_fromG2toE(i));
                q_fromG2toE(i)    = 0.5*rho0*(V_fromG2toE(i))^2;
                WBL_fromG2toE(i)  = q_fromG2toE(i)*S*CL_fromG2toE(i)*1e-1;
                CMCL_fromG2toE(i) = CLWB_contrib(obj1, CL_fromG2toE(i), deg2rad(alfa_fromG2toE(i)), XAC, XCG, bCG, MAC);
                CMCD_fromG2toE(i) = CDWB_contrib(obj1, CL_fromG2toE(i), deg2rad(alfa_fromG2toE(i)), XAC, XCG, bCG, MAC);
                CMCT_fromG2toE(i) = CT_contr(obj1, CD_fromG2toE(i), Thrust_axes, MAC);
                CMCG_fromG2toE(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromG2toE(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromG2toE(i) = CL_Tail(obj1, CMCL_fromG2toE(i), CMCD_fromG2toE(i), ...
                                                 CMCT_fromG2toE(i), CMCG_fromG2toE(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromG2toE(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromG2toE(i) = (0.5)*(V_fromG2toE(i)^2)*(S)*(rho0)*(CLHT_fromG2toE(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail          = CLHT_fromG2toE(i);
                CL_wb            = CL_fromG2toE(i);
                CL_new_fromG2toE = CL_wb - CL_tail;
                tol              = 1e-3; 
                n                = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromG2toE   = alpha_fullmodel_inv(CL_new_fromG2toE, CL_max_inv, CL0, CLalfa, alfa_0l);
                   CD_wb                = polyval(p_CD_wb1, alfa_new_fromG2toE);
                   CMCL_new             = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromG2toE), XAC, XCG, bCG, MAC);
                   CMCD_new             = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromG2toE), XAC, XCG, bCG, MAC);
                   CMCT_new             = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new             = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new             = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromG2toE));
                   CL_tail              = CLHT_new;
                   CL_new_fromG2toE     = CL_fromG2toE(i) - CL_tail;
                   CL_wb                = CL_wb + (CL_new_fromG2toE - CL_wb) * 1e-1;
                   n                    = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromG2toE_fullvehicle(i)   = CL_new_fromG2toE; 
                CLHT_fromG2toE(i) = CL_tail;
                alfa_fromG2toE_fullvehicle(i) = alfa_new_fromG2toE;
%                 alfa_fromG2toE(i) = alfa_new_fromG2toE;
                WBL_fromG2toE(i)  = q_fromG2toE(i) * S * CL_fromG2toE(i) * 1e-1;
                WBL_fromG2toE_new(i)  = q_fromG2toE(i) * S * CL_fromG2toE_fullvehicle(i) * 1e-1;
                LHT_fromG2toE(i)  = (0.5)*(V_fromG2toE(i)^2)*(S)*(rho0)*(CLHT_fromG2toE(i))*(1e-1);
                CMCL_fromG2toE(i) = CMCL_new;
                CMCD_fromG2toE(i) = CMCD_new;
                CMCG_fromG2toE(i) = CMCG_new;  
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG2toE.value = CL_fromG2toE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG2toE.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG2toE_new.value = CL_fromG2toE_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG2toE_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG2toE.value = alfa_fromG2toE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG2toE.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG2toE_rad.value = deg2rad(alfa_fromG2toE);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG2toE_rad.Attributes.unit = "rad";      
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG2toE_new.value = alfa_fromG2toE_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG2toE_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG2toE_new_rad.value = deg2rad(alfa_fromG2toE_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG2toE_new_rad.Attributes.unit = "rad";    
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromG2toE.value = CD_fromG2toE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromG2toE.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromG2toE.value = q_fromG2toE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromG2toE.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromG2toE.value = WBL_fromG2toE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromG2toE.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG2toE_new.value = WBL_fromG2toE_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG2toE_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromG2toE.value = CMCL_fromG2toE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromG2toE.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromG2toE.value = CMCD_fromG2toE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromG2toE.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromG2toE.value = CMCT_fromG2toE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromG2toE.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromG2toE.value = CMCG_fromG2toE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromG2toE.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromG2toE.value = CLHT_fromG2toE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromG2toE.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromG2toE.value = LHT_fromG2toE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromG2toE.Attributes.unit = "daN";
            % =============================================================
            % FROM E TO 0
            CL_fromEto0   = zeros(length(V_fromEto0), 1);
            alfa_fromEto0 = zeros(length(V_fromEto0), 1);
            CD_fromEto0   = zeros(length(V_fromEto0), 1);
            q_fromEto0    = zeros(length(V_fromEto0), 1);
            WBL_fromEto0  = zeros(length(V_fromEto0), 1);
            WBL_fromEto0_new = zeros(length(V_fromEto0), 1);
            CMCL_fromEto0 = zeros(length(V_fromEto0), 1);
            CMCD_fromEto0 = zeros(length(V_fromEto0), 1);
            CMCT_fromEto0 = zeros(length(V_fromEto0), 1);
            CMCG_fromEto0 = zeros(length(V_fromEto0), 1); 
            CLHT_fromEto0 = zeros(length(V_fromEto0), 1);
            LHT_fromEto0  = zeros(length(V_fromEto0), 1);
            % =================================================================
            CL_fromEto0_fullvehicle   = zeros(length(V_fromEto0), 1);
            alfa_fromEto0_fullvehicle = zeros(length(V_fromEto0), 1);
            % ================================================================= 
            for i = 1:length(V_fromEto0) 
                CL_fromEto0(i)   = CLmax_func(rho0, V_fromEto0(i), WS, n_fromEto0(i));
                alfa_fromEto0(i) = alpha_fullmodel_inv(CL_fromEto0(i), CL_max_inv, CL0, CLalfa, alfa_0l);
                CD_fromEto0(i)   = polyval(p_CD_wb1, alfa_fromEto0(i));
                q_fromEto0(i)    = 0.5*rho0*(V_fromEto0(i))^2;
                WBL_fromEto0(i)  = q_fromEto0(i)*S*CL_fromEto0(i)*1e-1;
                CMCL_fromEto0(i) = CLWB_contrib(obj1, CL_fromEto0(i), deg2rad(alfa_fromEto0(i)), XAC, XCG, bCG, MAC);
                CMCD_fromEto0(i) = CDWB_contrib(obj1, CL_fromEto0(i), deg2rad(alfa_fromEto0(i)), XAC, XCG, bCG, MAC);
                CMCT_fromEto0(i) = CT_contr(obj1, CD_fromEto0(i), Thrust_axes, MAC);
                CMCG_fromEto0(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromEto0(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromEto0(i) = CL_Tail(obj1, CMCL_fromEto0(i), CMCD_fromEto0(i), ...
                                                 CMCT_fromEto0(i), CMCG_fromEto0(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromEto0(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromEto0(i) = (0.5)*(V_fromEto0(i)^2)*(S)*(rho0)*(CLHT_fromEto0(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail         = CLHT_fromEto0(i);
                CL_wb           = CL_fromEto0(i);
                CL_new_fromEto0 = CL_wb - CL_tail;
                tol             = 1e-3; 
                n               = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromEto0    = alpha_fullmodel_inv(CL_new_fromEto0, CL_max_inv, CL0, CLalfa, alfa_0l);
                   CD_wb                = polyval(p_CD_wb1, alfa_new_fromEto0);
                   CMCL_new             = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromEto0), XAC, XCG, bCG, MAC);
                   CMCD_new             = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromEto0), XAC, XCG, bCG, MAC);
                   CMCT_new             = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new             = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new             = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromEto0));
                   CL_tail              = CLHT_new;
                   CL_new_fromEto0     = CL_fromEto0(i) - CL_tail;
                   CL_wb                = CL_wb + (CL_new_fromEto0 - CL_wb) * 1e-1;
                   n                    = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromEto0_fullvehicle(i)   = CL_new_fromEto0; 
                CLHT_fromEto0(i) = CL_tail;
                alfa_fromEto0_fullvehicle(i) = alfa_new_fromEto0;
%                 alfa_fromEto0(i) = alfa_new_fromEto0;
                WBL_fromEto0(i)  = q_fromEto0(i) * S * CL_fromEto0(i) * 1e-1;
                WBL_fromEto0_new(i)  = q_fromEto0(i) * S * CL_fromEto0(i) * 1e-1;
                LHT_fromEto0(i)  = (0.5)*(V_fromEto0(i)^2)*(S)*(rho0)*(CLHT_fromEto0(i))*(1e-1);
                CMCL_fromEto0(i) = CMCL_new;
                CMCD_fromEto0(i) = CMCD_new;
                CMCG_fromEto0(i) = CMCG_new;  
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0.value = CL_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.value = CL_fromEto0_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0.value = alfa_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_rad.value = deg2rad(alfa_fromEto0);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_rad.Attributes.unit = "rad";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_new.value = alfa_fromEto0_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_new_rad.value = deg2rad(alfa_fromEto0_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_new_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromEto0.value = CD_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromEto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromEto0.value = q_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromEto0.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromEto0.value = WBL_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromEto0.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.value = WBL_fromEto0_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromEto0.value = CMCL_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromEto0.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromEto0.value = CMCD_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromEto0.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromEto0.value = CMCT_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromEto0.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromEto0.value = CMCG_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromEto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromEto0.value = CLHT_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromEto0.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromEto0.value = LHT_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromEto0.Attributes.unit = "daN"; 
            % =================================================================
            % TAIL LOADS DIAGRAM - NEGATIVE SIDE 
            figure(5); 
            % hold on; grid on; grid minor; 
            % ---------------------------------------------------------------------
            plot(V_from0toSinv,  LHT_from0toSinv,  '-r', 'LineWidth', 1)
            plot(V_fromSinvtoG, LHT_fromSinvtoG,   '-r', 'LineWidth', 1)
            plot(V_fromGtoG1,    LHT_fromGtoG1,    '-r', 'LineWidth', 1)
            plot(V_fromG1toF,    LHT_fromG1toF,    '-r', 'LineWidth', 1)
            plot(V_fromFtoG2,    LHT_fromFtoG2,    '-r', 'LineWidth', 1)
            plot(V_fromG2toE,    LHT_fromG2toE,    '-r', 'LineWidth', 1)
            plot(V_fromEto0,     LHT_fromEto0,     '-r', 'LineWidth', 1)
            % ---------------------------------------------------------------------
            plot(V_from0toSinv(1),   LHT_from0toSinv(1),   'k.', 'MarkerSize', 10)
            plot(V_from0toSinv(end), LHT_from0toSinv(end), 'k.', 'MarkerSize', 10)
            plot(V_fromSinvtoG(end), LHT_fromSinvtoG(end), 'k.', 'MarkerSize', 10)
            plot(V_fromGtoG1(end),   LHT_fromGtoG1(end),   'k.', 'MarkerSize', 10)
            plot(V_fromG1toF(end),   LHT_fromG1toF(end),   'k.', 'MarkerSize', 10)
            plot(V_fromFtoG2(end),   LHT_fromFtoG2(end),   'k.', 'MarkerSize', 10)
            plot(V_fromG2toE(end),   LHT_fromG2toE(end),   'k.', 'MarkerSize', 10)
            % ---------------------------------------------------------------------
            text(V_from0toSinv(end),  LHT_from0toSinv(end), ' S inv.', 'FontSize', 6)
            text(V_fromSinvtoG(end), LHT_fromSinvtoG(end),  ' G',      'FontSize', 6)
            text(V_fromGtoG1(end),    LHT_fromGtoG1(end),   ' G1',     'FontSize', 6)
            text(V_fromG1toF(end),     LHT_fromG1toF(end),  ' F',      'FontSize', 6)
            text(V_fromFtoG2(end),     LHT_fromFtoG2(end),  ' G2',     'FontSize', 6)
            text(V_fromG2toE(end),     LHT_fromG2toE(end),  ' E',      'FontSize', 6)
            % =================================================================
%             % MAIN WING LOADS DIAGRAM        
%             CL_from0toSinv_new = CL_from0toSinv  + CLHT_from0toSinv;        
%             CL_fromSinvtoG_new = CL_fromSinvtoG + CLHT_fromSinvtoG;        
%             CL_fromGtoG1_new   = CL_fromGtoG1    + CLHT_fromGtoG1;        
%             CL_fromG1toF_new   = CL_fromG1toF    + CLHT_fromG1toF;        
%             CL_fromFtoG2_new   = CL_fromFtoG2    + CLHT_fromFtoG2;        
%             CL_fromG2toE_new   = CL_fromG2toE    + CLHT_fromG2toE;        
%             CL_fromEto0_new    = CL_fromEto0     + CLHT_fromEto0;

%             % MAIN WING LOADS DIAGRAM        
%             CL_from0toSinv_new = CL_from0toSinv  - CLHT_from0toSinv;        
%             CL_fromSinvtoG_new = CL_fromSinvtoG - CLHT_fromSinvtoG;        
%             CL_fromGtoG1_new   = CL_fromGtoG1    - CLHT_fromGtoG1;        
%             CL_fromG1toF_new   = CL_fromG1toF    - CLHT_fromG1toF;        
%             CL_fromFtoG2_new   = CL_fromFtoG2    - CLHT_fromFtoG2;        
%             CL_fromG2toE_new   = CL_fromG2toE    - CLHT_fromG2toE;        
%             CL_fromEto0_new    = CL_fromEto0     - CLHT_fromEto0;
%             LW_from0toSinv_new = zeros(length(CL_from0toSinv_new), 1);
%             LW_fromSinvtoG_new = zeros(length(CL_from0toSinv_new), 1);
%             LW_fromGtoG1_new   = zeros(length(CL_from0toSinv_new), 1);
%             LW_fromG1toF_new   = zeros(length(CL_from0toSinv_new), 1);
%             LW_fromFtoG2_new   = zeros(length(CL_from0toSinv_new), 1);
%             LW_fromG2toE_new   = zeros(length(CL_from0toSinv_new), 1);
%             LW_fromEto0_new    = zeros(length(CL_from0toSinv_new), 1);
%             for i = 1:length(CL_from0toSinv_new)
%                 LW_from0toSinv_new(i) = (0.5)*(V_from0toSinv(i)^2)* rho0 * S *(CL_from0toSinv_new(i))*(1e-1);  
%                 LW_fromSinvtoG_new(i) = (0.5)*(V_fromSinvtoG(i)^2)* rho0 * S *(CL_fromSinvtoG_new(i))*(1e-1);
%                 LW_fromGtoG1_new(i)   = (0.5)*(V_fromGtoG1(i)^2)* rho0 * S *(CL_fromGtoG1_new(i))*(1e-1);
%                 LW_fromG1toF_new(i)   = (0.5)*(V_fromG1toF(i)^2)* rho0 * S *(CL_fromG1toF_new(i))*(1e-1);
%                 LW_fromFtoG2_new(i)   = (0.5)*(V_fromFtoG2(i)^2)* rho0 * S *(CL_fromFtoG2_new(i))*(1e-1);
%                 LW_fromG2toE_new(i)   = (0.5)*(V_fromG2toE(i)^2)* rho0 * S *(CL_fromG2toE_new(i))*(1e-1);
%                 LW_fromEto0_new(i)    = (0.5)*(V_fromEto0(i)^2)* rho0 * S *(CL_fromEto0_new(i))*(1e-1);
% %                 LW_from0toSinv_new(i) = (0.5)*(V_from0toSinv(i)^2)* rho0 * S *(CL_from0toSinv(i))*(1e-1);  
% %                 LW_fromSinvtoG_new(i) = (0.5)*(V_fromSinvtoG(i)^2)* rho0 * S *(CL_fromSinvtoG(i))*(1e-1);
% %                 LW_fromGtoG1_new(i)   = (0.5)*(V_fromGtoG1(i)^2)* rho0 * S *(CL_fromGtoG1(i))*(1e-1);
% %                 LW_fromG1toF_new(i)   = (0.5)*(V_fromG1toF(i)^2)* rho0 * S *(CL_fromG1toF(i))*(1e-1);
% %                 LW_fromFtoG2_new(i)   = (0.5)*(V_fromFtoG2(i)^2)* rho0 * S *(CL_fromFtoG2(i))*(1e-1);
% %                 LW_fromG2toE_new(i)   = (0.5)*(V_fromG2toE(i)^2)* rho0 * S *(CL_fromG2toE(i))*(1e-1);
% %                 LW_fromEto0_new(i)    = (0.5)*(V_fromEto0(i)^2)* rho0 * S *(CL_fromEto0(i))*(1e-1);
%             end       
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.value = CL_from0toSinv_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG_new.value = CL_fromSinvtoG_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoG1_new.value = CL_fromGtoG1_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoG1_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF_new.value = CL_fromG1toF_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoG2_new.value = CL_fromFtoG2_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoG2_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG2toE_new.value = CL_fromG2toE_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG2toE_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.value = CL_fromEto0_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.Attributes.unit = "Non dimensional";
%             % ------------------------------------------------------------------------------------------------------------------------
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.value = LW_from0toSinv_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromSinvtoG_new.value = LW_fromSinvtoG_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromSinvtoG_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromGtoG1_new.value = LW_fromGtoG1_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromGtoG1_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG1toF_new.value = LW_fromG1toF_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG1toF_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoG2_new.value = LW_fromFtoG2_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoG2_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG2toE_new.value = LW_fromG2toE_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG2toE_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.value = LW_fromEto0_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.Attributes.unit = "daN"; 
            % =============================================================
            % WING LOADS DIAGRAM - NEGATIVE SIDE   
            figure(6)
%             % ---------------------------------------------------------------------
%             plot(V_from0toSinv, LW_from0toSinv_new, '-r', 'LineWidth', 1)
%             plot(V_fromSinvtoG, LW_fromSinvtoG_new, '-r', 'LineWidth', 1)
%             plot(V_fromGtoG1,   LW_fromGtoG1_new,   '-r', 'LineWidth', 1)
%             plot(V_fromG1toF,   LW_fromG1toF_new,   '-r', 'LineWidth', 1)
%             plot(V_fromFtoG2,   LW_fromFtoG2_new,   '-r', 'LineWidth', 1)
%             plot(V_fromG2toE,   LW_fromG2toE_new,   '-r', 'LineWidth', 1)
%             plot(V_fromEto0,    LW_fromEto0_new,    '-r', 'LineWidth', 1)
%             % ---------------------------------------------------------------------
%             plot(V_from0toSinv(1),   LW_from0toSinv_new(1),   'k.', 'MarkerSize', 10)
%             plot(V_from0toSinv(end), LW_from0toSinv_new(end), 'k.', 'MarkerSize', 10)
%             plot(V_fromSinvtoG(end), LW_fromSinvtoG_new(end), 'k.', 'MarkerSize', 10)
%             plot(V_fromGtoG1(end),   LW_fromGtoG1_new(end),   'k.', 'MarkerSize', 10)
%             plot(V_fromG1toF(end),   LW_fromG1toF_new(end),   'k.', 'MarkerSize', 10)
%             plot(V_fromFtoG2(end),   LW_fromFtoG2_new(end),   'k.', 'MarkerSize', 10)
%             plot(V_fromG2toE(end),   LW_fromG2toE_new(end),   'k.', 'MarkerSize', 10)
%             % ---------------------------------------------------------------------
%             text(V_from0toSinv(end),  LW_from0toSinv_new(end), ' S inv.', 'FontSize', 6)
%             text(V_fromSinvtoG(end),  LW_fromSinvtoG_new(end), ' G',      'FontSize', 6)
%             text(V_fromGtoG1(end),    LW_fromGtoG1_new(end),   ' G1',     'FontSize', 6)
%             text(V_fromG1toF(end),    LW_fromG1toF_new(end),   ' F',      'FontSize', 6) 
%             text(V_fromFtoG2(end),    LW_fromFtoG2_new(end),   ' G2',     'FontSize', 6)  
%             text(V_fromG2toE(end),    LW_fromG2toE_new(end),   ' E',      'FontSize', 6) 
            % ---------------------------------------------------------------------
            plot(V_from0toSinv, WBL_from0toSinv, '-r', 'LineWidth', 1)
            plot(V_fromSinvtoG, WBL_fromSinvtoG, '-r', 'LineWidth', 1)
            plot(V_fromGtoG1,   WBL_fromGtoG1,   '-r', 'LineWidth', 1)
            plot(V_fromG1toF,   WBL_fromG1toF,   '-r', 'LineWidth', 1)
            plot(V_fromFtoG2,   WBL_fromFtoG2,   '-r', 'LineWidth', 1)
            plot(V_fromG2toE,   WBL_fromG2toE,   '-r', 'LineWidth', 1)
            plot(V_fromEto0,    WBL_fromEto0,    '-r', 'LineWidth', 1)
            % ---------------------------------------------------------------------
            plot(V_from0toSinv(1),   WBL_from0toSinv(1),   'k.', 'MarkerSize', 10)
            plot(V_from0toSinv(end), WBL_from0toSinv(end), 'k.', 'MarkerSize', 10)
            plot(V_fromSinvtoG(end), WBL_fromSinvtoG(end), 'k.', 'MarkerSize', 10)
            plot(V_fromGtoG1(end),   WBL_fromGtoG1(end),   'k.', 'MarkerSize', 10)
            plot(V_fromG1toF(end),   WBL_fromG1toF(end),   'k.', 'MarkerSize', 10)
            plot(V_fromFtoG2(end),   WBL_fromFtoG2(end),   'k.', 'MarkerSize', 10)
            plot(V_fromG2toE(end),   WBL_fromG2toE(end),   'k.', 'MarkerSize', 10)
            % ---------------------------------------------------------------------
            text(V_from0toSinv(end),  WBL_from0toSinv(end), ' S inv.', 'FontSize', 6)
            text(V_fromSinvtoG(end),  WBL_fromSinvtoG(end), ' G',      'FontSize', 6)
            text(V_fromGtoG1(end),    WBL_fromGtoG1(end),   ' G1',     'FontSize', 6)
            text(V_fromG1toF(end),    WBL_fromG1toF(end),   ' F',      'FontSize', 6) 
            text(V_fromFtoG2(end),    WBL_fromFtoG2(end),   ' G2',     'FontSize', 6)  
            text(V_fromG2toE(end),    WBL_fromG2toE(end),   ' E',      'FontSize', 6)             
            ylim padded;
            xlim padded;
            % STORE INSIDE THE AIRCRAFT STRUCTURE VARIABLE

            % POINT S_INV
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.value = q_from0toSinv(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.Attributes.unit = "Pa";    
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv.value = WBL_from0toSinv(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv.Attributes.unit = "daN";           
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv_new.value = WBL_from0toSinv_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LHTS_inv.value = LHT_from0toSinv(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LHTS_inv.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv.value = alfa_from0toSinv(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_rad.value = deg2rad(alfa_from0toSinv(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_rad.Attributes.unit = "rad";  
            alfaS_inv = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_new.value = alfa_from0toSinv_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_new_rad.value = deg2rad(alfa_from0toSinv_fullvehicle(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_new_rad.Attributes.unit = "rad";  
            alfaS_inv_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_new.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv.value = CL_from0toSinv(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv.Attributes.unit = "Non dimensional";  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv_new.value = CL_from0toSinv_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CLHT_S_inv.value = CLHT_from0toSinv(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CLHT_S_inv.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CD_S_inv.value = CD_from0toSinv(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CD_S_inv.Attributes.unit = "Non dimensional";  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CM_S_inv.value = CMCG_from0toSinv(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CM_S_inv.Attributes.unit = "Non dimensional";    
            % POINT G
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value = q_fromSinvtoG(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG.value = WBL_fromSinvtoG(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG.Attributes.unit = "daN";     
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG_new.value = WBL_fromSinvtoG_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.value = LHT_fromSinvtoG(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.value = alfa_fromSinvtoG(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_rad.value = deg2rad(alfa_fromSinvtoG(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_rad.Attributes.unit = "rad";
            alfaG = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_new.value = alfa_fromSinvtoG_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_new_rad.value = deg2rad(alfa_fromSinvtoG_fullvehicle(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_new_rad.Attributes.unit = "rad";
            alfaG_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_new.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value = CL_fromSinvtoG(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.Attributes.unit = "Non dimensional";  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G_new.value = CL_fromSinvtoG_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CLHT_G.value = CLHT_fromSinvtoG(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CLHT_G.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.value = CD_fromSinvtoG(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CM_G.value = CMCG_fromSinvtoG(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CM_G.Attributes.unit = "Non dimensional";
            % POINT G1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.qG1.value = q_fromGtoG1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.qG1.Attributes.unit = "Pa";  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1.value = WBL_fromGtoG1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1.Attributes.unit = "daN";          
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1_new.value = WBL_fromGtoG1_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LHTG1.value = LHT_fromGtoG1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LHTG1.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1.value = alfa_fromGtoG1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_rad.value = deg2rad(alfa_fromGtoG1(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_rad.Attributes.unit = "rad";
            alfaG1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_new.value = alfa_fromGtoG1_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_new_rad.value = deg2rad(alfa_fromGtoG1_fullvehicle(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_new_rad.Attributes.unit = "rad";
            alfaG1_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_new.value; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CL_G1.value = CL_fromGtoG1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CL_G1.Attributes.unit = "Non dimensional";  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CL_G1_new.value = CL_fromGtoG1_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CL_G1_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CLHT_G1.value = CLHT_fromGtoG1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CLHT_G1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CD_G1.value = CD_fromGtoG1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CD_G1.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CM_G1.value = CMCG_fromGtoG1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CM_G1.Attributes.unit = "Non dimensional";         
            % POINT F
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value = q_fromG1toF(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.Attributes.unit = "Pa";             
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.value = WBL_fromG1toF(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.Attributes.unit = "daN";     
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF_new.value = WBL_fromG1toF_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.value = LHT_fromG1toF(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.value = alfa_fromG1toF(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_rad.value = deg2rad(alfa_fromG1toF(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_rad.Attributes.unit = "rad";
            alfaF = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_new.value = alfa_fromG1toF_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_new_rad.value = deg2rad(alfa_fromG1toF_fullvehicle(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_new_rad.Attributes.unit = "rad";
            alfaF_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_new.value; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value = CL_fromG1toF(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.Attributes.unit = "Non dimensional";  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F_new.value = CL_fromG1toF_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CLHT_F.value = CLHT_fromG1toF(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CLHT_F.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.value = CD_fromG1toF(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CM_F.value = CMCG_fromG1toF(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CM_F.Attributes.unit = "Non dimensional";
            % POINT G2
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.qG2.value = q_fromFtoG2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.qG2.Attributes.unit = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.LG2.value = WBL_fromFtoG2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.LG2.Attributes.unit = "daN";          
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.LG2_new.value = WBL_fromFtoG2_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.LG2_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.LHTG2.value = LHT_fromFtoG2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.LHTG2.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.alfaG2.value = alfa_fromFtoG2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.alfaG2.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.alfaG2_rad.value = deg2rad(alfa_fromFtoG2(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.alfaG2_rad.Attributes.unit = "rad";
            alfaG2 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.alfaG2.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.alfaG2_new.value = alfa_fromFtoG2_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.alfaG2_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.alfaG2_new_rad.value = deg2rad(alfa_fromFtoG2_fullvehicle(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.alfaG2_new_rad.Attributes.unit = "rad";
            alfaG2_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.alfaG2_new.value; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CL_G2.value = CL_fromFtoG2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CL_G2.Attributes.unit = "Non dimensional";  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CL_G2_new.value = CL_fromFtoG2_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CL_G2_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CLHT_G2.value = CLHT_fromFtoG2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CLHT_G2.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CD_G2.value = CD_fromFtoG2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CD_G2.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CM_G2.value = CMCG_fromFtoG2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CM_G2.Attributes.unit = "Non dimensional";
            % POINT E
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value = q_fromG2toE(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.Attributes.unit = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE.value = WBL_fromG2toE(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE_new.value = WBL_fromG2toE_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.value = LHT_fromG2toE(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.value = alfa_fromG2toE(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_rad.value = deg2rad(alfa_fromG2toE(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_rad.Attributes.unit = "rad";
            alfaE = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_new.value = alfa_fromG2toE_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_new_rad.value = deg2rad(alfa_fromG2toE_fullvehicle(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_new_rad.Attributes.unit = "rad";
            alfaE_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_new.value;  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.value = CL_fromG2toE(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.Attributes.unit = "Non dimensional";  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E_new.value = CL_fromG2toE_fullvehicle(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CLHT_E.value = CLHT_fromG2toE(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CLHT_E.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CD_E.value = CD_fromG2toE(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CD_E.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CM_E.value = CMCG_fromG2toE(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CM_E.Attributes.unit = "Non dimensional";           
            
        elseif abs(min(n_gust_cruise_neg)) < abs(nmin)       
            % S_inv - G - F - E           
            % =============================================================
            % FROM 0 TO S_INVERTED 
            CL_from0toSinv   = zeros(length(V_from0toSinv), 1);
            alfa_from0toSinv = zeros(length(V_from0toSinv), 1);
            CD_from0toSinv   = zeros(length(V_from0toSinv), 1);
            q_from0toSinv    = zeros(length(V_from0toSinv), 1);
            WBL_from0toSinv  = zeros(length(V_from0toSinv), 1);
            WBL_from0toSinv_new  = zeros(length(V_from0toSinv), 1);
            CMCL_from0toSinv = zeros(length(V_from0toSinv), 1);
            CMCD_from0toSinv = zeros(length(V_from0toSinv), 1);
            CMCT_from0toSinv = zeros(length(V_from0toSinv), 1);
            CMCG_from0toSinv = zeros(length(V_from0toSinv), 1); 
            CLHT_from0toSinv = zeros(length(V_from0toSinv), 1);
            LHT_from0toSinv  = zeros(length(V_from0toSinv), 1);
            % =================================================================
            CL_from0toSinv_fullvehicle   = zeros(length(V_from0toSinv), 1);
            alfa_from0toSinv_fullvehicle = zeros(length(V_from0toSinv), 1);
            % ================================================================= 
            for i = 1:length(V_from0toSinv)
                CL_from0toSinv(i)   = CLmax_func(rho0, V_from0toSinv(i), WS, n_from0toSinv(i));
                alfa_from0toSinv(i) = alpha_fullmodel_inv(CL_from0toSinv(i), CL_max_inv, CL0, CLalfa, alfa_0l);
                CD_from0toSinv(i)   = polyval(p_CD_wb1, alfa_from0toSinv(i));
                q_from0toSinv(i)    = 0.5*rho0*(V_from0toSinv(i))^2;
                WBL_from0toSinv(i)  = q_from0toSinv(i)*S*CL_from0toSinv(i)*1e-1;
                CMCL_from0toSinv(i) = CLWB_contrib(obj1, CL_from0toSinv(i), deg2rad(alfa_from0toSinv(i)), XAC, XCG, bCG, MAC);
                CMCD_from0toSinv(i) = CDWB_contrib(obj1, CL_from0toSinv(i), deg2rad(alfa_from0toSinv(i)), XAC, XCG, bCG, MAC);
                CMCT_from0toSinv(i) = CT_contr(obj1, CD_from0toSinv(i), Thrust_axes, MAC);
                CMCG_from0toSinv(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_from0toSinv(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_from0toSinv(i) = CL_Tail(obj1, CMCL_from0toSinv(i), CMCD_from0toSinv(i), ...
                                                 CMCT_from0toSinv(i), CMCG_from0toSinv(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_from0toSinv(i)));
                % HORIZONTAL TAIL LIFT
                LHT_from0toSinv(i) = (0.5)*(V_from0toSinv(i)^2)*(S)*(rho0)*(CLHT_from0toSinv(i))*(1e-1);

                % DRAFT VERSION OF ITERATION 
                CL_tail            = CLHT_from0toSinv(i);
                CL_wb              = CL_from0toSinv(i);
                CL_from0toSinv_sin = CL_wb - CL_tail;
                tol                = 1e-3; 
                n                  = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_from0toSinv = alpha_fullmodel_inv(CL_from0toSinv_sin, CL_max_inv, CL0, CLalfa, alfa_0l);
                   CD_wb                = polyval(p_CD_wb1, alfa_new_from0toSinv);
                   CMCL_new             = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_from0toSinv), XAC, XCG, bCG, MAC);
                   CMCD_new             = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_from0toSinv), XAC, XCG, bCG, MAC);
                   CMCT_new             = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new             = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new             = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_from0toSinv));
                   CL_tail              = CLHT_new;
                   CL_new_from0toSinv   = CL_from0toSinv(i) - CL_tail;
                   CL_wb                = CL_wb + (CL_new_from0toSinv - CL_wb) * 1e-1;
                   n                    = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_from0toSinv_fullvehicle(i)   = CL_new_from0toSinv; 
                CLHT_from0toSinv(i) = CL_tail;
                alfa_from0toSinv_fullvehicle(i) = alfa_new_from0toSinv;
%                 alfa_from0toSinv(i) = alfa_new_from0toSinv;
                WBL_from0toSinv(i)  = q_from0toSinv(i) * S * CL_from0toSinv(i) * 1e-1;
                WBL_from0toSinv_new(i)  = q_from0toSinv(i) * S * CL_from0toSinv_fullvehicle(i) * 1e-1;
                LHT_from0toSinv(i)  = (0.5)*(V_from0toSinv(i)^2)*(S)*(rho0)*(CLHT_from0toSinv(i))*(1e-1);
                CMCL_from0toSinv(i) = CMCL_new;
                CMCD_from0toSinv(i) = CMCD_new;
                CMCG_from0toSinv(i) = CMCG_new;                
                
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv.value = CL_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.value = CL_from0toSinv_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv.value = alfa_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_rad.value = deg2rad(alfa_from0toSinv);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_rad.Attributes.unit = "rad";     
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_new.value = alfa_from0toSinv_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_new_rad.value = deg2rad(alfa_from0toSinv_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_new_rad.Attributes.unit = "rad";     
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toSinv.value = CD_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toSinv.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toSinv.value = q_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toSinv.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toSinv.value = WBL_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toSinv.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.value = WBL_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_from0toSinv.value = CMCL_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_from0toSinv.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_from0toSinv.value = CMCD_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_from0toSinv.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_from0toSinv.value = CMCT_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_from0toSinv.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_from0toSinv.value = CMCG_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_from0toSinv.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_from0toSinv.value = CLHT_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_from0toSinv.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_from0toSinv.value = LHT_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_from0toSinv.Attributes.unit = "daN";
            % =============================================================
            % FROM S_INVERTED TO G
            CL_fromSinvtoG   = zeros(length(V_fromSinvtoG), 1);
            alfa_fromSinvtoG = zeros(length(V_fromSinvtoG), 1);
            CD_fromSinvtoG   = zeros(length(V_fromSinvtoG), 1);
            q_fromSinvtoG    = zeros(length(V_fromSinvtoG), 1);
            WBL_fromSinvtoG  = zeros(length(V_fromSinvtoG), 1);
            WBL_fromSinvtoG_new  = zeros(length(V_fromSinvtoG), 1);
            CMCL_fromSinvtoG = zeros(length(V_fromSinvtoG), 1);
            CMCD_fromSinvtoG = zeros(length(V_fromSinvtoG), 1);
            CMCT_fromSinvtoG = zeros(length(V_fromSinvtoG), 1);
            CMCG_fromSinvtoG = zeros(length(V_fromSinvtoG), 1); 
            CLHT_fromSinvtoG = zeros(length(V_fromSinvtoG), 1);
            LHT_fromSinvtoG  = zeros(length(V_fromSinvtoG), 1);
            % =================================================================
            CL_fromSinvtoG_fullvehicle   = zeros(length(V_fromSinvtoG), 1);
            alfa_fromSinvtoG_fullvehicle = zeros(length(V_fromSinvtoG), 1);
            % ================================================================= 
            for i = 1:length(V_fromSinvtoG)
                CL_fromSinvtoG(i)   = CLmax_func(rho0, V_fromSinvtoG(i), WS, n_fromSinvtoG(i));
                alfa_fromSinvtoG(i) = alpha_fullmodel_inv(CL_fromSinvtoG(i), CL_max_inv, CL0, CLalfa, alfa_0l);
                CD_fromSinvtoG(i)   = polyval(p_CD_wb1, alfa_fromSinvtoG(i));
                q_fromSinvtoG(i)    = 0.5*rho0*(V_fromSinvtoG(i))^2;
                WBL_fromSinvtoG(i)  = q_fromSinvtoG(i)*S*CL_fromSinvtoG(i)*1e-1;
                CMCL_fromSinvtoG(i) = CLWB_contrib(obj1, CL_fromSinvtoG(i), deg2rad(alfa_fromSinvtoG(i)), XAC, XCG, bCG, MAC);
                CMCD_fromSinvtoG(i) = CDWB_contrib(obj1, CL_fromSinvtoG(i), deg2rad(alfa_fromSinvtoG(i)), XAC, XCG, bCG, MAC);
                CMCT_fromSinvtoG(i) = CT_contr(obj1, CD_fromSinvtoG(i), Thrust_axes, MAC);
                CMCG_fromSinvtoG(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromSinvtoG(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromSinvtoG(i) = CL_Tail(obj1, CMCL_fromSinvtoG(i), CMCD_fromSinvtoG(i), ...
                                                 CMCT_fromSinvtoG(i), CMCG_fromSinvtoG(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromSinvtoG(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromSinvtoG(i) = (0.5)*(V_fromSinvtoG(i)^2)*(S)*(rho0)*(CLHT_fromSinvtoG(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail         = CLHT_fromSinvtoG(i);
                CL_wb           = CL_fromSinvtoG(i);
                CL_new_fromSinvtoG  = CL_wb - CL_tail;
                tol             = 1e-3; 
                n               = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromSinvtoG = alpha_fullmodel_inv(CL_new_fromSinvtoG, CL_max_inv, CL0, CLalfa, alfa_0l);
                   CD_wb                = polyval(p_CD_wb1, alfa_new_fromSinvtoG);
                   CMCL_new             = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromSinvtoG), XAC, XCG, bCG, MAC);
                   CMCD_new             = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromSinvtoG), XAC, XCG, bCG, MAC);
                   CMCT_new             = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new             = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new             = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromSinvtoG));
                   CL_tail              = CLHT_new;
                   CL_new_fromSinvtoG   = CL_fromSinvtoG(i) - CL_tail;
                   CL_wb                = CL_wb + (CL_new_fromSinvtoG - CL_wb) * 1e-1;
                   n                    = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromSinvtoG_fullvehicle(i)   = CL_new_fromSinvtoG; 
                CLHT_fromSinvtoG(i) = CL_tail;
                alfa_fromSinvtoG_fullvehicle(i) = alfa_new_fromSinvtoG;
%                 alfa_fromSinvtoG(i) = alfa_new_fromSinvtoG;
                WBL_fromSinvtoG(i)  = q_fromSinvtoG(i) * S * CL_fromSinvtoG(i) * 1e-1;
                WBL_fromSinvtoG_new(i)  = q_fromSinvtoG(i) * S * CL_fromSinvtoG_fullvehicle(i) * 1e-1;
                LHT_fromSinvtoG(i)  = (0.5)*(V_fromSinvtoG(i)^2)*(S)*(rho0)*(CLHT_fromSinvtoG(i))*(1e-1);
                CMCL_fromSinvtoG(i) = CMCL_new;
                CMCD_fromSinvtoG(i) = CMCD_new;
                CMCG_fromSinvtoG(i) = CMCG_new;     
                
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG.value = CL_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG_new.value = CL_fromSinvtoG_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG.value = alfa_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG_rad.value = deg2rad(alfa_fromSinvtoG);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG_rad.Attributes.unit = "rad"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG_new.value = alfa_fromSinvtoG_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG_new_rad.value = deg2rad(alfa_fromSinvtoG_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG_new_rad.Attributes.unit = "rad";       
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromSinvtoG.value = CD_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromSinvtoG.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromSinvtoG.value = q_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromSinvtoG.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromSinvtoG.value = WBL_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromSinvtoG.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromSinvtoG_new.value = WBL_fromSinvtoG_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromSinvtoG_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromSinvtoG.value = CMCL_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromSinvtoG.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromSinvtoG.value = CMCD_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromSinvtoG.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromSinvtoG.value = CMCT_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromSinvtoG.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromSinvtoG.value = CMCG_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromSinvtoG.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromSinvtoG.value = CLHT_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromSinvtoG.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromSinvtoG.value = LHT_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromSinvtoG.Attributes.unit = "daN";
            % =============================================================
            % FROM G TO F
            CL_fromGtoF   = zeros(length(V_fromGtoF), 1);
            alfa_fromGtoF = zeros(length(V_fromGtoF), 1);
            CD_fromGtoF   = zeros(length(V_fromGtoF), 1);
            q_fromGtoF    = zeros(length(V_fromGtoF), 1);
            WBL_fromGtoF  = zeros(length(V_fromGtoF), 1);
            WBL_fromGtoF_new  = zeros(length(V_fromGtoF), 1);
            CMCL_fromGtoF = zeros(length(V_fromGtoF), 1);
            CMCD_fromGtoF = zeros(length(V_fromGtoF), 1);
            CMCT_fromGtoF = zeros(length(V_fromGtoF), 1);
            CMCG_fromGtoF = zeros(length(V_fromGtoF), 1); 
            CLHT_fromGtoF = zeros(length(V_fromGtoF), 1);
            LHT_fromGtoF  = zeros(length(V_fromGtoF), 1);
            % =================================================================
            CL_fromGtoF_fullvehicle   = zeros(length(V_fromGtoF), 1);
            alfa_fromGtoF_fullvehicle = zeros(length(V_fromGtoF), 1);
            % ================================================================= 
            for i = 1:length(V_fromGtoF)
                CL_fromGtoF(i)   = CLmax_func(rho0, V_fromGtoF(i), WS, n_fromGtoF(i));
                alfa_fromGtoF(i) = alpha_fullmodel_inv(CL_fromGtoF(i), CL_max_inv, CL0, CLalfa, alfa_0l);
                CD_fromGtoF(i)   = polyval(p_CD_wb1, alfa_fromGtoF(i));
                q_fromGtoF(i)    = 0.5*rho0*(V_fromGtoF(i))^2;
                WBL_fromGtoF(i)  = q_fromGtoF(i)*S*CL_fromGtoF(i)*1e-1;
                CMCL_fromGtoF(i) = CLWB_contrib(obj1, CL_fromGtoF(i), deg2rad(alfa_fromGtoF(i)), XAC, XCG, bCG, MAC);
                CMCD_fromGtoF(i) = CDWB_contrib(obj1, CL_fromGtoF(i), deg2rad(alfa_fromGtoF(i)), XAC, XCG, bCG, MAC);
                CMCT_fromGtoF(i) = CT_contr(obj1, CD_fromGtoF(i), Thrust_axes, MAC);
                CMCG_fromGtoF(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromGtoF(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromGtoF(i) = CL_Tail(obj1, CMCL_fromGtoF(i), CMCD_fromGtoF(i), ...
                                                 CMCT_fromGtoF(i), CMCG_fromGtoF(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromGtoF(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromGtoF(i) = (0.5)*(V_fromGtoF(i)^2)*(S)*(rho0)*(CLHT_fromGtoF(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail         = CLHT_fromGtoF(i);
                CL_wb           = CL_fromGtoF(i);
                CL_new_fromGtoF = CL_wb - CL_tail;
                tol             = 1e-3; 
                n               = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromSinvtoG = alpha_fullmodel_inv(CL_new_fromGtoF, CL_max_inv, CL0, CLalfa, alfa_0l);
                   CD_fromGtoF(i)       = polyval(p_CD_wb1, alfa_new_fromSinvtoG);
                   CMCL_new             = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromGtoF), XAC, XCG, bCG, MAC);
                   CMCD_new             = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromGtoF), XAC, XCG, bCG, MAC);
                   CMCT_new             = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new             = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new             = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromGtoF));
                   CL_tail              = CLHT_new;
                   CL_new_fromGtoF      = CL_fromGtoF(i) - CL_tail;
                   CL_wb                = CL_wb + (CL_new_fromGtoF - CL_wb) * 1e-1;
                   n                    = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromGtoF_fullvehicle(i)   = CL_new_fromGtoF; 
                CLHT_fromGtoF(i) = CL_tail;
                alfa_fromGtoF_fullvehicle(i) = alfa_new_fromGtoF;
%                 alfa_fromGtoF(i) = alfa_new_fromGtoF;
                WBL_fromGtoF(i)  = q_fromGtoF(i) * S * CL_fromGtoF(i) * 1e-1;
                WBL_fromGtoF_new(i)  = q_fromGtoF(i) * S * CL_fromGtoF_fullvehicle(i) * 1e-1;
                LHT_fromGtoF(i)  = (0.5)*(V_fromGtoF(i)^2)*(S)*(rho0)*(CLHT_fromGtoF(i))*(1e-1);
                CMCL_fromGtoF(i) = CMCL_new;
                CMCD_fromGtoF(i) = CMCD_new;
                CMCG_fromGtoF(i) = CMCG_new;   
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoF.value = CL_fromGtoF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoF.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoF_new.value = CL_fromGtoF_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoF_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoF.value = alfa_fromGtoF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoF.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoF_rad.value = deg2rad(alfa_fromGtoF);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoF_rad.Attributes.unit = "rad";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoF_new.value = alfa_fromGtoF_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoF_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoF_new_rad.value = deg2rad(alfa_fromGtoF_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoF_new_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromGtoF.value = CD_fromGtoF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromGtoF.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromGtoF.value = q_fromGtoF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromGtoF.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromGtoF.value = WBL_fromGtoF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromGtoF.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromGtoF_new.value = WBL_fromGtoF_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromGtoF_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromGtoF.value = CMCL_fromGtoF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromGtoF.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromGtoF.value = CMCD_fromGtoF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromGtoF.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromGtoF.value = CMCT_fromGtoF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromGtoF.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromGtoF.value = CMCG_fromGtoF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromGtoF.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromGtoF.value = CLHT_fromGtoF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromGtoF.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromGtoF.value = LHT_fromGtoF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromGtoF.Attributes.unit = "daN";
            % =============================================================
            % FROM F TO E 
            CL_fromFtoE   = zeros(length(V_fromFtoE), 1);
            alfa_fromFtoE = zeros(length(V_fromFtoE), 1);
            CD_fromFtoE   = zeros(length(V_fromFtoE), 1);
            q_fromFtoE    = zeros(length(V_fromFtoE), 1);
            WBL_fromFtoE  = zeros(length(V_fromFtoE), 1);
            WBL_fromFtoE_new  = zeros(length(V_fromFtoE), 1);
            CMCL_fromFtoE = zeros(length(V_fromFtoE), 1);
            CMCD_fromFtoE = zeros(length(V_fromFtoE), 1);
            CMCT_fromFtoE = zeros(length(V_fromFtoE), 1);
            CMCG_fromFtoE = zeros(length(V_fromFtoE), 1); 
            CLHT_fromFtoE = zeros(length(V_fromFtoE), 1);
            LHT_fromFtoE  = zeros(length(V_fromFtoE), 1);
            % =================================================================
            CL_fromFtoE_fullvehicle   = zeros(length(V_fromFtoE), 1);
            alfa_fromFtoE_fullvehicle = zeros(length(V_fromFtoE), 1);
            % ================================================================= 
            for i = 1:length(V_fromFtoE)
                CL_fromFtoE(i)   = CLmax_func(rho0, V_fromFtoE(i), WS, n_fromFtoE(i));
                alfa_fromFtoE(i) = alpha_fullmodel_inv(CL_fromFtoE(i), CL_max_inv, CL0, CLalfa, alfa_0l);
                CD_fromFtoE(i)   = polyval(p_CD_wb1, alfa_fromFtoE(i));
                q_fromFtoE(i)    = 0.5*rho0*(V_fromFtoE(i))^2;
                WBL_fromFtoE(i)  = q_fromFtoE(i)*S*CL_fromFtoE(i)*1e-1;
                CMCL_fromFtoE(i) = CLWB_contrib(obj1, CL_fromFtoE(i), deg2rad(alfa_fromFtoE(i)), XAC, XCG, bCG, MAC);
                CMCD_fromFtoE(i) = CDWB_contrib(obj1, CL_fromFtoE(i), deg2rad(alfa_fromFtoE(i)), XAC, XCG, bCG, MAC);
                CMCT_fromFtoE(i) = CT_contr(obj1, CD_fromFtoE(i), Thrust_axes, MAC);
                CMCG_fromFtoE(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromFtoE(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromFtoE(i) = CL_Tail(obj1, CMCL_fromFtoE(i), CMCD_fromFtoE(i), ...
                                                 CMCT_fromFtoE(i), CMCG_fromFtoE(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromFtoE(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromFtoE(i) = (0.5)*(V_fromFtoE(i)^2)*(S)*(rho0)*(CLHT_fromFtoE(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail         = CLHT_fromFtoE(i);
                CL_wb           = CL_fromFtoE(i);
                CL_new_fromFtoE = CL_wb - CL_tail;
                tol             = 1e-3; 
                n               = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromFtoE = alpha_fullmodel_inv(CL_new_fromFtoE, CL_max_inv, CL0, CLalfa, alfa_0l);
                   CD_wb             = polyval(p_CD_wb1, alfa_new_fromFtoE);
                   CMCL_new          = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromFtoE), XAC, XCG, bCG, MAC);
                   CMCD_new          = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromFtoE), XAC, XCG, bCG, MAC);
                   CMCT_new          = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new          = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new          = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromFtoE));
                   CL_tail           = CLHT_new;
                   CL_new_fromFtoE   = CL_fromFtoE(i) - CL_tail;
                   CL_wb             = CL_wb + (CL_new_fromFtoE - CL_wb) * 1e-1;
                   n                 = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromFtoE_fullvehicle(i)   = CL_new_fromFtoE; 
                CLHT_fromFtoE(i) = CL_tail;
                alfa_fromFtoE_fullvehicle(i) = alfa_new_fromFtoE;
%                 alfa_fromFtoE(i) = alfa_new_fromFtoE;
                WBL_fromFtoE(i)  = q_fromFtoE(i) * S * CL_fromFtoE(i) * 1e-1;
                WBL_fromFtoE_new(i)  = q_fromFtoE(i) * S * CL_fromFtoE_fullvehicle(i) * 1e-1;
                LHT_fromFtoE(i)  = (0.5)*(V_fromFtoE(i)^2)*(S)*(rho0)*(CLHT_fromFtoE(i))*(1e-1);
                CMCL_fromFtoE(i) = CMCL_new;
                CMCD_fromFtoE(i) = CMCD_new;
                CMCG_fromFtoE(i) = CMCG_new;  
                
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value = CL_fromFtoE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value = CL_fromFtoE_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE.value = alfa_fromFtoE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE_rad.value = deg2rad(alfa_fromFtoE);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE_rad.Attributes.unit = "rad";       
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE_new.value = alfa_fromFtoE_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE_new_rad.value = deg2rad(alfa_fromFtoE_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE_new_rad.Attributes.unit = "rad";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value = CD_fromFtoE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromFtoE.value = q_fromFtoE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromFtoE.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE.value = WBL_fromFtoE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoE_new.value = WBL_fromFtoE_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoE_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromFtoE.value = CMCL_fromFtoE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromFtoE.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromFtoE.value = CMCD_fromFtoE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromFtoE.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromFtoE.value = CMCT_fromFtoE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromFtoE.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromFtoE.value = CMCG_fromFtoE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromFtoE.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromFtoE.value = CLHT_fromFtoE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromFtoE.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromFtoE.value = LHT_fromFtoE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromFtoE.Attributes.unit = "daN";   
            % =============================================================
            % FROM E TO 0
            CL_fromEto0   = zeros(length(V_fromEto0), 1);
            alfa_fromEto0 = zeros(length(V_fromEto0), 1);
            CD_fromEto0   = zeros(length(V_fromEto0), 1);
            q_fromEto0    = zeros(length(V_fromEto0), 1);
            WBL_fromEto0  = zeros(length(V_fromEto0), 1);
            WBL_fromEto0_new  = zeros(length(V_fromEto0), 1);
            CMCL_fromEto0 = zeros(length(V_fromEto0), 1);
            CMCD_fromEto0 = zeros(length(V_fromEto0), 1);
            CMCT_fromEto0 = zeros(length(V_fromEto0), 1);
            CMCG_fromEto0 = zeros(length(V_fromEto0), 1); 
            CLHT_fromEto0 = zeros(length(V_fromEto0), 1);
            LHT_fromEto0  = zeros(length(V_fromEto0), 1);
            % =================================================================
            CL_fromEto0_fullvehicle   = zeros(length(V_fromEto0), 1);
            alfa_fromEto0_fullvehicle = zeros(length(V_fromEto0), 1);
            % ================================================================= 
            for i = 1:length(V_fromEto0)
                CL_fromEto0(i)   = CLmax_func(rho0, V_fromEto0(i), WS, n_fromEto0(i));
                alfa_fromEto0(i) = alpha_fullmodel_inv(CL_fromEto0(i), CL_max_inv, CL0, CLalfa, alfa_0l);
                CD_fromEto0(i)   = polyval(p_CD_wb1, alfa_fromEto0(i));
                q_fromEto0(i)    = 0.5*rho0*(V_fromEto0(i))^2;
                WBL_fromEto0(i)  = q_fromEto0(i)*S*CL_fromEto0(i)*1e-1;
                CMCL_fromEto0(i) = CLWB_contrib(obj1, CL_fromEto0(i), deg2rad(alfa_fromEto0(i)), XAC, XCG, bCG, MAC);
                CMCD_fromEto0(i) = CDWB_contrib(obj1, CL_fromEto0(i), deg2rad(alfa_fromEto0(i)), XAC, XCG, bCG, MAC);
                CMCT_fromEto0(i) = CT_contr(obj1, CD_fromEto0(i), Thrust_axes, MAC);
                CMCG_fromEto0(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromEto0(i));
                % HORIZONTAL TAIL LIFT COEFFICIENT
                CLHT_fromEto0(i) = CL_Tail(obj1, CMCL_fromEto0(i), CMCD_fromEto0(i), ...
                                                 CMCT_fromEto0(i), CMCG_fromEto0(i), ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_fromEto0(i)));
                % HORIZONTAL TAIL LIFT
                LHT_fromEto0(i) = (0.5)*(V_fromEto0(i)^2)*(S)*(rho0)*(CLHT_fromEto0(i))*(1e-1);
                
                % DRAFT VERSION OF ITERATION 
                CL_tail      = CLHT_fromEto0(i);
                CL_wb        = CL_fromEto0(i);
                CL_new_fromEto0  = CL_wb - CL_tail;
                tol          = 1e-3; 
                n            = 1;
                while abs(CL_wb - CL_tail) > tol
                   alfa_new_fromEto0 = alpha_fullmodel_inv(CL_new_fromEto0, CL_max_inv, CL0, CLalfa, alfa_0l);
                   CD_wb             = polyval(p_CD_wb1, alfa_new_fromEto0);
                   CMCL_new          = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromEto0), XAC, XCG, bCG, MAC);
                   CMCD_new          = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromEto0), XAC, XCG, bCG, MAC);
                   CMCT_new          = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
                   CMCG_new          = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
                   CLHT_new          = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                                 l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromEto0));
                   CL_tail           = CLHT_new;
                   CL_new_fromEto0   = CL_fromEto0(i) - CL_tail;
                   CL_wb             = CL_wb + (CL_new_fromEto0 - CL_wb) * 1e-1;
                   n                 = n + 1;
                   if n == 15
                       break
                   end
                end
                CL_fromEto0_fullvehicle(i)   = CL_new_fromEto0; 
                CLHT_fromEto0(i) = CL_tail;
                alfa_fromEto0_fullvehicle(i) = alfa_new_fromEto0;
%                 alfa_fromEto0(i) = alfa_new_fromEto0;
                WBL_fromEto0(i)  = q_fromEto0(i) * S * CL_fromEto0(i) * 1e-1;
                WBL_fromEto0_new(i)  = q_fromEto0(i) * S * CL_fromEto0_fullvehicle(i) * 1e-1;
                LHT_fromEto0(i)  = (0.5)*(V_fromEto0(i)^2)*(S)*(rho0)*(CLHT_fromEto0(i))*(1e-1);
                CMCL_fromEto0(i) = CMCL_new;
                CMCD_fromEto0(i) = CMCD_new;
                CMCG_fromEto0(i) = CMCG_new;
                
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0.value = CL_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.value = CL_fromEto0_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0.value = alfa_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_rad.value = deg2rad(alfa_fromEto0);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_rad.Attributes.unit = "rad";    
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_new.value = alfa_fromEto0_fullvehicle;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_new.Attributes.unit = "deg";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_new_rad.value = deg2rad(alfa_fromEto0_fullvehicle);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_new_rad.Attributes.unit = "rad";      
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromEto0.value = CD_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromEto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromEto0.value = q_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromEto0.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromEto0.value = WBL_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromEto0.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.value = WBL_fromEto0_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromEto0.value = CMCL_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromEto0.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromEto0.value = CMCD_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromEto0.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromEto0.value = CMCT_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromEto0.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromEto0.value = CMCG_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromEto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromEto0.value = CLHT_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromEto0.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromEto0.value = LHT_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromEto0.Attributes.unit = "daN"; 
            % =================================================================
            % TAIL LOADS DIAGRAM - NEGATIVE SIDE 
            figure(5); 
            % hold on; grid on; grid minor; 
            % ---------------------------------------------------------------------
            plot(V_from0toSinv, LHT_from0toSinv, '-r', 'LineWidth', 1)
            plot(V_fromSinvtoG, LHT_fromSinvtoG, '-r', 'LineWidth', 1)
            plot(V_fromGtoF,    LHT_fromGtoF,    '-r', 'LineWidth', 1)
            plot(V_fromFtoE,    LHT_fromFtoE,    '-r', 'LineWidth', 1)
            plot(V_fromEto0,    LHT_fromEto0,    '-r', 'LineWidth', 1)
            % ---------------------------------------------------------------------
            plot(V_from0toSinv(1),   LHT_from0toSinv(1),   'k.', 'MarkerSize', 10)
            plot(V_from0toSinv(end), LHT_from0toSinv(end), 'k.', 'MarkerSize', 10)
            plot(V_fromSinvtoG(end), LHT_fromSinvtoG(end), 'k.', 'MarkerSize', 10)
            plot(V_fromGtoF(end),    LHT_fromGtoF(end),    'k.', 'MarkerSize', 10)
            plot(V_fromFtoE(end),    LHT_fromFtoE(end),    'k.', 'MarkerSize', 10)
            % ---------------------------------------------------------------------
            text(V_from0toSinv(end), LHT_from0toSinv(end), ' S inv.', 'FontSize', 6)
            text(V_fromSinvtoG(end), LHT_fromSinvtoG(end), ' G',      'FontSize', 6)
            text(V_fromGtoF(end),    LHT_fromGtoF(end),    ' F',      'FontSize', 6)
            text(V_fromFtoE(end),    LHT_fromFtoE(end),    ' E',      'FontSize', 6)
            % =================================================================
            ylim padded;
            xlim padded;
% %             % MAIN WING LOADS DIAGRAM        
% %             CL_from0toSinv_new = CL_from0toSinv + CLHT_from0toSinv;        
% %             CL_fromSinvtoG_new = CL_fromSinvtoG + CLHT_fromSinvtoG;        
% %             CL_fromGtoF_new    = CL_fromGtoF    + CLHT_fromGtoF;        
% %             CL_fromFtoE_new    = CL_fromFtoE    + CLHT_fromFtoE;       
% %             CL_fromEto0_new    = CL_fromEto0    + CLHT_fromEto0;
%             % MAIN WING LOADS DIAGRAM        
%             CL_from0toSinv_new = CL_from0toSinv - CLHT_from0toSinv;        
%             CL_fromSinvtoG_new = CL_fromSinvtoG - CLHT_fromSinvtoG;        
%             CL_fromGtoF_new    = CL_fromGtoF    - CLHT_fromGtoF;        
%             CL_fromFtoE_new    = CL_fromFtoE    - CLHT_fromFtoE;       
%             CL_fromEto0_new    = CL_fromEto0    - CLHT_fromEto0;
%             LW_from0toSinv_new = zeros(length(CL_from0toSinv_new), 1);
%             LW_fromSinvtoG_new = zeros(length(CL_from0toSinv_new), 1);
%             LW_fromGtoF_new    = zeros(length(CL_from0toSinv_new), 1);
%             LW_fromFtoE_new    = zeros(length(CL_from0toSinv_new), 1);
%             LW_fromEto0_new    = zeros(length(CL_from0toSinv_new), 1);
%             for i = 1:length(CL_from0toSinv_new)
%                 LW_from0toSinv_new(i)  = (0.5)*(V_from0toSinv(i)^2)* rho0 * S *(CL_from0toSinv_new(i))*(1e-1);  
%                 LW_fromSinvtoG_new(i) = (0.5)*(V_fromSinvtoG(i)^2)* rho0 * S *(CL_fromSinvtoG_new(i))*(1e-1);
%                 LW_fromGtoF_new(i)    = (0.5)*(V_fromGtoF(i)^2)* rho0 * S *(CL_fromGtoF_new(i))*(1e-1);
%                 LW_fromFtoE_new(i)     = (0.5)*(V_fromFtoE(i)^2)* rho0 * S *(CL_fromFtoE_new(i))*(1e-1);
%                 LW_fromEto0_new(i)     = (0.5)*(V_fromEto0(i)^2)* rho0 * S *(CL_fromEto0_new(i))*(1e-1);
% %                 LW_from0toSinv_new(i)  = (0.5)*(V_from0toSinv(i)^2)* rho0 * S *(CL_from0toSinv(i))*(1e-1);  
% %                 LW_fromSinvtoG_new(i) = (0.5)*(V_fromSinvtoG(i)^2)* rho0 * S *(CL_fromSinvtoG(i))*(1e-1);
% %                 LW_fromGtoF_new(i)    = (0.5)*(V_fromGtoF(i)^2)* rho0 * S *(CL_fromGtoF(i))*(1e-1);
% %                 LW_fromFtoE_new(i)     = (0.5)*(V_fromFtoE(i)^2)* rho0 * S *(CL_fromFtoE(i))*(1e-1);
% %                 LW_fromEto0_new(i)     = (0.5)*(V_fromEto0(i)^2)* rho0 * S *(CL_fromEto0(i))*(1e-1);
%             end       
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.value = CL_from0toSinv_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG_new.value = CL_fromSinvtoG_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoF_new.value = CL_fromGtoF_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoF_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value = CL_fromFtoE_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.Attributes.unit = "Non dimensional";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.value = CL_fromEto0_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.Attributes.unit = "Non dimensional";
%             % ------------------------------------------------------------------------------------------------------------------------
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.value = LW_from0toSinv_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromSinvtoG_new.value = LW_fromSinvtoG_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromSinvtoG_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromGtoF_new.value = LW_fromGtoF_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromGtoF_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoE_new.value = LW_fromFtoE_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoE_new.Attributes.unit = "daN";
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.value = LW_fromEto0_new;
%             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.Attributes.unit = "daN"; 
            % =============================================================
            % WING LOADS DIAGRAM - NEGATIVE SIDE   
            figure(6)
%             % ---------------------------------------------------------------------
%             plot(V_from0toSinv, LW_from0toSinv_new, '-r', 'LineWidth', 1)
%             plot(V_fromSinvtoG, LW_fromSinvtoG_new, '-r', 'LineWidth', 1)
%             plot(V_fromGtoF,    LW_fromGtoF_new,    '-r', 'LineWidth', 1)
%             plot(V_fromFtoE,    LW_fromFtoE_new,    '-r', 'LineWidth', 1)
%             plot(V_fromEto0,    LW_fromEto0_new,    '-r', 'LineWidth', 1)
%             % ---------------------------------------------------------------------
%             plot(V_from0toSinv(1),   LW_from0toSinv_new(1),   'k.', 'MarkerSize', 10)
%             plot(V_from0toSinv(end), LW_from0toSinv_new(end), 'k.', 'MarkerSize', 10)
%             plot(V_fromSinvtoG(end), LW_fromSinvtoG_new(end), 'k.', 'MarkerSize', 10)
%             plot(V_fromGtoF(end),    LW_fromGtoF_new(end),    'k.', 'MarkerSize', 10)
%             plot(V_fromFtoE(end),    LW_fromFtoE_new(end),    'k.', 'MarkerSize', 10)
%             plot(V_fromEto0(end),    LW_fromEto0_new(end),    'k.', 'MarkerSize', 10)
%             % ---------------------------------------------------------------------
%             text(V_from0toSinv(end),  LW_from0toSinv_new(end), ' S inv.', 'FontSize', 6)
%             text(V_fromSinvtoG(end), LW_fromSinvtoG_new(end),  ' G',      'FontSize', 6)
%             text(V_fromGtoF(end),    LW_fromGtoF_new(end),     ' F',      'FontSize', 6)
%             text(V_fromFtoE(end),    LW_fromFtoE_new(end),     ' E',      'FontSize', 6)  

            % ---------------------------------------------------------------------
            plot(V_from0toSinv, WBL_from0toSinv, '-r', 'LineWidth', 1)
            plot(V_fromSinvtoG, WBL_fromSinvtoG, '-r', 'LineWidth', 1)
            plot(V_fromGtoF,    WBL_fromGtoF,    '-r', 'LineWidth', 1)
            plot(V_fromFtoE,    WBL_fromFtoE,    '-r', 'LineWidth', 1)
            plot(V_fromEto0,    WBL_fromEto0,    '-r', 'LineWidth', 1)
            % ---------------------------------------------------------------------
            plot(V_from0toSinv(1),   WBL_from0toSinv(1),   'k.', 'MarkerSize', 10)
            plot(V_from0toSinv(end), WBL_from0toSinv(end), 'k.', 'MarkerSize', 10)
            plot(V_fromSinvtoG(end), WBL_fromSinvtoG(end), 'k.', 'MarkerSize', 10)
            plot(V_fromGtoF(end),    WBL_fromGtoF(end),    'k.', 'MarkerSize', 10)
            plot(V_fromFtoE(end),    WBL_fromFtoE(end),    'k.', 'MarkerSize', 10)
            plot(V_fromEto0(end),    WBL_fromEto0(end),    'k.', 'MarkerSize', 10)
            % ---------------------------------------------------------------------
            text(V_from0toSinv(end), WBL_from0toSinv(end),  ' S inv.', 'FontSize', 6)
            text(V_fromSinvtoG(end), WBL_fromSinvtoG(end),  ' G',      'FontSize', 6)
            text(V_fromGtoF(end),    WBL_fromGtoF(end),     ' F',      'FontSize', 6)
            text(V_fromFtoE(end),    WBL_fromFtoE(end),     ' E',      'FontSize', 6)      
            ylim padded;
            xlim padded;
        % STORE INSIDE THE AIRCRAFT STRUCTURE VARIABLE
        
        % POINT S_INV
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.value = q_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.Attributes.unit = "Pa";    
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv.value = WBL_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv.Attributes.unit = "daN";           
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv_new.value = WBL_from0toSinv_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LHTS_inv.value = LHT_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LHTS_inv.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv.value = alfa_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_rad.value = deg2rad(alfa_from0toSinv(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_rad.Attributes.unit = "rad";  
        alfaS_inv = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_new.value = alfa_from0toSinv_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_new_rad.value = deg2rad(alfa_from0toSinv_fullvehicle(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_new_rad.Attributes.unit = "rad";  
        alfaS_inv_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_new.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv.value = CL_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv.Attributes.unit = "Non dimensional";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv_new.value = CL_from0toSinv_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv_new.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CLHT_S_inv.value = CLHT_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CLHT_S_inv.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CD_S_inv.value = CD_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CD_S_inv.Attributes.unit = "Non dimensional";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CM_S_inv.value = CMCG_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CM_S_inv.Attributes.unit = "Non dimensional"; 
        % POINT G
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value = q_fromSinvtoG(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG.value = WBL_fromSinvtoG(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG.Attributes.unit = "daN";     
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG_new.value = WBL_fromSinvtoG_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.value = LHT_fromSinvtoG(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.value = alfa_fromSinvtoG(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_rad.value = deg2rad(alfa_fromSinvtoG(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_rad.Attributes.unit = "rad";
        alfaG = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_new.value = alfa_fromSinvtoG_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_new_rad.value = deg2rad(alfa_fromSinvtoG_fullvehicle(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_new_rad.Attributes.unit = "rad";
        alfaG_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_new.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value = CL_fromSinvtoG(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.Attributes.unit = "Non dimensional";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G_new.value = CL_fromSinvtoG_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G_new.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CLHT_G.value = CLHT_fromSinvtoG(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CLHT_G.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.value = CD_fromSinvtoG(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CM_G.value = CMCG_fromSinvtoG(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CM_G.Attributes.unit = "Non dimensional";
        % POINT F
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value = q_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.Attributes.unit = "Pa";             
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.value = WBL_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.Attributes.unit = "daN";     
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF_new.value = WBL_fromG1toF_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.value = LHT_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.value = alfa_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_rad.value = deg2rad(alfa_fromG1toF(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_rad.Attributes.unit = "rad";
        alfaF = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaG1.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_new.value = alfa_fromG1toF_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_new_rad.value = deg2rad(alfa_fromG1toF_fullvehicle(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_new_rad.Attributes.unit = "rad";
        alfaF_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_new.value; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value = CL_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.Attributes.unit = "Non dimensional";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F_new.value = CL_fromG1toF_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F_new.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CLHT_F.value = CLHT_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CLHT_F.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.value = CD_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CM_F.value = CMCG_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CM_F.Attributes.unit = "Non dimensional";
        % POINT E
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value = q_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.Attributes.unit = "Pa"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE.value = WBL_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE_new.value = WBL_fromFtoE_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.value = LHT_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.value = alfa_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_rad.value = deg2rad(alfa_fromFtoE(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_rad.Attributes.unit = "rad";
        alfaE = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_new.value = alfa_fromFtoE_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_new_rad.value = deg2rad(alfa_fromFtoE_fullvehicle(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_new_rad.Attributes.unit = "rad";
        alfaE_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_new.value;  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.value = CL_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.Attributes.unit = "Non dimensional";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E_new.value = CL_fromFtoE_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E_new.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CLHT_E.value = CLHT_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CLHT_E.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CD_E.value = CD_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CD_E.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CM_E.value = CMCG_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CM_E.Attributes.unit = "Non dimensional";  
        
        end
        
    % CASE 2: Real solutions of the intercept     
    case 'Case 2'
        % S_inv - G1 - F - E 
        % CL, ALFA, CD CALCULATIONS - POSITIVE LOAD FACTORS
        % ----------------------------------------------------------------- 
        % Calculation of the CL for the wing - body configuration. It must be
        % noticed that the calculations relative to lift coefficient in this
        % particular range of airspeed V and load factor n are going to give a
        % constant as a results, which is the maximum lift coefficient following
        % the lift curve of the aircrat. The previously defined CLMAX is a little
        % bit greater than the one obtained from V - n diagram flight condition to
        % take into account the extra lift produced by the horizontal empennage.
        % WING BODY LIFT CALCULATION
        % ----------------------------------------------------------------- 
        % In this section of the code the following formula is applied to evaluate
        % the wing-body lift with the formula 
        %        
        %   L = 0.5*rho*(V^2)*Sw*CL 
        %
        % using the previously calculated lift coefficients, coming from the 
        % various curves of the final envelope diagram.
        % =================================================================       
        % FROM 0 TO S_INVERTED 
        CL_from0toSinv      = zeros(length(V_from0toSinv), 1);
        alfa_from0toSinv    = zeros(length(V_from0toSinv), 1);
        CD_from0toSinv      = zeros(length(V_from0toSinv), 1);
        q_from0toSinv       = zeros(length(V_from0toSinv), 1);
        WBL_from0toSinv     = zeros(length(V_from0toSinv), 1);
        WBL_from0toSinv_new = zeros(length(V_from0toSinv), 1);
        CMCL_from0toSinv    = zeros(length(V_from0toSinv), 1);
        CMCD_from0toSinv    = zeros(length(V_from0toSinv), 1);
        CMCT_from0toSinv    = zeros(length(V_from0toSinv), 1);
        CMCG_from0toSinv    = zeros(length(V_from0toSinv), 1); 
        CLHT_from0toSinv    = zeros(length(V_from0toSinv), 1);
        LHT_from0toSinv     = zeros(length(V_from0toSinv), 1);
        % =================================================================
        CL_from0toSinv_fullvehicle   = zeros(length(V_from0toSinv), 1);
        alfa_from0toSinv_fullvehicle = zeros(length(V_from0toSinv), 1);
        % ================================================================= 
        for i = 1:length(V_from0toSinv)
            CL_from0toSinv(i)   = CLmax_func(rho0, V_from0toSinv(i), WS, n_from0toSinv(i));
            if CL_from0toSinv(i) < CL_max_inv
                CL_from0toSinv(i) = CL_from0toSinv(i);
            elseif CL_from0toSinv(i) > CL_max_inv
                CL_from0toSinv(i) = CL_max_inv;
            end
            alfa_from0toSinv(i) = alpha_fullmodel_inv(CL_from0toSinv(i), CL_max_inv, CL0, CLalfa, alfa_0l);
            CD_from0toSinv(i)   = polyval(p_CD_wb1, alfa_from0toSinv(i));
            q_from0toSinv(i)    = 0.5*rho0*(V_from0toSinv(i))^2;
            WBL_from0toSinv(i)  = q_from0toSinv(i)*S*CL_from0toSinv(i)*1e-1;
            CMCL_from0toSinv(i) = CLWB_contrib(obj1, CL_from0toSinv(i), deg2rad(alfa_from0toSinv(i)), XAC, XCG, bCG, MAC);
            CMCD_from0toSinv(i) = CDWB_contrib(obj1, CL_from0toSinv(i), deg2rad(alfa_from0toSinv(i)), XAC, XCG, bCG, MAC);
            CMCT_from0toSinv(i) = CT_contr(obj1, CD_from0toSinv(i), Thrust_axes, MAC);
            CMCG_from0toSinv(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_from0toSinv(i));
            % HORIZONTAL TAIL LIFT COEFFICIENT
            CLHT_from0toSinv(i) = CL_Tail(obj1, CMCL_from0toSinv(i), CMCD_from0toSinv(i), ...
                                             CMCT_from0toSinv(i), CMCG_from0toSinv(i), ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_from0toSinv(i)));
            % HORIZONTAL TAIL LIFT
            LHT_from0toSinv(i) = (0.5)*(V_from0toSinv(i)^2)*(S)*(rho0)*(CLHT_from0toSinv(i))*(1e-1);
            
            % DRAFT VERSION OF ITERATION 
            CL_tail            = CLHT_from0toSinv(i);
            CL_wb              = CL_from0toSinv(i);
            CL_new_from0toSinv = CL_wb - CL_tail;
            tol                = 1e-3; 
            n                  = 1;
            while abs(CL_wb - CL_tail) > tol
               alfa_new_from0toSinv = alpha_fullmodel_inv(CL_new_from0toSinv, CL_max_inv, CL0, CLalfa, alfa_0l);
               CD_wb              = polyval(p_CD_wb1, alfa_new_from0toSinv);
               CMCL_new           = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_from0toSinv), XAC, XCG, bCG, MAC);
               CMCD_new           = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_from0toSinv), XAC, XCG, bCG, MAC);
               CMCT_new           = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
               CMCG_new           = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
               CLHT_new           = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_new_from0toSinv));
               CL_tail            = CLHT_new;
               CL_new_from0toSinv = CL_from0toSinv(i) - CL_tail;
               CL_wb              = CL_wb + (CL_new_from0toSinv - CL_wb) * 1e-1;
               n                  = n + 1;
               if n == 15
                   break
               end
            end
            CL_from0toSinv_fullvehicle(i)   = CL_new_from0toSinv; 
            CLHT_from0toSinv(i) = CL_tail;
            alfa_from0toSinv_fullvehicle(i) = alfa_new_from0toSinv;
%             alfa_from0toSinv(i) = alfa_new_from0toSinv;
            WBL_from0toSinv(i)  = q_from0toSinv(i) * S * CL_from0toSinv(i) * 1e-1;
            WBL_from0toSinv_new = q_from0toSinv(i) * S * CL_from0toSinv_fullvehicle(i) * 1e-1;
            LHT_from0toSinv(i)  = (0.5)*(V_from0toSinv(i)^2)*(S)*(rho0)*(CLHT_from0toSinv(i))*(1e-1);
            CMCL_from0toSinv(i) = CMCL_new;
            CMCD_from0toSinv(i) = CMCD_new;
            CMCG_from0toSinv(i) = CMCG_new;             
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv.value = CL_from0toSinv;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.value = CL_from0toSinv_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv.value = alfa_from0toSinv;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_rad.value = deg2rad(alfa_from0toSinv);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_rad.Attributes.unit = "rad";     
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_new.value = alfa_from0toSinv_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_new_rad.value = deg2rad(alfa_from0toSinv_fullvehicle);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_new_rad.Attributes.unit = "rad";         
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toSinv.value = CD_from0toSinv;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toSinv.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toSinv.value = q_from0toSinv;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toSinv.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toSinv.value = WBL_from0toSinv;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toSinv.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.value = WBL_from0toSinv_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_from0toSinv.value = CMCL_from0toSinv;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_from0toSinv.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_from0toSinv.value = CMCD_from0toSinv;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_from0toSinv.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_from0toSinv.value = CMCT_from0toSinv;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_from0toSinv.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_from0toSinv.value = CMCG_from0toSinv;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_from0toSinv.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_from0toSinv.value = CLHT_from0toSinv;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_from0toSinv.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_from0toSinv.value = LHT_from0toSinv;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_from0toSinv.Attributes.unit = "daN";   
        
        % =================================================================  
        % POINT G 
        CL_G   = CLmax_func(rho0, VG, WS, nG);
            if CL_G < CL_max_inv
                CL_G = CL_G;
            elseif CL_G > CL_max_inv
                CL_G = CL_max_inv;
            end
        alfa_G = alpha_fullmodel_inv(CL_G, CL_max_inv, CL0, CLalfa, alfa_0l);
        CD_G   = polyval(p_CD_wb1, alfa_G);
        q_G    = 0.5*rho0*(VG)^2;
        WBL_G  = q_G * S * CL_G * 1e-1; 
        CMCL_G = CLWB_contrib(obj1, CL_G, deg2rad(alfa_G), XAC, XCG, bCG, MAC);
        CMCD_G = CDWB_contrib(obj1, CL_G, deg2rad(alfa_G), XAC, XCG, bCG, MAC);
        CMCT_G = CT_contr(obj1, CD_G, Thrust_axes, MAC);
        CMCG_G = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_G);
        % HORIZONTAL TAIL LIFT COEFFICIENT
        CLHT_G = CL_Tail(obj1, CMCL_G, CMCD_G, CMCT_G, CMCG_G, l_ht, MAC, XAC, XCG, deg2rad(alfa_G)); 
        % HORIZONTAL TAIL LIFT
        LHT_G = (0.5)*(VG^2)*(S)*(rho0)*(CLHT_G)*(1e-1);
        
        % DRAFT VERSION OF ITERATION 
        CL_tail  = CLHT_G;
        CL_wb    = CL_G;
        CL_new_G = CL_wb - CL_tail;
        tol      = 1e-3; 
        n        = 1;
        while abs(CL_wb - CL_tail) > tol
           alfa_new_G = alpha_fullmodel_inv(CL_new_G, CL_max_inv, CL0, CLalfa, alfa_0l);
           CD_wb     = polyval(p_CD_wb1, alfa_new_G);
           CMCL_new  = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_G), XAC, XCG, bCG, MAC);
           CMCD_new  = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_G), XAC, XCG, bCG, MAC);
           CMCT_new  = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
           CMCG_new  = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
           CLHT_new  = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                               l_ht, MAC, XAC, XCG, deg2rad(alfa_new_G));
           CL_tail   = CLHT_new;
           CL_new_G  = CL_G - CL_tail;
           CL_wb     = CL_wb + (CL_new_G - CL_wb) * 1e-1;
           n         = n + 1;
           if n == 15
               break
           end
        end
        CL_G_fullvehicle = CL_new_G; 
        CLHT_G           = CL_tail;
%         alfa_G = alfa_new_G;
        WBL_G     = q_G * S * CL_G * 1e-1;
        WBL_G_new = q_G * S * CL_G_fullvehicle * 1e-1;
        LHT_G     = (0.5)*(VG^2)*(S)*(rho0)*(CLHT_G)*(1e-1);
        CMCL_G    = CMCL_new;
        CMCD_G    = CMCD_new;
        CMCG_G    = CMCG_new;          
        
        
        % STORE ALL THE DATA INSIDE THE STRUCT VARIABLE
        % POINT G
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value = q_G;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.Attributes.unit = "Pa";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG.value = WBL_G;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG_new.value = WBL_G_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.value = LHT_G;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.value = alfa_G;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_rad.value = deg2rad(alfa_G);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfa_new_G.value = alfa_new_G;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfa_new_G.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfa_new_G_rad.value = deg2rad(alfa_new_G);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfa_new_G_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value = CL_G;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G_new.value = CL_G_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CLHT_G.value = CLHT_G;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CLHT_G.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.value = CD_G;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CM_G.value = CMCG_G;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CM_G.Attributes.unit = "Non dimensional";          
        % =================================================================  
        % FROM Sinv TO G1  
        CL_fromSinvtoG1      = zeros(length(V_fromSinvtoG1), 1);
        alfa_fromSinvtoG1    = zeros(length(V_fromSinvtoG1), 1);
        CD_fromSinvtoG1      = zeros(length(V_fromSinvtoG1), 1);
        q_fromSinvtoG1       = zeros(length(V_fromSinvtoG1), 1);
        WBL_fromSinvtoG1     = zeros(length(V_fromSinvtoG1), 1);
        WBL_fromSinvtoG1_new = zeros(length(V_fromSinvtoG1), 1);
        CMCL_fromSinvtoG1    = zeros(length(V_fromSinvtoG1), 1);
        CMCD_fromSinvtoG1    = zeros(length(V_fromSinvtoG1), 1);
        CMCT_fromSinvtoG1    = zeros(length(V_fromSinvtoG1), 1);
        CMCG_fromSinvtoG1    = zeros(length(V_fromSinvtoG1), 1); 
        CLHT_fromSinvtoG1    = zeros(length(V_fromSinvtoG1), 1);
        LHT_fromSinvtoG1     = zeros(length(V_fromSinvtoG1), 1);
        % =================================================================
        CL_fromSinvtoG1_fullvehicle   = zeros(length(V_fromSinvtoG1), 1);
        alfa_fromSinvtoG1_fullvehicle = zeros(length(V_fromSinvtoG1), 1);
        % ================================================================= 
        for i = 1:length(V_fromSinvtoG1)
            CL_fromSinvtoG1(i)   = CLmax_func(rho0, V_fromSinvtoG1(i), WS, n_fromSinvtoG1(i));
            if CL_fromSinvtoG1(i) < CL_max_inv
                CL_fromSinvtoG1(i) = CL_fromSinvtoG1(i);
            elseif CL_fromSinvtoG1(i) > CL_max_inv
                CL_fromSinvtoG1(i) = CL_max_inv;
            end
            alfa_fromSinvtoG1(i) = alpha_fullmodel_inv(CL_fromSinvtoG1(i), CL_max_inv, CL0, CLalfa, alfa_0l);
            CD_fromSinvtoG1(i)   = polyval(p_CD_wb1, alfa_fromSinvtoG1(i));
            q_fromSinvtoG1(i)    = 0.5*rho0*(V_fromSinvtoG1(i))^2;
            WBL_fromSinvtoG1(i)  = q_fromSinvtoG1(i)*S*CL_fromSinvtoG1(i)*1e-1;
            CMCL_fromSinvtoG1(i) = CLWB_contrib(obj1, CL_fromSinvtoG1(i), deg2rad(alfa_fromSinvtoG1(i)), XAC, XCG, bCG, MAC);
            CMCD_fromSinvtoG1(i) = CDWB_contrib(obj1, CL_fromSinvtoG1(i), deg2rad(alfa_fromSinvtoG1(i)), XAC, XCG, bCG, MAC);
            CMCT_fromSinvtoG1(i) = CT_contr(obj1, CD_fromSinvtoG1(i), Thrust_axes, MAC);
            CMCG_fromSinvtoG1(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromSinvtoG1(i));
            % HORIZONTAL TAIL LIFT COEFFICIENT
            CLHT_fromSinvtoG1(i) = CL_Tail(obj1, CMCL_fromSinvtoG1(i), CMCD_fromSinvtoG1(i), ...
                                             CMCT_fromSinvtoG1(i), CMCG_fromSinvtoG1(i), ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_fromSinvtoG1(i)));
            % HORIZONTAL TAIL LIFT
            LHT_fromSinvtoG1(i) = (0.5)*(V_fromSinvtoG1(i)^2)*(S)*(rho0)*(CLHT_fromSinvtoG1(i))*(1e-1);
            
            % DRAFT VERSION OF ITERATION 
            CL_tail            = CLHT_fromSinvtoG1(i);
            CL_wb              = CL_fromSinvtoG1(i);
            CL_new_fromSinvtoG1 = CL_wb - CL_tail;
            tol                = 1e-3; 
            n                  = 1;
            while abs(CL_wb - CL_tail) > tol
               alfa_new_fromSinvtoG1 = alpha_fullmodel_inv(CL_new_fromSinvtoG1, CL_max_inv, CL0, CLalfa, alfa_0l);
               CD_wb              = polyval(p_CD_wb1, alfa_new_fromSinvtoG1);
               CMCL_new           = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromSinvtoG1), XAC, XCG, bCG, MAC);
               CMCD_new           = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromSinvtoG1), XAC, XCG, bCG, MAC);
               CMCT_new           = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
               CMCG_new           = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
               CLHT_new           = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromSinvtoG1));
               CL_tail             = CLHT_new;
               CL_new_fromSinvtoG1 = CL_fromSinvtoG1(i) - CL_tail;
               CL_wb               = CL_wb + (CL_new_fromSinvtoG1 - CL_wb) * 1e-1;
               n                   = n + 1;
               if n == 15
                   break
               end
            end
            CL_fromSinvtoG1_fullvehicle(i)   = CL_new_fromSinvtoG1; 
            CLHT_fromSinvtoG1(i) = CL_tail;
            alfa_fromSinvtoG1_fullvehicle(i) = alfa_new_fromSinvtoG1;
%             alfa_fromSinvtoG1(i) = alfa_new_fromSinvtoG1;
            WBL_fromSinvtoG1(i)  = q_fromSinvtoG1(i) * S * CL_fromSinvtoG1(i) * 1e-1;
            WBL_fromSinvtoG1_new = q_fromSinvtoG1(i) * S * CL_fromSinvtoG1_fullvehicle(i) * 1e-1;
            LHT_fromSinvtoG1(i)  = (0.5)*(V_fromSinvtoG1(i)^2)*(S)*(rho0)*(CLHT_fromSinvtoG1(i))*(1e-1);
            CMCL_fromSinvtoG1(i) = CMCL_new;
            CMCD_fromSinvtoG1(i) = CMCD_new;
            CMCG_fromSinvtoG1(i) = CMCG_new;              
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG1.value = CL_fromSinvtoG1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG1.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG1_new.value = CL_fromSinvtoG1_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG1_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG1.value = alfa_fromSinvtoG1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG1.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG1_rad.value = deg2rad(alfa_fromSinvtoG1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG1_rad.Attributes.unit = "rad";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG1_new.value = alfa_fromSinvtoG1_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG1_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG1_new_rad.value = deg2rad(alfa_fromSinvtoG1_fullvehicle);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG1_new_rad.Attributes.unit = "rad";       
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromSinvtoG1.value = CD_fromSinvtoG1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromSinvtoG1.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromSinvtoG1.value = q_fromSinvtoG1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromSinvtoG1.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromSinvtoG1.value = WBL_fromSinvtoG1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromSinvtoG1.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromSinvtoG1.value = CMCL_fromSinvtoG1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromSinvtoG1.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromSinvtoG1.value = CMCD_fromSinvtoG1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromSinvtoG1.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromSinvtoG1.value = CMCT_fromSinvtoG1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromSinvtoG1.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromSinvtoG1.value = CMCG_fromSinvtoG1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromSinvtoG1.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromSinvtoG1.value = CLHT_fromSinvtoG1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromSinvtoG1.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromSinvtoG1.value = LHT_fromSinvtoG1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromSinvtoG1.Attributes.unit = "daN";  
        % =================================================================
        % FROM G1 TO F 
        CL_fromG1toF      = zeros(length(V_fromG1toF), 1);
        alfa_fromG1toF    = zeros(length(V_fromG1toF), 1);
        CD_fromG1toF      = zeros(length(V_fromG1toF), 1);
        q_fromG1toF       = zeros(length(V_fromG1toF), 1);
        WBL_fromG1toF     = zeros(length(V_fromG1toF), 1);
        WBL_fromG1toF_new = zeros(length(V_fromG1toF), 1);
        CMCL_fromG1toF    = zeros(length(V_fromG1toF), 1);
        CMCD_fromG1toF    = zeros(length(V_fromG1toF), 1);
        CMCT_fromG1toF    = zeros(length(V_fromG1toF), 1);
        CMCG_fromG1toF    = zeros(length(V_fromG1toF), 1); 
        CLHT_fromG1toF    = zeros(length(V_fromG1toF), 1);
        LHT_fromG1toF     = zeros(length(V_fromG1toF), 1);
        % =================================================================
        CL_fromG1toF_fullvehicle   = zeros(length(V_fromG1toF), 1);
        alfa_fromG1toF_fullvehicle = zeros(length(V_fromG1toF), 1);
        % ================================================================= 
        for i = 1:length(V_fromG1toF)
            CL_fromG1toF(i)   = CLmax_func(rho0, V_fromG1toF(i), WS, n_fromG1toF(i));
            if CL_fromG1toF(i) < CL_max_inv
                CL_fromG1toF(i) = CL_fromG1toF(i);
            elseif CL_fromG1toF(i) > CL_max_inv
                CL_fromG1toF(i) = CL_max_inv;
            end
            alfa_fromG1toF(i) = alpha_fullmodel_inv(CL_fromG1toF(i), CL_max_inv, CL0, CLalfa, alfa_0l);
            CD_fromG1toF(i)   = polyval(p_CD_wb1, alfa_fromG1toF(i));
            q_fromG1toF(i)    = 0.5*rho0*(V_fromG1toF(i))^2;
            WBL_fromG1toF(i)  = q_fromG1toF(i)*S*CL_fromG1toF(i)*1e-1;
            CMCL_fromG1toF(i) = CLWB_contrib(obj1, CL_fromG1toF(i), deg2rad(alfa_fromG1toF(i)), XAC, XCG, bCG, MAC);
            CMCD_fromG1toF(i) = CDWB_contrib(obj1, CL_fromG1toF(i), deg2rad(alfa_fromG1toF(i)), XAC, XCG, bCG, MAC);
            CMCT_fromG1toF(i) = CT_contr(obj1, CD_fromG1toF(i), Thrust_axes, MAC);
            CMCG_fromG1toF(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromG1toF(i));
            % HORIZONTAL TAIL LIFT COEFFICIENT
            CLHT_fromG1toF(i) = CL_Tail(obj1, CMCL_fromG1toF(i), CMCD_fromG1toF(i), ...
                                             CMCT_fromG1toF(i), CMCG_fromG1toF(i), ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_fromG1toF(i)));
            % HORIZONTAL TAIL LIFT
            LHT_fromG1toF(i) = (0.5)*(V_fromG1toF(i)^2)*(S)*(rho0)*(CLHT_fromG1toF(i))*(1e-1);
            
            % DRAFT VERSION OF ITERATION 
            CL_tail            = CLHT_fromG1toF(i);
            CL_wb              = CL_fromG1toF(i);
            CL_new_fromG1toF = CL_wb - CL_tail;
            tol                = 1e-3; 
            n                  = 1;
            while abs(CL_wb - CL_tail) > tol
               alfa_new_fromG1toF = alpha_fullmodel_inv(CL_new_fromG1toF, CL_max_inv, CL0, CLalfa, alfa_0l);
               CD_wb              = polyval(p_CD_wb1, alfa_new_fromG1toF);
               CMCL_new           = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromG1toF), XAC, XCG, bCG, MAC);
               CMCD_new           = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromG1toF), XAC, XCG, bCG, MAC);
               CMCT_new           = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
               CMCG_new           = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
               CLHT_new           = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromG1toF));
               CL_tail            = CLHT_new;
               CL_new_fromG1toF = CL_fromG1toF(i) - CL_tail;
               CL_wb              = CL_wb + (CL_new_fromG1toF - CL_wb) * 1e-1;
               n                  = n + 1;
               if n == 15
                   break
               end
            end
            CL_fromG1toF_fullvehicle(i)   = CL_new_fromG1toF; 
            CLHT_fromG1toF(i) = CL_tail;
            alfa_fromG1toF_fullvehicle(i) = alfa_new_fromG1toF;
%             alfa_fromG1toF(i) = alfa_new_fromG1toF;
            WBL_fromG1toF(i)     = q_fromG1toF(i) * S * CL_fromG1toF(i) * 1e-1;
            WBL_fromG1toF_new(i) = q_fromG1toF(i) * S * CL_fromG1toF_fullvehicle(i) * 1e-1;
            LHT_fromG1toF(i)     = (0.5)*(V_fromG1toF(i)^2)*(S)*(rho0)*(CLHT_fromG1toF(i))*(1e-1);
            CMCL_fromG1toF(i)    = CMCL_new;
            CMCD_fromG1toF(i)    = CMCD_new;
            CMCG_fromG1toF(i)    = CMCG_new;              
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF.value = CL_fromG1toF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF_new.value = CL_fromG1toF_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF.value = alfa_fromG1toF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF_rad.value = deg2rad(alfa_fromG1toF);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF_rad.Attributes.unit = "rad";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF_new.value = alfa_fromG1toF_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF_new_rad.value = deg2rad(alfa_fromG1toF_fullvehicle);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF_new_rad.Attributes.unit = "rad";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromG1toF.value = CD_fromG1toF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromG1toF.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromG1toF.value = q_fromG1toF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromG1toF.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromG1toF.value = WBL_fromG1toF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromG1toF.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG1toF_new.value = WBL_fromG1toF_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG1toF_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromG1toF.value = CMCL_fromG1toF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromG1toF.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromG1toF.value = CMCD_fromG1toF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromG1toF.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromG1toF.value = CMCT_fromG1toF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromG1toF.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromG1toF.value = CMCG_fromG1toF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromG1toF.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromG1toF.value = CLHT_fromG1toF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromG1toF.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromG1toF.value = LHT_fromG1toF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromG1toF.Attributes.unit = "daN";         
        % ================================================================= 
        % FROM F TO E
        CL_fromFtoE      = zeros(length(V_fromFtoE), 1);
        alfa_fromFtoE    = zeros(length(V_fromFtoE), 1);
        CD_fromFtoE      = zeros(length(V_fromFtoE), 1);
        q_fromFtoE       = zeros(length(V_fromFtoE), 1);
        WBL_fromFtoE     = zeros(length(V_fromFtoE), 1);
        WBL_fromFtoE_new = zeros(length(V_fromFtoE), 1);
        CMCL_fromFtoE    = zeros(length(V_fromFtoE), 1);
        CMCD_fromFtoE    = zeros(length(V_fromFtoE), 1);
        CMCT_fromFtoE    = zeros(length(V_fromFtoE), 1);
        CMCG_fromFtoE    = zeros(length(V_fromFtoE), 1); 
        CLHT_fromFtoE    = zeros(length(V_fromFtoE), 1);
        LHT_fromFtoE     = zeros(length(V_fromFtoE), 1);
        % =================================================================
        CL_fromFtoE_fullvehicle   = zeros(length(V_fromFtoE), 1);
        alfa_fromFtoE_fullvehicle = zeros(length(V_fromFtoE), 1);
        % ================================================================= 
        for i = 1:length(V_fromFtoE)
            CL_fromFtoE(i)   = CLmax_func(rho0, V_fromFtoE(i), WS, n_fromFtoE(i));
            if CL_fromFtoE(i) < CL_max_inv
                CL_fromFtoE(i) = CL_fromFtoE(i);
            elseif CL_fromFtoE(i) > CL_max_inv
                CL_fromFtoE(i) = CL_max_inv;
            end
            alfa_fromFtoE(i) = alpha_fullmodel_inv(CL_fromFtoE(i), CL_max_inv, CL0, CLalfa, alfa_0l);
            CD_fromFtoE(i)   = polyval(p_CD_wb1, alfa_fromFtoE(i));
            q_fromFtoE(i)    = 0.5*rho0*(V_fromFtoE(i))^2;
            WBL_fromFtoE(i)  = q_fromFtoE(i)*S*CL_fromFtoE(i)*1e-1;
            CMCL_fromFtoE(i) = CLWB_contrib(obj1, CL_fromFtoE(i), deg2rad(alfa_fromFtoE(i)), XAC, XCG, bCG, MAC);
            CMCD_fromFtoE(i) = CDWB_contrib(obj1, CL_fromFtoE(i), deg2rad(alfa_fromFtoE(i)), XAC, XCG, bCG, MAC);
            CMCT_fromFtoE(i) = CT_contr(obj1, CD_fromFtoE(i), Thrust_axes, MAC);
            CMCG_fromFtoE(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromFtoE(i));
            % HORIZONTAL TAIL LIFT COEFFICIENT
            CLHT_fromFtoE(i) = CL_Tail(obj1, CMCL_fromFtoE(i), CMCD_fromFtoE(i), ...
                                             CMCT_fromFtoE(i), CMCG_fromFtoE(i), ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_fromFtoE(i)));
            % HORIZONTAL TAIL LIFT
            LHT_fromFtoE(i) = (0.5)*(V_fromFtoE(i)^2)*(S)*(rho0)*(CLHT_fromFtoE(i))*(1e-1);
            
            % DRAFT VERSION OF ITERATION 
            CL_tail         = CLHT_fromFtoE(i);
            CL_wb           = CL_fromFtoE(i);
            CL_new_fromFtoE = CL_wb - CL_tail;
            tol             = 1e-3; 
            n               = 1;
            while abs(CL_wb - CL_tail) > tol
               alfa_new_fromFtoE = alpha_fullmodel_inv(CL_new_fromFtoE, CL_max_inv, CL0, CLalfa, alfa_0l);
               CD_wb             = polyval(p_CD_wb1, alfa_new_fromFtoE);
               CMCL_new          = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromFtoE), XAC, XCG, bCG, MAC);
               CMCD_new          = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromFtoE), XAC, XCG, bCG, MAC);
               CMCT_new          = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
               CMCG_new          = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
               CLHT_new          = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromFtoE));
               CL_tail           = CLHT_new;
               CL_new_fromFtoE   = CL_fromFtoE(i) - CL_tail;
               CL_wb             = CL_wb + (CL_new_fromFtoE - CL_wb) * 1e-1;
               n                 = n + 1;
               if n == 15
                   break
               end
            end
            CL_fromFtoE_fullvehicle(i)   = CL_new_fromFtoE; 
            CLHT_fromFtoE(i) = CL_tail;
            alfa_fromFtoE_fullvehicle(i) = alfa_new_fromFtoE;
%             alfa_fromFtoE(i) = alfa_new_fromFtoE;
            WBL_fromFtoE(i)  = q_fromFtoE(i) * S * CL_fromFtoE(i) * 1e-1;
            WBL_fromFtoE_new = q_fromFtoE(i) * S * CL_fromFtoE_fullvehicle(i) * 1e-1;
            LHT_fromFtoE(i)  = (0.5)*(V_fromFtoE(i)^2)*(S)*(rho0)*(CLHT_fromFtoE(i))*(1e-1);
            CMCL_fromFtoE(i) = CMCL_new;
            CMCD_fromFtoE(i) = CMCD_new;
            CMCG_fromFtoE(i) = CMCG_new;              
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value = CL_fromFtoE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value = CL_fromFtoE_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE.value = alfa_fromFtoE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE_rad.value = deg2rad(alfa_fromFtoE);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE_rad.Attributes.unit = "rad";     
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE_new.value = alfa_fromFtoE_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE_new_rad.value = deg2rad(alfa_fromFtoE_fullvehicle);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE_new_rad.Attributes.unit = "rad";      
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value = CD_fromFtoE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromFtoE.value = q_fromFtoE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromFtoE.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE.value = WBL_fromFtoE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoE_new.value = WBL_fromFtoE_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoE_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromFtoE.value = CMCL_fromFtoE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromFtoE.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromFtoE.value = CMCD_fromFtoE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromFtoE.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromFtoE.value = CMCT_fromFtoE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromFtoE.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromFtoE.value = CMCG_fromFtoE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromFtoE.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromFtoE.value = CLHT_fromFtoE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromFtoE.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromFtoE.value = LHT_fromFtoE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromFtoE.Attributes.unit = "daN";        
        % =================================================================
        % FROM E TO 0
        CL_fromEto0      = zeros(length(V_fromEto0), 1);
        alfa_fromEto0    = zeros(length(V_fromEto0), 1);
        CD_fromEto0      = zeros(length(V_fromEto0), 1);
        q_fromEto0       = zeros(length(V_fromEto0), 1);
        WBL_fromEto0     = zeros(length(V_fromEto0), 1);
        WBL_fromEto0_new = zeros(length(V_fromEto0), 1);
        CMCL_fromEto0    = zeros(length(V_fromEto0), 1);
        CMCD_fromEto0    = zeros(length(V_fromEto0), 1);
        CMCT_fromEto0    = zeros(length(V_fromEto0), 1);
        CMCG_fromEto0    = zeros(length(V_fromEto0), 1); 
        CLHT_fromEto0    = zeros(length(V_fromEto0), 1);
        LHT_fromEto0     = zeros(length(V_fromEto0), 1);
        % =================================================================
        CL_fromEto0_fullvehicle   = zeros(length(V_fromEto0), 1);
        alfa_fromEto0_fullvehicle = zeros(length(V_fromEto0), 1);
        % ================================================================= 
        for i = 1:length(V_fromEto0)
            CL_fromEto0(i)   = CLmax_func(rho0, V_fromEto0(i), WS, n_fromEto0(i));
            if CL_fromEto0(i) < CL_max_inv
                CL_fromEto0(i) = CL_fromEto0(i);
            elseif CL_fromEto0(i) > CL_max_inv
                CL_fromEto0(i) = CL_max_inv;
            end
            alfa_fromEto0(i) = alpha_fullmodel_inv(CL_fromEto0(i), CL_max_inv, CL0, CLalfa, alfa_0l);
            CD_fromEto0(i)   = polyval(p_CD_wb1, alfa_fromEto0(i));
            q_fromEto0(i)    = 0.5*rho0*(V_fromEto0(i))^2;
            WBL_fromEto0(i)  = q_fromEto0(i)*S*CL_fromEto0(i)*1e-1;
            CMCL_fromEto0(i) = CLWB_contrib(obj1, CL_fromEto0(i), deg2rad(alfa_fromEto0(i)), XAC, XCG, bCG, MAC);
            CMCD_fromEto0(i) = CDWB_contrib(obj1, CL_fromEto0(i), deg2rad(alfa_fromEto0(i)), XAC, XCG, bCG, MAC);
            CMCT_fromEto0(i) = CT_contr(obj1, CD_fromEto0(i), Thrust_axes, MAC);
            CMCG_fromEto0(i) = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_fromEto0(i));
            % HORIZONTAL TAIL LIFT COEFFICIENT
            CLHT_fromEto0(i) = CL_Tail(obj1, CMCL_fromEto0(i), CMCD_fromEto0(i), ...
                                             CMCT_fromEto0(i), CMCG_fromEto0(i), ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_fromEto0(i)));
            % HORIZONTAL TAIL LIFT
            LHT_fromEto0(i) = (0.5)*(V_fromEto0(i)^2)*(S)*(rho0)*(CLHT_fromEto0(i))*(1e-1);
            
            % DRAFT VERSION OF ITERATION 
            CL_tail         = CLHT_fromEto0(i);
            CL_wb           = CL_fromEto0(i);
            CL_new_fromEto0 = CL_wb - CL_tail;
            tol             = 1e-3; 
            n               = 1;
            while abs(CL_wb - CL_tail) > tol
               alfa_new_fromEto0 = alpha_fullmodel_inv(CL_new_fromEto0, CL_max_inv, CL0, CLalfa, alfa_0l);
               CD_wb             = polyval(p_CD_wb1, alfa_new_fromEto0);
               CMCL_new          = CLWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromEto0), XAC, XCG, bCG, MAC);
               CMCD_new          = CDWB_contrib(obj1, CL_wb, deg2rad(alfa_new_fromEto0), XAC, XCG, bCG, MAC);
               CMCT_new          = CT_contr(obj1, CD_wb, Thrust_axes, MAC);
               CMCG_new          = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_wb);
               CLHT_new          = CL_Tail(obj1, CMCL_new, CMCD_new, CMCT_new, CMCG_new, ...
                                             l_ht, MAC, XAC, XCG, deg2rad(alfa_new_fromEto0));
               CL_tail           = CLHT_new;
               CL_new_fromEto0   = CL_fromEto0(i) - CL_tail;
               CL_wb             = CL_wb + (CL_new_fromEto0 - CL_wb) * 1e-1;
               n                 = n + 1;
               if n == 15
                   break
               end
            end
            CL_fromEto0_fullvehicle(i)   = CL_new_fromEto0; 
            CLHT_fromEto0(i) = CL_tail;
            alfa_fromEto0_fullvehicle(i) = alfa_new_fromEto0;
%             alfa_fromEto0(i) = alfa_new_fromEto0;
            WBL_fromEto0(i)     = q_fromEto0(i) * S * CL_fromEto0(i) * 1e-1;
            WBL_fromEto0_new(i) = q_fromEto0(i) * S * CL_fromEto0_fullvehicle(i) * 1e-1;
            LHT_fromEto0(i)     = (0.5)*(V_fromEto0(i)^2)*(S)*(rho0)*(CLHT_fromEto0(i))*(1e-1);
            CMCL_fromEto0(i)    = CMCL_new;
            CMCD_fromEto0(i)    = CMCD_new;
            CMCG_fromEto0(i)    = CMCG_new;              
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0.value = CL_fromEto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.value = CL_fromEto0_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0.value = alfa_fromEto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_rad.value = deg2rad(alfa_fromEto0);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_rad.Attributes.unit = "rad";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_new.value = alfa_fromEto0_fullvehicle;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_new_rad.value = deg2rad(alfa_fromEto0_fullvehicle);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_new_rad.Attributes.unit = "rad";       
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromEto0.value = CD_fromEto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromEto0.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromEto0.value = q_fromEto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromEto0.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromEto0.value = WBL_fromEto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromEto0.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.value = WBL_fromEto0_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromEto0.value = CMCL_fromEto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCL_fromEto0.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromEto0.value = CMCD_fromEto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCD_fromEto0.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromEto0.value = CMCT_fromEto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCT_fromEto0.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromEto0.value = CMCG_fromEto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMCG_fromEto0.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromEto0.value = CLHT_fromEto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_fromEto0.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromEto0.value = LHT_fromEto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_fromEto0.Attributes.unit = "daN"; 
        % =================================================================
        % TAIL LOADS DIAGRAM - NEGATIVE SIDE 
        figure(5); 
        % hold on; grid on; grid minor; 
        % ---------------------------------------------------------------------
        plot(V_from0toSinv,  LHT_from0toSinv,  '-r', 'LineWidth', 1)
        plot(V_fromSinvtoG1, LHT_fromSinvtoG1, '-r', 'LineWidth', 1)
        plot(V_fromG1toF,    LHT_fromG1toF,    '-r', 'LineWidth', 1)
        plot(V_fromFtoE,     LHT_fromFtoE,     '-r', 'LineWidth', 1)
        plot(V_fromEto0,     LHT_fromEto0,     '-r', 'LineWidth', 1)
        % ---------------------------------------------------------------------
        plot(V_from0toSinv(1),    LHT_from0toSinv(1),    'k.', 'MarkerSize', 10)
        plot(V_from0toSinv(end),  LHT_from0toSinv(end),  'k.', 'MarkerSize', 10)
        plot(VG,                  LHT_G,                 'k.', 'MarkerSize', 10)
        plot(V_fromSinvtoG1(end), LHT_fromSinvtoG1(end), 'k.', 'MarkerSize', 10)
        plot(V_fromG1toF(end),    LHT_fromG1toF(end),    'k.', 'MarkerSize', 10)
        plot(V_fromFtoE(end),     LHT_fromFtoE(end),     'k.', 'MarkerSize', 10)
        plot(V_fromEto0(end),     LHT_fromEto0(end),     'k.', 'MarkerSize', 10)
        % ---------------------------------------------------------------------
        text(V_from0toSinv(end),  LHT_from0toSinv(end),  ' S inv.', 'FontSize', 6)
        text(VG,                  LHT_G,                 ' G',      'FontSize', 6)
        text(V_fromSinvtoG1(end), LHT_fromSinvtoG1(end), ' G1',     'FontSize', 6)
        text(V_fromG1toF(end),    LHT_fromG1toF(end),    ' F',      'FontSize', 6)
        text(V_fromFtoE(end),     LHT_fromFtoE(end),     ' E',      'FontSize', 6)
        ylim padded;
        xlim padded;
%         % MAIN WING LOADS DIAGRAM        
%         CL_from0toSinv_new  = CL_from0toSinv  + CLHT_from0toSinv;        
%         CL_fromSinvtoG1_new = CL_fromSinvtoG1 + CLHT_fromSinvtoG1;        
%         CL_fromG1toF_new    = CL_fromG1toF + CLHT_fromG1toF;        
%         CL_fromFtoE_new     = CL_fromFtoE + CLHT_fromFtoE;        
%         CL_fromEto0_new     = CL_fromEto0 + CLHT_fromEto0;

%         % MAIN WING LOADS DIAGRAM        
%         CL_from0toSinv_new  = CL_from0toSinv  - CLHT_from0toSinv;        
%         CL_fromSinvtoG1_new = CL_fromSinvtoG1 - CLHT_fromSinvtoG1;        
%         CL_fromG1toF_new    = CL_fromG1toF - CLHT_fromG1toF;        
%         CL_fromFtoE_new     = CL_fromFtoE - CLHT_fromFtoE;        
%         CL_fromEto0_new     = CL_fromEto0 - CLHT_fromEto0;
%         LW_from0toSinv_new  = zeros(length(CL_from0toSinv_new), 1);
%         LW_fromSinvtoG1_new = zeros(length(CL_from0toSinv_new), 1);
%         LW_fromG1toF_new    = zeros(length(CL_from0toSinv_new), 1);
%         LW_fromFtoE_new     = zeros(length(CL_from0toSinv_new), 1);
%         LW_fromEto0_new     = zeros(length(CL_from0toSinv_new), 1);
%         for i = 1:length(CL_from0toSinv_new)
%             LW_from0toSinv_new(i)  = (0.5)*(V_from0toSinv(i)^2)* rho0 * S *(CL_from0toSinv_new(i))*(1e-1);  
%             LW_fromSinvtoG1_new(i) = (0.5)*(V_fromSinvtoG1(i)^2)* rho0 * S *(CL_fromSinvtoG1_new(i))*(1e-1);
%             LW_fromG1toF_new(i)    = (0.5)*(V_fromG1toF(i)^2)* rho0 * S *(CL_fromG1toF_new(i))*(1e-1);
%             LW_fromFtoE_new(i)     = (0.5)*(V_fromFtoE(i)^2)* rho0 * S *(CL_fromFtoE_new(i))*(1e-1);
%             LW_fromEto0_new(i)     = (0.5)*(V_fromEto0(i)^2)* rho0 * S *(CL_fromEto0_new(i))*(1e-1);
% %             LW_from0toSinv_new(i)  = (0.5)*(V_from0toSinv(i)^2)* rho0 * S *(CL_from0toSinv(i))*(1e-1);  
% %             LW_fromSinvtoG1_new(i) = (0.5)*(V_fromSinvtoG1(i)^2)* rho0 * S *(CL_fromSinvtoG1(i))*(1e-1);
% %             LW_fromG1toF_new(i)    = (0.5)*(V_fromG1toF(i)^2)* rho0 * S *(CL_fromG1toF(i))*(1e-1);
% %             LW_fromFtoE_new(i)     = (0.5)*(V_fromFtoE(i)^2)* rho0 * S *(CL_fromFtoE(i))*(1e-1);
% %             LW_fromEto0_new(i)     = (0.5)*(V_fromEto0(i)^2)* rho0 * S *(CL_fromEto0(i))*(1e-1); 
%         end       
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.value = CL_from0toSinv_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.Attributes.unit = "Non dimensional";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG1_new.value = CL_fromSinvtoG1_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG1_new.Attributes.unit = "Non dimensional";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF_new.value = CL_fromG1toF_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF_new.Attributes.unit = "Non dimensional";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value = CL_fromFtoE_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.Attributes.unit = "Non dimensional";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.value = CL_fromEto0_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.Attributes.unit = "Non dimensional";
%         % ------------------------------------------------------------------------------------------------------------------------
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.value = LW_from0toSinv_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.Attributes.unit = "daN";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromSinvtoG1_new.value = LW_fromSinvtoG1_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromSinvtoG1_new.Attributes.unit = "daN";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG1toF_new.value = LW_fromG1toF_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG1toF_new.Attributes.unit = "daN";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoE_new.value = LW_fromFtoE_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoE_new.Attributes.unit = "daN";
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.value = LW_fromEto0_new;
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.Attributes.unit = "daN";
        
        % WING LOADS DIAGRAM - NEGATIVE SIDE   
        figure(6)
%         % ---------------------------------------------------------------------
%         plot(V_from0toSinv,  LW_from0toSinv_new,  '-r', 'LineWidth', 1)
%         plot(V_fromSinvtoG1, LW_fromSinvtoG1_new, '-r', 'LineWidth', 1)
%         plot(V_fromG1toF,    LW_fromG1toF_new,    '-r', 'LineWidth', 1)
%         plot(V_fromFtoE,     LW_fromFtoE_new,     '-r', 'LineWidth', 1)
%         plot(V_fromEto0,     LW_fromEto0_new,     '-r', 'LineWidth', 1)
%         % ---------------------------------------------------------------------
%         plot(V_from0toSinv(1),    LW_from0toSinv_new(1),    'k.', 'MarkerSize', 10)
%         plot(V_from0toSinv(end),  LW_from0toSinv_new(end),  'k.', 'MarkerSize', 10)
% %         plot(VG,                  WBL_G,                    'k.', 'MarkerSize', 10)
%         plot(V_fromSinvtoG1(end), LW_fromSinvtoG1_new(end), 'k.', 'MarkerSize', 10)
%         plot(V_fromG1toF(end),    LW_fromG1toF_new(end),    'k.', 'MarkerSize', 10)
%         plot(V_fromFtoE(end),     LW_fromFtoE_new(end),     'k.', 'MarkerSize', 10)
%         plot(V_fromEto0(end),     LW_fromEto0_new(end),     'k.', 'MarkerSize', 10)
%         % ---------------------------------------------------------------------
%         text(V_from0toSinv(end),  LW_from0toSinv_new(end),  ' S inv.', 'FontSize', 6)
% %         text(VG,                  WBL_G,                    ' G',      'FontSize', 6)
%         text(V_fromSinvtoG1(end), LW_fromSinvtoG1_new(end), ' G1',     'FontSize', 6)
%         text(V_fromG1toF(end),    LW_fromG1toF_new(end),    ' F',      'FontSize', 6)
%         text(V_fromFtoE(end),     LW_fromFtoE_new(end),     ' E',      'FontSize', 6) 

        % ---------------------------------------------------------------------
        plot(V_from0toSinv,  WBL_from0toSinv,  '-r', 'LineWidth', 1)
        plot(V_fromSinvtoG1, WBL_fromSinvtoG1, '-r', 'LineWidth', 1)
        plot(V_fromG1toF,    WBL_fromG1toF,    '-r', 'LineWidth', 1)
        plot(V_fromFtoE,     WBL_fromFtoE,     '-r', 'LineWidth', 1)
        plot(V_fromEto0,     WBL_fromEto0,     '-r', 'LineWidth', 1)
        % ---------------------------------------------------------------------
        plot(V_from0toSinv(1),    WBL_from0toSinv(1),    'k.', 'MarkerSize', 10)
        plot(V_from0toSinv(end),  WBL_from0toSinv(end),  'k.', 'MarkerSize', 10)
%         plot(VG,                  WBL_G,                    'k.', 'MarkerSize', 10)
        plot(V_fromSinvtoG1(end), WBL_fromSinvtoG1(end), 'k.', 'MarkerSize', 10)
        plot(V_fromG1toF(end),    WBL_fromG1toF(end),    'k.', 'MarkerSize', 10)
        plot(V_fromFtoE(end),     WBL_fromFtoE(end),     'k.', 'MarkerSize', 10)
        plot(V_fromEto0(end),     WBL_fromEto0(end),     'k.', 'MarkerSize', 10)
        % ---------------------------------------------------------------------
        text(V_from0toSinv(end),  WBL_from0toSinv(end),  ' S inv.', 'FontSize', 6)
%         text(VG,                  WBL_G,                    ' G',      'FontSize', 6)
        text(V_fromSinvtoG1(end), WBL_fromSinvtoG1(end), ' G1',     'FontSize', 6)
        text(V_fromG1toF(end),    WBL_fromG1toF(end),    ' F',      'FontSize', 6)
        text(V_fromFtoE(end),     WBL_fromFtoE(end),     ' E',      'FontSize', 6) 
        
        ylim padded;
        xlim padded;
        
        % STORE INSIDE THE AIRCRAFT STRUCTURE VARIABLE
        
        % POINT S_INV
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.value = q_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.Attributes.unit = "Pa";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv.value = WBL_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv.Attributes.unit = "daN";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv_new.value = WBL_from0toSinv_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LHTS_inv.value = LHT_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LHTS_inv.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv.value = alfa_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_rad.value = deg2rad(alfa_from0toSinv(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_rad.Attributes.unit = "rad"; 
        alfaS_inv = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv.value;  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_new.value = alfa_from0toSinv_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_new_rad.value = deg2rad(alfa_from0toSinv_fullvehicle(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_new_rad.Attributes.unit = "rad";
        alfaS_inv_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_new.value;  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv.value = CL_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv.Attributes.unit = "Non dimensional";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv_new.value = CL_from0toSinv_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv_new.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CLHT_S_inv.value = CLHT_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CLHT_S_inv.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CD_S_inv.value = CD_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CD_S_inv.Attributes.unit = "Non dimensional";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CM_S_inv.value = CMCG_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CM_S_inv.Attributes.unit = "Non dimensional";     
        % POINT G1
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.qG1.value = q_fromSinvtoG1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.qG1.Attributes.unit = "Pa";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1.value = WBL_fromSinvtoG1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1.Attributes.unit = "daN";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1_new.value = WBL_fromSinvtoG1_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LHTG1.value = LHT_fromSinvtoG1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LHTG1.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1.value = alfa_fromSinvtoG1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_rad.value = deg2rad(alfa_fromSinvtoG1(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_rad.Attributes.unit = "rad";
        alfaG1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1.value;  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_new.value = alfa_fromSinvtoG1_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_new_rad.value = deg2rad(alfa_fromSinvtoG1_fullvehicle(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_new_rad.Attributes.unit = "rad";
        alfaG1_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_new.value;  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CL_G1.value = CL_fromSinvtoG1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CL_G1.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CL_G1_new.value = CL_fromSinvtoG1_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CL_G1_new.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CLHT_G1.value = CLHT_fromSinvtoG1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CLHT_G1.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CD_G1.value = CD_fromSinvtoG1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CD_G1.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CM_G1.value = CMCG_fromSinvtoG1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CM_G1.Attributes.unit = "Non dimensional";
        % POINT F
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value = q_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.Attributes.unit = "Pa"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.value = WBL_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.Attributes.unit = "daN";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF_new.value = WBL_fromG1toF_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.value = LHT_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.value = alfa_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_rad.value = deg2rad(alfa_fromG1toF(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_rad.Attributes.unit = "rad";
        alfaF = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.value;  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_new.value = alfa_fromG1toF_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_new_rad.value = deg2rad(alfa_fromG1toF_fullvehicle(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_new_rad.Attributes.unit = "rad";
        alfaF_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_new.value;  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value = CL_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F_new.value = CL_fromG1toF_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F_new.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CLHT_F.value = CLHT_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CLHT_F.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.value = CD_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CM_F.value = CMCG_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CM_F.Attributes.unit = "Non dimensional";
        % POINT E
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value = q_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.Attributes.unit = "Pa"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE.value = WBL_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE_new.value = WBL_fromFtoE_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.value = LHT_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.value = alfa_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_rad.value = deg2rad(alfa_fromFtoE(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_rad.Attributes.unit = "rad";
        alfaE = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_new.value = alfa_fromFtoE_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_new.Attributes.unit = "deg";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_new_rad.value = deg2rad(alfa_fromFtoE_fullvehicle(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_new_rad.Attributes.unit = "rad";
        alfaE_new = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_new.value;  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.value = CL_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.Attributes.unit = "Non dimensional";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E_new.value = CL_fromFtoE_fullvehicle(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E_new.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CLHT_E.value = CLHT_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CLHT_E.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CD_E.value = CD_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CD_E.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CM_E.value = CMCG_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CM_E.Attributes.unit = "Non dimensional";   
end
% -------------------------------------------------------------------------
% SAVING FIGURE 5 - HORIZONTAL TAIL BALANCING LOADS
% -------------------------------------------------------------------------
exportgraphics(HT_balancing_loads, 'Balancingloads.pdf', 'ContentType', 'vector')
exportgraphics(HT_balancing_loads, 'Balancingloads.png', 'ContentType', 'vector')   

% Saving figures inside correct folder
fprintf('Saving Balancingloads.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile Balancingloads.pdf Output
movefile Balancingloads.png Output 
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% SAVING FIGURE 6 - WING LOADS
% -------------------------------------------------------------------------
exportgraphics(Wing_balancing_loads, 'Wingairloads.pdf', 'ContentType', 'vector')
exportgraphics(Wing_balancing_loads, 'Wingairloads.png', 'ContentType', 'vector')   

% Saving figures inside correct folder
fprintf('Saving Wingairloads.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile Wingairloads.pdf Output
movefile Wingairloads.png Output 

Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Attributes.cs = " 421 ";
