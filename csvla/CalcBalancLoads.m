% =========================================================================
% Alpha_star and Alpha_max 
% =========================================================================
a = -0.021837866;
b = 0.436064773;
c = -0.56312855;
Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a.value = a;
Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b.value = b;
Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c.value = c;
% =========================================================================
% ZERO LIFT COEFFICIENT
CL0 = Aircraft.Certification.Aerodynamic_data.CL0.value;

% ==== USEFUL FUNCTION DEFINED LOCALLY ====
% -------------------------------------------------------------------------
% CLMAX FUNCTION
CLmax_func = @(rho, S, V, WS, n) (2 / rho) * (1 / V.^2) * (WS) * n;
% -------------------------------------------------------------------------
% CLMAX NON LINEAR
CLmax_non_lin = @(alfa) a*alfa^2 + b*alfa + c;
% -------------------------------------------------------------------------
% ALFA FUNCTION
alfa_func = @(rho, S, V, WS, n, CLalfa, alfa_0lift) (2 / rho) * (1 / V.^2) * (1/CLalfa) * (WS) * n + alfa_0lift;
% -------------------------------------------------------------------------
% GUST LOAD FACTOR - POSITIVE FLIGHT
nGust  = @(rho, V, a, kG, Ude, WS) 1 + (0.5 * rho * V * a * kG * Ude)/(WS); 
% -------------------------------------------------------------------------
% GUST LOAD FACTOR - INVERTED FLIGHT
nGust_inverted  = @(rho, V, a, kG, Ude, WS) 1 - (0.5 * rho * V * a * kG * Ude)/(WS); 
% -------------------------------------------------------------------------
% STALL SPEED FUNCTION
Vstall = @(WS, rho, CLmax, n) sqrt(WS * (2/rho) * (1/CLmax).*n); 
% -------------------------------------------------------------------------

% =========================================================================
% DRAG POLYNOMIAL COEFFICIENTS DEFINITION 
% =========================================================================
Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value = 0.079;                       % Coefficient inside an expression for the CD in polynomial form
Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value = 0.365;                       % Coefficient inside an expressione for the CD in polynomial form
Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.Attributes.unit = "Non dimensional"; 
k1 = Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value;
k2 = Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value;
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
CLalfa = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;

% =========================================================================

%% Script to evaluate balancing horizontal tail loads
%   DESCRIPTION
%    In this script, all the methods relative to the calculation of
%    aerodynamic coefficients, necessary to the correct estimation of the
%    airloads acting on the horizontal empennage. 

%% RETURN INSIDE UTILITIES
cd .. 
cd utilities
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);
obj1 = aero_model; 

%% MAX LIFT COEFFICIENT - WING BODY CONFIGURATION
% -------------------------------------------------------------------------
% Maximum lift coefficient for the wing - body configuration 
% A complete documentation of the function CLMax(...) used here is inside
% the class file aero_model.m, which can be found inside the 'utilities'
% folder of this library.
Aircraft.Certification.Aerodynamic_data.CLMAX_wingbody.value = CLMax(obj1, ...
                                                      Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a.value, ...
                                                      Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b.value, ...
                                                      Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c.value); 
Aircraft.Certification.Aerodynamic_data.CLMAX_wingbody.Attributes.unit = "Non dimensional";

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

%% LIFT CHARACTERISTIC
% Before to start the balancing loads analysis, it is required to plot lift
% coefficient values calculated with a non-linear formulation. Then, linear
% and non-linear calculation are compared with raw data available. The
% selected range is -2.0<CL<2.0. To find a complete documentation of the
% function 'CL_Non_linear_model(...)' search inside aero_model.m file. 

% NUMBER OF ELEMENTS
numb = 1e3;

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CL_WB_model = @(alpha) a*alpha.^2 + b*alpha + c; 
alpha_plus  = @(CL) (-b + sqrt(b^2 - 4*a*(c - CL)))/(2*a);
alpha_meno  = @(CL) (-b - sqrt(b^2 - 4*a*(c - CL)))/(2*a);
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Alpha zero lift
Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value = alpha_calc_lin(obj1, ...
                                                                               0.0, ...
                                                                               Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                                                               Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.Attributes.unit = "degree";

alpha_zerol = Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
alpha_star  = alpha_plus(Aircraft.Certification.Aerodynamic_data.CL_star.value);
alpha_max   = alpha_meno(Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value);

% Lift coefficient curve
Aircraft.Certification.Aerodynamic_data.AOA_aux.value = linspace(-8.0, 25*alpha_star*1e-3 - alpha_star, numb)';
Aircraft.Certification.Aerodynamic_data.AOA_aux.Attributes.unit = "degree";
Aircraft.Certification.Aerodynamic_data.AOA_aux1.value = linspace(alpha_star + 50*alpha_star*1e-3, 13.0 , numb)';
Aircraft.Certification.Aerodynamic_data.AOA_aux1.Attributes.unit = "degrees";

% Linear and non linear lift curve
Aircraft.Certification.Aerodynamic_data.CL_Linear.value = (Aircraft.Certification.Aerodynamic_data.CL0.value+0.01) + Aircraft.Certification.Aerodynamic_data.AOA_aux.value*Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;
Aircraft.Certification.Aerodynamic_data.CL_Linear.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.CL_Non_Linear.value = CL_WB_model(Aircraft.Certification.Aerodynamic_data.AOA_aux1.value);
Aircraft.Certification.Aerodynamic_data.CL_Non_Linear.Attributes.unit = "Non dimensional";

% Full lift model 
Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value = [Aircraft.Certification.Aerodynamic_data.CL_Linear.value; Aircraft.Certification.Aerodynamic_data.CL_Non_Linear.value];
Aircraft.Certification.Aerodynamic_data.CL_fullmodel.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.value = [Aircraft.Certification.Aerodynamic_data.AOA_aux.value; Aircraft.Certification.Aerodynamic_data.AOA_aux1.value];
Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.Attributes.unit = "degrees";
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
disp(" ")
disp(" ++++ FIGURE 6 - LIFT CURVE INTERPOLATION ++++ ");
Aircraft.Certification.Aerodynamic_data.CL_fullmodel_diagram.value = Lift_fullmodel_curve(Aircraft.Certification.Aerodynamic_data.AOA_aux_fullmodel.value, ...
    Aircraft.Certification.Aerodynamic_data.CL_fullmodel.value, ...
    str2num(Aircraft.Certification.Aerodynamic_data.CL.value), ...
    str2num(Aircraft.Certification.Aerodynamic_data.alpha.value));

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

% INSTANTIATION OF SOME CONSTANTS IMPORTANT TO MOMENT CONTRIB. CALCULATIONS
XAC         = Aircraft.Certification.Aerodynamic_data.XAC_nondim.value;
XCG         = Aircraft.Certification.Aerodynamic_data.XCG_nondim.value;
bCG         = Aircraft.Certification.Aerodynamic_data.bcg.value;
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
% 


%% BALANCING LOADS CALCULATIONS - POSITIVE LOAD FACTORS
% ------------------------------------------------------------------------- 
% Calculation of the CL for the wing - body configuration. It must be
% noticed that the calculations relative to lift coefficient in this
% particular range of airspeed V and load factor n are going to give a
% constant as a results, which is the maximum lift coefficient following
% the lift curve of the aircrat. The previously defined CLMAX is a little
% bit greater than the one obtained from V - n diagram flight condition to
% take into account the extra lift produced by the horizontal empennage.

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
        CMCL_from0toS = zeros(length(V_from0toS), 1);
        CMCD_from0toS = zeros(length(V_from0toS), 1);
        CMCT_from0toS = zeros(length(V_from0toS), 1);
        CMCG_from0toS = zeros(length(V_from0toS), 1); 
        CLHT_from0toS = zeros(length(V_from0toS), 1);
        LHT_from0toS  = zeros(length(V_from0toS), 1);
        for i = 1:length(V_from0toS)
            alfa_from0toS(i) = alfa_func(rho0, S, V_from0toS(i), WS, n_from0toS(i), CLalfa, alpha_zerol);
            CL_from0toS(i)   = CL_calc(obj1, n_from0toS(i), Mass, g, V_from0toS(i), rho0, S);
            if CL_from0toS(i) < CL_star
                CL_from0toS(i) = CL_calc(obj1, n_from0toS(i), Mass, g, V_from0toS(i), rho0, S);
            elseif CL_from0toS(i) > CL_star
                CL_from0toS(i) = CLmax_non_lin(alfa_from0toS(i));
            end
            CD_from0toS(i)   = cd_calc(obj1, CD0, CL_from0toS(i), AR, e, k1, k2);
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
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS.value = CL_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS.value = alfa_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS_rad.value = deg2rad(alfa_from0toS);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS_rad.Attributes.unit = "rad";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toS.value = CD_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toS.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toS.value = q_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toS.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toS.value = WBL_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toS.Attributes.unit = "daN";
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
        CMCL_fromStoA1 = zeros(length(V_fromStoA1), 1);
        CMCD_fromStoA1 = zeros(length(V_fromStoA1), 1);
        CMCT_fromStoA1 = zeros(length(V_fromStoA1), 1);
        CMCG_fromStoA1 = zeros(length(V_fromStoA1), 1); 
        CLHT_fromStoA1 = zeros(length(V_fromStoA1), 1);
        LHT_fromStoA1  = zeros(length(V_fromStoA1), 1);
        for i = 1:length(V_fromStoA1)
            alfa_fromStoA1(i) = alfa_func(rho0, S, V_fromStoA1(i), WS, n_fromStoA1(i), CLalfa, alpha_zerol);
            CL_fromStoA1(i)   = CL_calc(obj1, n_fromStoA1(i), Mass, g, V_fromStoA1(i), rho0, S);
            if CL_fromStoA1(i) < CL_star
                CL_fromStoA1(i) = CL_calc(obj1, n_fromStoA1(i), Mass, g, V_fromStoA1(i), rho0, S);
            elseif CL_fromStoA1(i) > CL_star
                CL_fromStoA1(i) = CLmax_non_lin(alfa_fromStoA1(i));
            end
            CD_fromStoA1(i)   = cd_calc(obj1, CD0, CL_fromStoA1(i), AR, e, k1, k2);
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
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1.value = CL_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1.value = alfa_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1_rad.value = deg2rad(alfa_fromStoA1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1_rad.Attributes.unit = "rad";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromStoA1.value = CD_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromStoA1.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromStoA1.value = q_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromStoA1.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromStoA1.value = WBL_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromStoA1.Attributes.unit = "daN";
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
            CMCL_fromA1toC1 = zeros(length(V_fromA1toC1), 1);
            CMCD_fromA1toC1 = zeros(length(V_fromA1toC1), 1);
            CMCT_fromA1toC1 = zeros(length(V_fromA1toC1), 1);
            CMCG_fromA1toC1 = zeros(length(V_fromA1toC1), 1); 
            CLHT_fromA1toC1 = zeros(length(V_fromA1toC1), 1);
            LHT_fromA1toC1  = zeros(length(V_fromA1toC1), 1);
            for i = 1:length(V_fromA1toC1)
                alfa_fromA1toC1(i) = alfa_func(rho0, S, V_fromA1toC1(i), WS, n_fromA1toC1(i), CLalfa, alpha_zerol);
                CL_fromA1toC1(i)   = CL_calc(obj1, n_fromA1toC1(i), Mass, g, V_fromA1toC1(i), rho0, S);
                if CL_fromA1toC1(i) < CL_star
                    CL_fromA1toC1(i) = CL_calc(obj1, n_fromA1toC1(i), Mass, g, V_fromA1toC1(i), rho0, S);
                elseif CL_fromA1toC1(i) > CL_star
                    CL_fromA1toC1(i) = CLmax_non_lin(alfa_fromA1toC1(i));
                end
                CD_fromA1toC1(i)   = cd_calc(obj1, CD0, CL_fromA1toC1(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC1.value = CL_fromA1toC1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC1.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC1.value = alfa_fromA1toC1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC1.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC1_rad.value = deg2rad(alfa_fromA1toC1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC1_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromA1toC1.value = CD_fromA1toC1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromA1toC1.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromA1toC1.value = q_fromA1toC1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromA1toC1.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromA1toC1.value = WBL_fromA1toC1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromA1toC1.Attributes.unit = "daN";
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
            CMCL_fromC1toC = zeros(length(V_fromC1toC), 1);
            CMCD_fromC1toC = zeros(length(V_fromC1toC), 1);
            CMCT_fromC1toC = zeros(length(V_fromC1toC), 1);
            CMCG_fromC1toC = zeros(length(V_fromC1toC), 1); 
            CLHT_fromC1toC = zeros(length(V_fromC1toC), 1);
            LHT_fromC1toC  = zeros(length(V_fromC1toC), 1);
            for i = 1:length(V_fromC1toC)
                alfa_fromC1toC(i) = alfa_func(rho0, S, V_fromC1toC(i), WS, n_fromC1toC(i), CLalfa, alpha_zerol);
                CL_fromC1toC(i)   = CL_calc(obj1, n_fromC1toC(i), Mass, g, V_fromC1toC(i), rho0, S);
                if CL_fromC1toC(i) < CL_star
                    CL_fromC1toC(i) = CL_calc(obj1, n_fromC1toC(i), Mass, g, V_fromC1toC(i), rho0, S);
                elseif CL_fromC1toC(i) > CL_star
                    CL_fromC1toC(i) = CLmax_non_lin(alfa_fromC1toC(i));
                end
                CD_fromC1toC(i)   = cd_calc(obj1, CD0, CL_fromC1toC(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC1toC.value = CL_fromC1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC1toC.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC1toC.value = alfa_fromC1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC1toC.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC1toC_rad.value = deg2rad(alfa_fromC1toC);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC1toC_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromC1toC.value = CD_fromC1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromC1toC.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromC1toC.value = q_fromC1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromC1toC.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromC1toC.value = WBL_fromC1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromC1toC.Attributes.unit = "daN";
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
            CMCL_fromCtoC2 = zeros(length(V_fromCtoC2), 1);
            CMCD_fromCtoC2 = zeros(length(V_fromCtoC2), 1);
            CMCT_fromCtoC2 = zeros(length(V_fromCtoC2), 1);
            CMCG_fromCtoC2 = zeros(length(V_fromCtoC2), 1); 
            CLHT_fromCtoC2 = zeros(length(V_fromCtoC2), 1);
            LHT_fromCtoC2  = zeros(length(V_fromCtoC2), 1);
            for i = 1:length(V_fromCtoC2)
                alfa_fromCtoC2(i) = alfa_func(rho0, S, V_fromCtoC2(i), WS, n_fromCtoC2(i), CLalfa, alpha_zerol);
                CL_fromCtoC2(i)   = CL_calc(obj1, n_fromCtoC2(i), Mass, g, V_fromCtoC2(i), rho0, S);
                if CL_fromCtoC2(i) < CL_star
                    CL_fromCtoC2(i) = CL_calc(obj1, n_fromCtoC2(i), Mass, g, V_fromCtoC2(i), rho0, S);
                elseif CL_fromCtoC2(i) > CL_star
                    CL_fromCtoC2(i) = CLmax_non_lin(alfa_fromCtoC2(i));
                end
                CD_fromCtoC2(i)   = cd_calc(obj1, CD0, CL_fromCtoC2(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoC2.value = CL_fromCtoC2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoC2.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoC2.value = alfa_fromCtoC2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoC2.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoC2_rad.value = deg2rad(alfa_fromCtoC2);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoC2_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtoC2.value = CD_fromCtoC2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtoC2.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromCtoC2.value = q_fromCtoC2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromCtoC2.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtoC2.value = WBL_fromCtoC2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtoC2.Attributes.unit = "daN";
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
            CMCL_fromC2toD = zeros(length(V_fromC2toD), 1);
            CMCD_fromC2toD = zeros(length(V_fromC2toD), 1);
            CMCT_fromC2toD = zeros(length(V_fromC2toD), 1);
            CMCG_fromC2toD = zeros(length(V_fromC2toD), 1); 
            CLHT_fromC2toD = zeros(length(V_fromC2toD), 1);
            LHT_fromC2toD  = zeros(length(V_fromC2toD), 1);
            for i = 1:length(V_fromC2toD)
                alfa_fromC2toD(i) = alfa_func(rho0, S, V_fromC2toD(i), WS, n_fromC2toD(i), CLalfa, alpha_zerol);
                CL_fromC2toD(i)   = CL_calc(obj1, n_fromC2toD(i), Mass, g, V_fromC2toD(i), rho0, S);
                if CL_fromC2toD(i) < CL_star
                    CL_fromC2toD(i) = CL_calc(obj1, n_fromC2toD(i), Mass, g, V_fromC2toD(i), rho0, S);
                elseif CL_fromC2toD(i) > CL_star
                    CL_fromC2toD(i) = CLmax_non_lin(alfa_fromC2toD(i));
                end
                CD_fromC2toD(i)   = cd_calc(obj1, CD0, CL_fromC2toD(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC2toD.value = CL_fromC2toD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC2toD.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC2toD.value = alfa_fromC2toD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC2toD.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC2toD_rad.value = deg2rad(alfa_fromC2toD);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromC2toD_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromC2toD.value = CD_fromC2toD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromC2toD.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromC2toD.value = q_fromC2toD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromC2toD.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromC2toD.value = WBL_fromC2toD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromC2toD.Attributes.unit = "daN";
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
            CMCL_fromDto0 = zeros(length(V_fromDto0), 1);
            CMCD_fromDto0 = zeros(length(V_fromDto0), 1);
            CMCT_fromDto0 = zeros(length(V_fromDto0), 1);
            CMCG_fromDto0 = zeros(length(V_fromDto0), 1); 
            CLHT_fromDto0 = zeros(length(V_fromDto0), 1);
            LHT_fromDto0  = zeros(length(V_fromDto0), 1);
            for i = 1:length(V_fromDto0)
                alfa_fromDto0(i) = alfa_func(rho0, S, V_fromDto0(i), WS, n_fromDto0(i), CLalfa, alpha_zerol);
                CL_fromDto0(i)   = CL_calc(obj1, n_fromDto0(i), Mass, g, V_fromDto0(i), rho0, S);
                if CL_fromDto0(i) < CL_star
                    CL_fromDto0(i) = CL_calc(obj1, n_fromDto0(i), Mass, g, V_fromDto0(i), rho0, S);
                elseif CL_fromDto0(i) > CL_star
                    CL_fromDto0(i) = CLmax_non_lin(alfa_fromDto0(i));
                end
                CD_fromDto0(i)   = cd_calc(obj1, CD0, CL_fromDto0(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0.value = CL_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0.value = alfa_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_rad.value = deg2rad(alfa_fromDto0);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDto0.value = CD_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromDto0.value = q_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromDto0.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDto0.value = WBL_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDto0.Attributes.unit = "daN";
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
            end

            % WING BODY LIFT COEFFICIENT
            Aircraft_weight = ones(length(V_unit_load_factor), 1)*Mass*g;
            CLWB_unit_load_factor = zeros(length(V_unit_load_factor), 1);
            for i = 1:length(V_unit_load_factor)
                Aero_force = q_unit_load_factor(i)*S;
                CLWB_unit_load_factor(i) = ((Aircraft_weight(i))*(1/Aero_force));
            end

            % DRAG LIFT COEFFICIENT
            CD_unit_load_factor = zeros(length(V_unit_load_factor), 1);
            for i = 1:length(V_unit_load_factor)
                CD_unit_load_factor(i) = cd_calc(obj1, CD0, CLWB_unit_load_factor(i), ...
                                                 AR, e, k1, k2);
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value = V_unit_load_factor;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.Attributes.unit = "m/s";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.value = q_unit_load_factor;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value           = CLWB_unit_load_factor;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value = CD_unit_load_factor;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.Attributes.unit = "Non dimensional"; 

            % Interpolation coefficient from actual aerodynamic data
            CL_supp    = str2num(Aircraft.Certification.Aerodynamic_data.CL.value);
            alpha_supp = str2num(Aircraft.Certification.Aerodynamic_data.alpha.value);
            x          = 0.027*ones(length(CL_supp), 1);
            p = polyfit(CL_supp + x, alpha_supp, 2);
            CLalfa_deg = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;

            % ALPHA CALCULATION 
            alpha_unit_load_factor = zeros(length(V_unit_load_factor), 1);
            for i = 1:length(alpha_unit_load_factor)
                alpha_unit_load_factor(i) = alpha_calc(obj1, ...
                                                       CLWB_unit_load_factor(i), CL0, CL_star, CLalfa_deg, p);
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.value = alpha_unit_load_factor;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value = deg2rad(alpha_unit_load_factor);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.Attributes.unit = "radians"; 

            % WING BODY LIFT 
            WBLift_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
            for i = 1:length(alpha_unit_load_factor)   
                WBLift_unit_load_factor(i) = (0.5) * (V_unit_load_factor(i)^2) * rho0 * S * (CLWB_unit_load_factor(i))*(1E-1);   
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.value = WBLift_unit_load_factor;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.Attributes.unit = "daN";       

            % EVALUATION OF TAIL LOADS CONTRIBUTION     
            CMCL_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
            CMCD_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
            CMCT_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
            CMCG_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
            CLHT_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
            LHT_unit_load_factor  = zeros(length(alpha_unit_load_factor), 1);
            LW_unit_load_factor   = zeros(length(alpha_unit_load_factor), 1);
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
                LHT_unit_load_factor(i) = (0.5)*(V_unit_load_factor(i)^2)* S * rho0 * (CLHT_unit_load_factor(i))*(1e-1);
                LW_unit_load_factor(i)  = WBLift_unit_load_factor(i) - LHT_unit_load_factor(i);                                                                                     
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
            text(V_fromCtoC2(end), V_fromCtoC2(end),     '  C2', 'FontSize', 6)
            text(V_fromC2toD(end), V_fromC2toD(end),     '  D',  'FontSize', 6)
            text(40.75, -18, 'n = 1')
            % ---------------------------------------------------------------------
            xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
            ylabel("Horizontal tail lift - $L_{ht}$ (daN)", "Interpreter", "latex")
            title("Horizontal empennage airloads per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct      
            
            % MAIN WING LOADS DIAGRAM        
            CL_from0toS_new   = CL_from0toS   - CLHT_from0toS;        
            CL_fromStoA1_new  = CL_fromStoA1  - CLHT_fromStoA1;        
            CL_fromA1toC1_new = CL_fromA1toC1 - CLHT_fromA1toC1;        
            CL_fromC1toC_new  = CL_fromC1toC  - CLHT_fromC1toC;      
            CL_fromCtoC2_new  = CL_fromCtoC2  - CLHT_fromCtoC2;
            CL_fromC2toD_new  = CL_fromC2toD  - CLHT_fromC2toD;
            CL_fromDto0_new   = CL_fromDto0   - CLHT_fromDto0;
            LW_from0toS_new   = zeros(length(CL_from0toS_new), 1);
            LW_fromStoA1_new  = zeros(length(CL_from0toS_new), 1);
            LW_fromA1toC1_new = zeros(length(CL_from0toS_new), 1);
            LW_fromC1toC_new  = zeros(length(CL_from0toS_new), 1);
            LW_fromCtoC2_new  = zeros(length(CL_from0toS_new), 1);
            LW_fromC2toD_new  = zeros(length(CL_from0toS_new), 1);
            LW_fromDto0_new   = zeros(length(CL_from0toS_new), 1);
            for i = 1:length(CL_from0toS_new)
                LW_from0toS_new(i)   = (0.5)*(V_from0toS(i)^2)* rho0 * S *(CL_from0toS_new(i))*(1e-1);  
                LW_fromStoA1_new(i)  = (0.5)*(V_fromStoA1(i)^2)* rho0 * S *(CL_fromStoA1_new(i))*(1e-1);
                LW_fromA1toC1_new(i) = (0.5)*(V_fromA1toC1(i)^2)* rho0 * S *(CL_fromA1toC1_new(i))*(1e-1);
                LW_fromC1toC_new(i)  = (0.5)*(V_fromC1toC(i)^2)* rho0 * S *(CL_fromC1toC_new(i))*(1e-1);
                LW_fromCtoC2_new(i)  = (0.5)*(V_fromCtoC2(i)^2)* rho0 * S *(CL_fromCtoC2_new(i))*(1e-1);
                LW_fromC2toD_new(i)  = (0.5)*(V_fromC2toD(i)^2)* rho0 * S *(CL_fromC2toD_new(i))*(1e-1);
                LW_fromDto0_new(i)   = (0.5)*(V_fromDto0(i)^2)* rho0 * S *(CL_fromDto0_new(i))*(1e-1);
            end       
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS_new.value = CL_from0toS_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1_new.value = CL_fromStoA1_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC1_new.value = CL_fromA1toC1_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC1_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC1toC_new.value = CL_fromC1toC_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC1toC_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoC2_new.value = CL_fromCtoC2_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoC2_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC2toD_new.value = CL_fromC2toD_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromC2toD_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.value = CL_fromDto0_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.Attributes.unit = "Non dimensional";
            % ------------------------------------------------------------------------------------------------------------------------
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.value = LW_from0toS_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromStoA1_new.value = LW_fromStoA1_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromStoA1_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC1_new.value = LW_fromA1toC1_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC1_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromC1toC_new.value = LW_fromC1toC_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromC1toC_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoC2_new.value = LW_fromCtoC2_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoC2_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromC2toD_new.value = LW_fromC2toD_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromC2toD_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.value = LW_fromDto0_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.Attributes.unit = "daN";  
            
            % TAIL LOADS DIAGRAM - POSITIVE SIDE   
            Wing_balancing_loads = figure(6); 
            hold on; grid on; grid minor; 

    %         ylim([LHT_from0toS(1) LHT_fromDto0(1)])
    %         xlim([0.0 V_fromDto0(1)])

            plot(V_from0toS,   LW_from0toS_new,   '-r', 'LineWidth', 1)
            plot(V_fromStoA1,  LW_fromStoA1_new,  '-r', 'LineWidth', 1)
            plot(V_fromA1toC1, LW_fromA1toC1_new, '-r', 'LineWidth', 1)
            plot(V_fromC1toC,  LW_fromC1toC_new,  '-r', 'LineWidth', 1)
            plot(V_fromCtoC2,  LW_fromCtoC2_new,  '-r', 'LineWidth', 1)
            plot(V_fromC2toD,  LW_fromC2toD_new,  '-r', 'LineWidth', 1)
            plot(V_fromDto0,   LW_fromDto0_new,   '-r', 'LineWidth', 1)
            % UNIT LOAD DISTR.
            plot(V_unit_load_factor, LW_unit_load_factor,  '--k', 'LineWidth', 2)
            % ---------------------------------------------------------------------
            plot(V_from0toS(1),     LW_from0toS_new(1),     'k.', 'MarkerSize', 10)
            plot(V_from0toS(end),   LW_from0toS_new(end),   'k.', 'MarkerSize', 10)
            plot(V_fromStoA1(end),  LW_fromStoA1_new(end),  'k.', 'MarkerSize', 10)
            plot(V_fromA1toC1(end), LW_fromA1toC1_new(end), 'k.', 'MarkerSize', 10)
            plot(V_fromC1toC(end),  LW_fromC1toC_new(end),  'k.', 'MarkerSize', 10)
            plot(V_fromCtoC2(end),  LW_fromCtoC2_new(end),  'k.', 'MarkerSize', 10)
            plot(V_fromC2toD(end),  LW_fromC2toD_new(end),  'k.', 'MarkerSize', 10)
            plot(V_fromDto0(end),   LW_fromDto0_new(end),   'k.', 'MarkerSize', 10)
            % ---------------------------------------------------------------------
            text(V_fromStoA1(1),   LW_fromStoA1_new(1),     '  S',  'FontSize', 6)
            text(V_fromStoA1(end), LW_fromStoA1_new(end),   '  A1', 'FontSize', 6)
            text(V_fromA1toC1(end), LW_fromA1toC1_new(end), '  C1', 'FontSize', 6)
            text(V_fromC1toC(end),  LW_fromC1toC_new(end),  '  C',  'FontSize', 6)
            text(V_fromCtoC2(end),  LW_fromCtoC2_new(end),  '  C2', 'FontSize', 6)
            text(V_fromC2toD(end),  LW_fromC2toD_new(end),  '  D',  'FontSize', 6)
            % text(40.75, -18, 'n = 1')
            % ---------------------------------------------------------------------
            xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
            ylabel("Full body lift - $L_{ht}$ (daN)", "Interpreter", "latex")
            title("Full body airloads per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct 
            
            % STORE INSIDE THE AIRCRAFT STRUCTURE VARIABLE

            % POINT S
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.value = q_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS_new.value = LW_from0toS_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTS.value = LHT_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTS.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.value = alfa_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_rad.value = deg2rad(alfa_from0toS(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_rad.Attributes.unit = "rad"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CLHT_S.value = CLHT_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CLHT_S.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CD_S.value = CD_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CD_S.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CM_S.value = CMCG_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CM_S.Attributes.unit = "Non dimensional";
            % POINT A1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.value = q_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1_new.value = LW_fromStoA1_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.value = LHT_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.value = alfa_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_rad.value = deg2rad(alfa_fromStoA1(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_rad.Attributes.unit = "rad"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CLHT_A1.value = CLHT_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CLHT_A1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CD_A1.value = CD_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CD_A1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CM_A1.value = CMCG_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CM_A1.Attributes.unit = "Non dimensional";
            % POINT C1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.qC1.value = q_fromA1toC1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.qC1.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.LC1_new.value = LW_fromA1toC1_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.LC1_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.LHTC1.value = LHT_fromA1toC1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.LHTC1.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.alfaC1.value = alfa_fromA1toC1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.alfaC1.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.alfaC1_rad.value = deg2rad(alfa_fromA1toC1(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.alfaC1_rad.Attributes.unit = "rad"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CLHT_C1.value = CLHT_fromA1toC1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CLHT_C1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CD_C1.value = CD_fromA1toC1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CD_C1.Attributes.unit = "Non dimensional";     
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CM_C1.value = CMCG_fromA1toC1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CM_C1.Attributes.unit = "Non dimensional";    
            % POINT C
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value = q_fromC1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.value = LW_fromC1toC_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.value = LHT_fromC1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value = alfa_fromC1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_rad.value = deDrad(alfa_fromC1toC(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_rad.Attributes.unit = "rad";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CLHT_C.value = CLHT_fromC1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CLHT_C.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CD_C.value = CD_fromC1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CD_C.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CM_C.value = CMCG_fromC1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CM_C.Attributes.unit = "Non dimensional";              
            % POINT C2
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.qC2.value = q_fromCtoC2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.qC2.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.LC2_new.value = LW_fromCtoC2_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.LC2_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.LHTC2.value = LHT_fromCtoC2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.LHTC2.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.alfaC2.value = alfa_fromCtoC2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.alfaC2.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.alfaC2_rad.value = deDrad(alfa_fromCtoC2(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.alfaC2_rad.Attributes.unit = "rad"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CLHT_C2.value = CLHT_fromCtoC2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CLHT_C2.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CD_C2.value = CD_fromCtoC2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CD_C2.Attributes.unit = "Non dimensional";  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CM_C2.value = CMCG_fromCtoC2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CM_C2.Attributes.unit = "Non dimensional";
            % POINT D
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value = q_fromC2toD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.value = LW_fromC2toD_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value = LHT_fromC2toD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value = alfa_fromC2toD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_rad.value = deDrad(alfa_fromC2toD(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_rad.Attributes.unit = "rad";  
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
            CMCL_fromA1toC = zeros(length(V_fromA1toC), 1);
            CMCD_fromA1toC = zeros(length(V_fromA1toC), 1);
            CMCT_fromA1toC = zeros(length(V_fromA1toC), 1);
            CMCG_fromA1toC = zeros(length(V_fromA1toC), 1); 
            CLHT_fromA1toC = zeros(length(V_fromA1toC), 1);
            LHT_fromA1toC  = zeros(length(V_fromA1toC), 1);
            for i = 1:length(V_fromA1toC)
                alfa_fromA1toC(i) = alfa_func(rho0, S, V_fromA1toC(i), WS, n_fromA1toC(i), CLalfa, alpha_zerol);
                CL_fromA1toC(i)   = CL_calc(obj1, n_fromA1toC(i), Mass, g, V_fromA1toC(i), rho0, S);
                if CL_fromA1toC(i) < CL_star
                    CL_fromA1toC(i) = CL_calc(obj1, n_fromA1toC(i), Mass, g, V_fromA1toC(i), rho0, S);
                elseif CL_fromA1toC(i) > CL_star
                    CL_fromA1toC(i) = CLmax_non_lin(alfa_fromA1toC(i));
                end
                CD_fromA1toC(i)   = cd_calc(obj1, CD0, CL_fromA1toC(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC.value = CL_fromA1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC.value = alfa_fromA1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC_rad.value = deg2rad(alfa_fromA1toC);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromA1toC.value = CD_fromA1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromA1toC.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromA1toC.value = q_fromA1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromA1toC.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromA1toC.value = WBL_fromA1toC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromA1toC.Attributes.unit = "daN";
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
            CMCL_fromCtoD = zeros(length(V_fromCtoD), 1);
            CMCD_fromCtoD = zeros(length(V_fromCtoD), 1);
            CMCT_fromCtoD = zeros(length(V_fromCtoD), 1);
            CMCG_fromCtoD = zeros(length(V_fromCtoD), 1); 
            CLHT_fromCtoD = zeros(length(V_fromCtoD), 1);
            LHT_fromCtoD  = zeros(length(V_fromCtoD), 1);
            for i = 1:length(V_fromCtoD)
                alfa_fromCtoD(i) = alfa_func(rho0, S, V_fromCtoD(i), WS, n_fromCtoD(i), CLalfa, alpha_zerol);
                CL_fromCtoD(i)   = CL_calc(obj1, n_fromCtoD(i), Mass, g, V_fromCtoD(i), rho0, S);
                if CL_fromCtoD(i) < CL_star
                    CL_fromCtoD(i) = CL_calc(obj1, n_fromCtoD(i), Mass, g, V_fromCtoD(i), rho0, S);
                elseif CL_fromCtoD(i) > CL_star
                    CL_fromCtoD(i) = CLmax_non_lin(alfa_fromCtoD(i));
                end
                CD_fromCtoD(i)   = cd_calc(obj1, CD0, CL_fromCtoD(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoD.value = CL_fromCtoD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoD.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoD.value = alfa_fromCtoD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoD.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoD_rad.value = deg2rad(alfa_fromCtoD);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoD_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtoD.value = CD_fromCtoD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtoD.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromCtoD.value = q_fromCtoD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromCtoD.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtoD.value = WBL_fromCtoD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtoD.Attributes.unit = "daN";
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
            CMCL_fromDto0 = zeros(length(V_fromDto0), 1);
            CMCD_fromDto0 = zeros(length(V_fromDto0), 1);
            CMCT_fromDto0 = zeros(length(V_fromDto0), 1);
            CMCG_fromDto0 = zeros(length(V_fromDto0), 1); 
            CLHT_fromDto0 = zeros(length(V_fromDto0), 1);
            LHT_fromDto0  = zeros(length(V_fromDto0), 1);
            for i = 1:length(V_fromDto0)
                alfa_fromDto0(i) = alfa_func(rho0, S, V_fromDto0(i), WS, n_fromDto0(i), CLalfa, alpha_zerol);
                CL_fromDto0(i)   = CL_calc(obj1, n_fromDto0(i), Mass, g, V_fromDto0(i), rho0, S);
                if CL_fromDto0(i) < CL_star
                    CL_fromDto0(i) = CL_calc(obj1, n_fromDto0(i), Mass, g, V_fromDto0(i), rho0, S);
                elseif CL_fromDto0(i) > CL_star
                    CL_fromDto0(i) = CLmax_non_lin(alfa_fromDto0(i));
                end
                CD_fromDto0(i)   = cd_calc(obj1, CD0, CL_fromDto0(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0.value = CL_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0.value = alfa_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_rad.value = deg2rad(alfa_fromDto0);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDto0.value = CD_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromDto0.value = q_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromDto0.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDto0.value = WBL_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDto0.Attributes.unit = "daN";
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
            end

            % WING BODY LIFT COEFFICIENT
            Aircraft_weight = ones(length(V_unit_load_factor), 1)*Mass*g;
            CLWB_unit_load_factor = zeros(length(V_unit_load_factor), 1);
            for i = 1:length(V_unit_load_factor)
                Aero_force = q_unit_load_factor(i)*S;
                CLWB_unit_load_factor(i) = ((Aircraft_weight(i))*(1/Aero_force));
            end

            % DRAG LIFT COEFFICIENT
            CD_unit_load_factor = zeros(length(V_unit_load_factor), 1);
            for i = 1:length(V_unit_load_factor)
                CD_unit_load_factor(i) = cd_calc(obj1, CD0, CLWB_unit_load_factor(i), ...
                                                 AR, e, k1, k2);
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value = V_unit_load_factor;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.Attributes.unit = "m/s";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.value = q_unit_load_factor;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value           = CLWB_unit_load_factor;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value = CD_unit_load_factor;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.Attributes.unit = "Non dimensional"; 

            % Interpolation coefficient from actual aerodynamic data
            CL_supp    = str2num(Aircraft.Certification.Aerodynamic_data.CL.value);
            alpha_supp = str2num(Aircraft.Certification.Aerodynamic_data.alpha.value);
            x          = 0.027*ones(length(CL_supp), 1);
            p = polyfit(CL_supp + x, alpha_supp, 2);
            CLalfa_deg = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;

            % ALPHA CALCULATION 
            alpha_unit_load_factor = zeros(length(V_unit_load_factor), 1);
            for i = 1:length(alpha_unit_load_factor)
                alpha_unit_load_factor(i) = alpha_calc(obj1, ...
                                                       CLWB_unit_load_factor(i), CL0, CL_star, CLalfa_deg, p);
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.value = alpha_unit_load_factor;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value = deg2rad(alpha_unit_load_factor);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.Attributes.unit = "radians"; 

            % WING BODY LIFT 
            WBLift_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
            for i = 1:length(alpha_unit_load_factor)   
                WBLift_unit_load_factor(i) = (0.5) * (V_unit_load_factor(i)^2) * rho0 * S * (CLWB_unit_load_factor(i))*(1E-1);   
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.value = WBLift_unit_load_factor;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.Attributes.unit = "daN";       

            % EVALUATION OF TAIL LOADS CONTRIBUTION     
            CMCL_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
            CMCD_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
            CMCT_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
            CMCG_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
            CLHT_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
            LHT_unit_load_factor  = zeros(length(alpha_unit_load_factor), 1);
            LW_unit_load_factor   = zeros(length(alpha_unit_load_factor), 1);
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
                LHT_unit_load_factor(i) = (0.5)*(V_unit_load_factor(i)^2)* S * rho0 * (CLHT_unit_load_factor(i))*(1e-1);
                LW_unit_load_factor(i)  = WBLift_unit_load_factor(i) - LHT_unit_load_factor(i);                                                                                     
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
            
            % MAIN WING LOADS DIAGRAM        
            CL_from0toS_new  = CL_from0toS  - CLHT_from0toS;        
            CL_fromStoA1_new = CL_fromStoA1 - CLHT_fromStoA1;        
            CL_fromA1toC_new = CL_fromA1toC - CLHT_fromA1toC;        
            CL_fromCtoD_new  = CL_fromCtoD  - CLHT_fromCtoD;      
            CL_fromDto0_new  = CL_fromDto0  - CLHT_fromDto0;
            LW_from0toS_new  = zeros(length(CL_from0toS_new), 1);
            LW_fromStoA1_new = zeros(length(CL_from0toS_new), 1);
            LW_fromA1toC_new = zeros(length(CL_from0toS_new), 1);
            LW_fromCtoD_new  = zeros(length(CL_from0toS_new), 1);
            LW_fromDto0_new  = zeros(length(CL_from0toS_new), 1);
            for i = 1:length(CL_from0toS_new)
                LW_from0toS_new(i)  = (0.5)*(V_from0toS(i)^2)* rho0 * S *(CL_from0toS_new(i))*(1e-1);  
                LW_fromStoA1_new(i) = (0.5)*(V_fromStoA1(i)^2)* rho0 * S *(CL_fromStoA1_new(i))*(1e-1);
                LW_fromA1toC_new(i) = (0.5)*(V_fromA1toC(i)^2)* rho0 * S *(CL_fromA1toC_new(i))*(1e-1);
                LW_fromCtoD_new(i) = (0.5)*(V_fromCtoD(i)^2)* rho0 * S *(CL_fromCtoD_new(i))*(1e-1);
                LW_fromDto0_new(i)  = (0.5)*(V_fromDto0(i)^2)* rho0 * S *(CL_fromDto0_new(i))*(1e-1);
            end       
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS_new.value = CL_from0toS_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1_new.value = CL_fromStoA1_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC_new.value = CL_fromA1toC_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoD_new.value = CL_fromCtoD_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoD_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.value = CL_fromDto0_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.Attributes.unit = "Non dimensional";
            % ------------------------------------------------------------------------------------------------------------------------
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.value = LW_from0toS_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromStoA1_new.value = LW_fromStoA1_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromStoA1_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC_new.value = LW_fromA1toC_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoD_new.value = LW_fromCtoD_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoD_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.value = LW_fromDto0_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.Attributes.unit = "daN";    
            
            % TAIL LOADS DIAGRAM - POSITIVE SIDE   
            Wing_balancing_loads = figure(6); 
            hold on; grid on; grid minor; 

    %         ylim([LHT_from0toS(1) LHT_fromDto0(1)])
    %         xlim([0.0 V_fromDto0(1)])

            plot(V_from0toS,  LW_from0toS_new,  '-r', 'LineWidth', 1)
            plot(V_fromStoA1, LW_fromStoA1_new, '-r', 'LineWidth', 1)
            plot(V_fromA1toC, LW_fromA1toC_new, '-r', 'LineWidth', 1)
            plot(V_fromCtoD,  LW_fromCtoD_new, '-r', 'LineWidth', 1)
            plot(V_fromDto0,  LW_fromDto0_new,  '-r', 'LineWidth', 1)
            % UNIT LOAD DISTR.
            plot(V_unit_load_factor, LW_unit_load_factor,  '--k', 'LineWidth', 2)
            % ---------------------------------------------------------------------
            plot(V_from0toS(1),    LW_from0toS_new(1),    'k.', 'MarkerSize', 10)
            plot(V_from0toS(end),  LW_from0toS_new(end),  'k.', 'MarkerSize', 10)
            plot(V_fromStoA1(end), LW_fromStoA1_new(end), 'k.', 'MarkerSize', 10)
            plot(V_fromA1toC(end), LW_fromA1toC_new(end), 'k.', 'MarkerSize', 10)
            plot(V_fromCtoD(end),  LW_fromCtoD_new(end),  'k.', 'MarkerSize', 10)
            plot(V_fromDto0(end),  LW_fromDto0_new(end),  'k.', 'MarkerSize', 10)
            % ---------------------------------------------------------------------
            text(V_fromStoA1(1),   LW_fromStoA1_new(1),   '  S',  'FontSize', 6)
            text(V_fromStoA1(end), LW_fromStoA1_new(end), '  A1', 'FontSize', 6)
            text(V_fromA1toC(end), LW_fromA1toC_new(end), '  C',  'FontSize', 6)
            text(V_fromCtoD(end),  LW_fromCtoD_new(end),  '  D', 'FontSize', 6)
            % text(40.75, -18, 'n = 1')
            % ---------------------------------------------------------------------
            xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
            ylabel("Full body lift - $L_{ht}$ (daN)", "Interpreter", "latex")
            title("Full body airloads per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct  
            
            % STORE INSIDE THE AIRCRAFT STRUCTURE VARIABLE

            % POINT S
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.value = q_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS_new.value = LW_from0toS_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTS.value = LHT_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTS.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.value = alfa_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_rad.value = deg2rad(alfa_from0toS(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_rad.Attributes.unit = "rad"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CLHT_S.value = CLHT_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CLHT_S.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CD_S.value = CD_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CD_S.Attributes.unit = "Non dimensional";  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CM_S.value = CMCG_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CM_S.Attributes.unit = "Non dimensional";     
            % POINT A1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.value = q_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1_new.value = LW_fromStoA1_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.value = LHT_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.value = alfa_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_rad.value = deg2rad(alfa_fromStoA1(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_rad.Attributes.unit = "rad";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CLHT_A1.value = CLHT_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CLHT_A1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CD_A1.value = CD_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CD_A1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CM_A1.value = CMCG_fromStoA1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CM_A1.Attributes.unit = "Non dimensional"; 
            % POINT C
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value = q_fromA1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.value = LW_fromA1toC_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.value = LHT_fromA1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value = alfa_fromA1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_rad.value = deg2rad(alfa_fromA1toC(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_rad.Attributes.unit = "rad"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CLHT_C.value = CLHT_fromA1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CLHT_C.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CD_C.value = CD_fromA1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CD_C.Attributes.unit = "Non dimensional";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CM_C.value = CMCG_fromA1toC(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CM_C.Attributes.unit = "Non dimensional";        
            % POINT D
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value = q_fromCtoD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.Attributes.unit = "Pa";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.value = LW_fromCtoD_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value = LHT_fromCtoD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value = alfa_fromCtoD(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_rad.value = deDrad(alfa_fromCtoD(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_rad.Attributes.unit = "rad";  
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
        CL_from0toS   = zeros(length(V_from0toS), 1);
        alfa_from0toS = zeros(length(V_from0toS), 1);
        CD_from0toS   = zeros(length(V_from0toS), 1);
        q_from0toS    = zeros(length(V_from0toS), 1);
        WBL_from0toS  = zeros(length(V_from0toS), 1);
        CMCL_from0toS = zeros(length(V_from0toS), 1);
        CMCD_from0toS = zeros(length(V_from0toS), 1);
        CMCT_from0toS = zeros(length(V_from0toS), 1);
        CMCG_from0toS = zeros(length(V_from0toS), 1); 
        CLHT_from0toS = zeros(length(V_from0toS), 1);
        LHT_from0toS  = zeros(length(V_from0toS), 1);
        for i = 1:length(V_from0toS)
            alfa_from0toS(i) = alfa_func(rho0, S, V_from0toS(i), WS, n_from0toS(i), CLalfa, alpha_zerol);
            CL_from0toS(i)   = CL_calc(obj1, n_from0toS(i), Mass, g, V_from0toS(i), rho0, S);
            if CL_from0toS(i) < CL_star
                CL_from0toS(i) = CL_calc(obj1, n_from0toS(i), Mass, g, V_from0toS(i), rho0, S);
            elseif CL_from0toS(i) > CL_star
                CL_from0toS(i) = CLmax_non_lin(alfa_from0toS(i));
            end
            CD_from0toS(i)   = cd_calc(obj1, CD0, CL_from0toS(i), AR, e, k1, k2);
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
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS.value = CL_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS.value = alfa_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS_rad.value = deg2rad(alfa_from0toS);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toS_rad.Attributes.unit = "rad";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toS.value = CD_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toS.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toS.value = q_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toS.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toS.value = WBL_from0toS;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toS.Attributes.unit = "daN";
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
        CMCL_fromStoA1 = zeros(length(V_fromStoA1), 1);
        CMCD_fromStoA1 = zeros(length(V_fromStoA1), 1);
        CMCT_fromStoA1 = zeros(length(V_fromStoA1), 1);
        CMCG_fromStoA1 = zeros(length(V_fromStoA1), 1);    
        CLHT_fromStoA1 = zeros(length(V_fromStoA1), 1); 
        LHT_fromStoA1  = zeros(length(V_fromStoA1), 1);    
        for i = 1:length(V_fromStoA1)
            alfa_fromStoA1(i) = alfa_func(rho0, S, V_fromStoA1(i), WS, n_fromStoA1(i), CLalfa, alpha_zerol);
            CL_fromStoA1(i)   = CL_calc(obj1, n_fromStoA1(i), Mass, g, V_fromStoA1(i), rho0, S);
            if CL_fromStoA1(i) < CL_star
                CL_fromStoA1(i) = CL_calc(obj1, n_fromStoA1(i), Mass, g, V_fromStoA1(i), rho0, S);
            elseif CL_fromStoA1(i) > CL_star
                CL_fromStoA1(i) = CLmax_non_lin(alfa_fromStoA1(i));
            end
            CD_fromStoA1(i)   = cd_calc(obj1, CD0, CL_fromStoA1(i), AR, e, k1, k2);
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
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1.value = CL_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1.value = alfa_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1_rad.value = deg2rad(alfa_fromStoA1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromStoA1_rad.Attributes.unit = "rad";         
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromStoA1.value = CD_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromStoA1.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromStoA1.value = q_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromStoA1.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromStoA1.value = WBL_fromStoA1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromStoA1.Attributes.unit = "daN";
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
        alfa_A = alfa_func(rho0, S, VA, WS, nA, CLalfa, alpha_zerol);
        CL_A   = CL_calc(obj1, nA, Mass, g, VA, rho0, S);
            if CL_A < CL_star - 0.03
                CL_A = CL_calc(obj1, nA, Mass, g, VA, rho0, S);
            elseif CL_A > CL_star + 0.03
                CL_A = CLmax_non_lin(alfa_A);
            end
        CD_A   = cd_calc(obj1, CD0, CL_A, AR, e, k1, k2);
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
        % NEW LIFT COEFFICIENT
        CL_A_new = CL_A - CLHT_A;
        % NEW LIFT COEFFICIENT
        LW_A_new = q_A*S*CL_A_new*1e-1;
        
        % STORE ALL THE DATA INSIDE THE STRUCT VARIABLE
        % POINT A
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value = q_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.Attributes.unit = "Pa";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA_new.value = LW_A_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.value = LHT_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.value = alfa_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA_rad.value = deg2rad(alfa_A);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CLHT_A.value = CLHT_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CLHT_A.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CD_A.value = CD_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CD_A.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CM_A.value = CMCG_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CM_A.Attributes.unit = "Non dimensional";          
        % =================================================================  
        % FROM A1 TO C
        CL_fromA1toC   = zeros(length(V_fromA1toC), 1);
        alfa_fromA1toC = zeros(length(V_fromA1toC), 1);
        CD_fromA1toC   = zeros(length(V_fromA1toC), 1);
        q_fromA1toC    = zeros(length(V_fromA1toC), 1);
        WBL_fromA1toC  = zeros(length(V_fromA1toC), 1);   
        CMCL_fromA1toC = zeros(length(V_fromA1toC), 1);
        CMCD_fromA1toC = zeros(length(V_fromA1toC), 1);
        CMCT_fromA1toC = zeros(length(V_fromA1toC), 1);
        CMCG_fromA1toC = zeros(length(V_fromA1toC), 1);    
        CLHT_fromA1toC = zeros(length(V_fromA1toC), 1); 
        LHT_fromA1toC  = zeros(length(V_fromA1toC), 1);
        for i = 1:length(V_fromA1toC)
            alfa_fromA1toC(i) = alfa_func(rho0, S, V_fromA1toC(i), WS, n_fromA1toC(i), CLalfa, alpha_zerol);
            CL_fromA1toC(i)   = CL_calc(obj1, n_fromA1toC(i), Mass, g, V_fromA1toC(i), rho0, S);
            if CL_fromA1toC(i) < CL_star - 0.03
                CL_fromA1toC(i) = CL_calc(obj1, n_fromA1toC(i), Mass, g, V_fromA1toC(i), rho0, S);
            elseif CL_fromA1toC(i) > CL_star + 0.03
                CL_fromA1toC(i) = CLmax_non_lin(alfa_fromA1toC(i));
            end
            CD_fromA1toC(i)   = cd_calc(obj1, CD0, CL_fromA1toC(i), AR, e, k1, k2);
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
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC.value = CL_fromA1toC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC.value = alfa_fromA1toC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC.Attributes.unit = "degrees"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC_rad.value = deg2rad(alfa_fromA1toC);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA1toC_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromA1toC.value = CD_fromA1toC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromA1toC.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromA1toC.value = q_fromA1toC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromA1toC.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromA1toC.value = WBL_fromA1toC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromA1toC.Attributes.unit = "daN";    
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
        CL_fromCtoA2   = zeros(length(V_fromCtoA2), 1);
        alfa_fromCtoA2 = zeros(length(V_fromCtoA2), 1);
        CD_fromCtoA2   = zeros(length(V_fromCtoA2), 1);
        q_fromCtoA2    = zeros(length(V_fromCtoA2), 1);
        WBL_fromCtoA2  = zeros(length(V_fromCtoA2), 1);   
        CMCL_fromCtoA2 = zeros(length(V_fromCtoA2), 1);
        CMCD_fromCtoA2 = zeros(length(V_fromCtoA2), 1);
        CMCT_fromCtoA2 = zeros(length(V_fromCtoA2), 1);
        CMCG_fromCtoA2 = zeros(length(V_fromCtoA2), 1);    
        CLHT_fromCtoA2 = zeros(length(V_fromCtoA2), 1); 
        LHT_fromCtoA2  = zeros(length(V_fromCtoA2), 1);
        for i = 1:length(V_fromCtoA2)
            alfa_fromCtoA2(i) = alfa_func(rho0, S, V_fromCtoA2(i), WS, n_fromCtoA2(i), CLalfa, alpha_zerol);
            CL_fromCtoA2(i)   = CL_calc(obj1, n_fromCtoA2(i), Mass, g, V_fromCtoA2(i), rho0, S);
            if CL_fromCtoA2(i) < CL_star
                CL_fromCtoA2(i) = CL_calc(obj1, n_fromCtoA2(i), Mass, g, V_fromCtoA2(i), rho0, S);
            elseif CL_fromCtoA2(i) > CL_star
                CL_fromCtoA2(i) = CLmax_non_lin(alfa_fromCtoA2(i));
            end
            CD_fromCtoA2(i)   = cd_calc(obj1, CD0, CL_fromCtoA2(i), AR, e, k1, k2);
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
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoA2.value = CL_fromCtoA2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoA2.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoA2.value = alfa_fromCtoA2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoA2.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoA2_rad.value = deg2rad(alfa_fromCtoA2);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromCtoA2_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtoA2.value = CD_fromCtoA2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtoA2.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromCtoA2.value = q_fromCtoA2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromCtoA2.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromCtoA2.value = WBL_fromCtoA2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromCtoA2.Attributes.unit = "daN";    
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
        CL_fromA2toD   = zeros(length(V_fromA2toD), 1);
        alfa_fromA2toD = zeros(length(V_fromA2toD), 1);
        CD_fromA2toD   = zeros(length(V_fromA2toD), 1);
        q_fromA2toD    = zeros(length(V_fromA2toD), 1);
        WBL_fromA2toD  = zeros(length(V_fromA2toD), 1);   
        CMCL_fromA2toD = zeros(length(V_fromA2toD), 1);
        CMCD_fromA2toD = zeros(length(V_fromA2toD), 1);
        CMCT_fromA2toD = zeros(length(V_fromA2toD), 1);
        CMCG_fromA2toD = zeros(length(V_fromA2toD), 1);      
        CLHT_fromA2toD = zeros(length(V_fromA2toD), 1);    
        LHT_fromA2toD  = zeros(length(V_fromA2toD), 1);    
        for i = 1:length(V_fromA2toD)
            alfa_fromA2toD(i) = alfa_func(rho0, S, V_fromA2toD(i), WS, n_fromA2toD(i), CLalfa, alpha_zerol);
            CL_fromA2toD(i)   = CL_calc(obj1, n_fromA2toD(i), Mass, g, V_fromA2toD(i), rho0, S);
            if CL_fromA2toD(i) < CL_star
                CL_fromA2toD(i) = CL_calc(obj1, n_fromA2toD(i), Mass, g, V_fromA2toD(i), rho0, S);
            elseif CL_fromA2toD(i) > CL_star
                CL_fromA2toD(i) = CLmax_non_lin(alfa_fromA2toD(i));
            end
            CD_fromA2toD(i) = cd_calc(obj1, CD0, CL_fromA2toD(i), AR, e, k1, k2);
            q_fromA2toD(i)   = 0.5*rho0*(V_fromA2toD(i))^2;
            WBL_fromA2toD(i) = q_fromA2toD(i)*S*CL_fromA2toD(i)*1e-1;   
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
        end        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA2toD.value = CL_fromA2toD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA2toD.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA2toD.value = alfa_fromA2toD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA2toD.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA2toD_rad.value = deg2rad(alfa_fromA2toD);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromA2toD_rad.Attributes.unit = "rad";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromA2toD.value = CD_fromA2toD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromA2toD.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromA2toD.value = q_fromA2toD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromA2toD.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromA2toD.value = WBL_fromA2toD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromA2toD.Attributes.unit = "daN"; 
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
        CL_fromDto0   = zeros(length(V_fromDto0), 1);
        alfa_fromDto0 = zeros(length(V_fromDto0), 1);
        CD_fromDto0   = zeros(length(V_fromDto0), 1);
        q_fromDto0    = zeros(length(V_fromDto0), 1);
        WBL_fromDto0  = zeros(length(V_fromDto0), 1);   
        CMCL_fromDto0 = zeros(length(V_fromDto0), 1);  
        CMCD_fromDto0 = zeros(length(V_fromDto0), 1);
        CMCT_fromDto0 = zeros(length(V_fromDto0), 1); 
        CMCG_fromDto0 = zeros(length(V_fromDto0), 1);        
        CLHT_fromDto0 = zeros(length(V_fromDto0), 1);      
        LHT_fromDto0  = zeros(length(V_fromDto0), 1);  
        for i = 1:length(V_fromDto0)
            alfa_fromDto0(i) = alfa_func(rho0, S, V_fromDto0(i), WS, n_fromDto0(i), CLalfa, alpha_zerol);
            CL_fromDto0(i)   = CL_calc(obj1, n_fromDto0(i), Mass, g, V_fromDto0(i), rho0, S);
            if CL_fromDto0(i) < CL_star
                CL_fromDto0(i) = CL_calc(obj1, n_fromDto0(i), Mass, g, V_fromDto0(i), rho0, S);
            elseif CL_fromA2toD(i) > CL_star
                CL_fromDto0(i) = CLmax_non_lin(alfa_fromDto0(i));
            end
            CD_fromDto0(i) = cd_calc(obj1, CD0, CL_fromDto0(i), AR, e, k1, k2);
            q_fromDto0(i)   = 0.5*rho0*(V_fromDto0(i))^2;
            WBL_fromDto0(i) = q_fromDto0(i)*S*CL_fromDto0(i)*1e-1;     
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
        end           
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0.value = CL_fromDto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0.value = alfa_fromDto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_rad.value = deg2rad(alfa_fromDto0);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromDto0_rad.Attributes.unit = "rad";         
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDto0.value = CD_fromDto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDto0.Attributes.unit = "Non dimensional";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromDto0.value = q_fromDto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromDto0.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromDto0.value = WBL_fromDto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBL_fromDto0.Attributes.unit = "daN";  
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
        end
        
        % WING BODY LIFT COEFFICIENT
        Aircraft_weight = ones(length(V_unit_load_factor), 1)*Mass*g;
        CLWB_unit_load_factor = zeros(length(V_unit_load_factor), 1);
        for i = 1:length(V_unit_load_factor)
            Aero_force = q_unit_load_factor(i)*S;
            CLWB_unit_load_factor(i) = ((Aircraft_weight(i))*(1/Aero_force));
        end
        
        % DRAG LIFT COEFFICIENT
        CD_unit_load_factor = zeros(length(V_unit_load_factor), 1);
        for i = 1:length(V_unit_load_factor)
            CD_unit_load_factor(i) = cd_calc(obj1, CD0, CLWB_unit_load_factor(i), ...
                                             AR, e, k1, k2);
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value = V_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.Attributes.unit = "m/s";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.value = q_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.Attributes.unit = "Pa";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value           = CLWB_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value = CD_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.Attributes.unit = "Non dimensional"; 
        
        % Interpolation coefficient from actual aerodynamic data
        CL_supp    = str2num(Aircraft.Certification.Aerodynamic_data.CL.value);
        alpha_supp = str2num(Aircraft.Certification.Aerodynamic_data.alpha.value);
        x          = 0.027*ones(length(CL_supp), 1);
        p = polyfit(CL_supp + x, alpha_supp, 2);
        CLalfa_deg = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;
        
        % ALPHA CALCULATION 
        alpha_unit_load_factor = zeros(length(V_unit_load_factor), 1);
        for i = 1:length(alpha_unit_load_factor)
            alpha_unit_load_factor(i) = alpha_calc(obj1, ...
                                                   CLWB_unit_load_factor(i), CL0, CL_star, CLalfa_deg, p);
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.value = alpha_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value = deg2rad(alpha_unit_load_factor);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.Attributes.unit = "radians"; 
        
        % WING BODY LIFT 
        WBLift_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
        for i = 1:length(alpha_unit_load_factor)   
            WBLift_unit_load_factor(i) = (0.5) * (V_unit_load_factor(i)^2) * rho0 * S * (CLWB_unit_load_factor(i))*(1E-1);   
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.value = WBLift_unit_load_factor;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.Attributes.unit = "daN";       
        
        % EVALUATION OF TAIL LOADS CONTRIBUTION     
        CMCL_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
        CMCD_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
        CMCT_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
        CMCG_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
        CLHT_unit_load_factor = zeros(length(alpha_unit_load_factor), 1);
        LHT_unit_load_factor  = zeros(length(alpha_unit_load_factor), 1);
        LW_unit_load_factor   = zeros(length(alpha_unit_load_factor), 1);
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
            LHT_unit_load_factor(i) = (0.5)*(V_unit_load_factor(i)^2)* S * rho0 * (CLHT_unit_load_factor(i))*(1e-1);
            LW_unit_load_factor(i)  = WBLift_unit_load_factor(i) - LHT_unit_load_factor(i);                                                                                     
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
        
        % TAIL LOADS DIAGRAM - POSITIVE SIDE   
        HT_balancing_loads = figure(5); 
        hold on; grid on; grid minor; 
        
%         ylim([LHT_from0toS(1) LHT_fromDto0(1)])
%         xlim([0.0 V_fromDto0(1)])
        
        plot(V_from0toS,  LHT_from0toS,  '-r', 'LineWidth', 1)
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
        xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
        ylabel("Horizontal tail lift - $L_{ht}$ (daN)", "Interpreter", "latex")
        title("Horizontal empennage airloads per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
        
        % MAIN WING LOADS DIAGRAM        
        CL_from0toS_new  = CL_from0toS  - CLHT_from0toS;        
        CL_fromStoA1_new = CL_fromStoA1 - CLHT_fromStoA1;        
        CL_fromA1toC_new = CL_fromA1toC - CLHT_fromA1toC;        
        CL_fromCtoA2_new = CL_fromCtoA2 - CLHT_fromCtoA2;        
        CL_fromA2toD_new = CL_fromA2toD - CLHT_fromA2toD;        
        CL_fromDto0_new  = CL_fromDto0  - CLHT_fromDto0;
        LW_from0toS_new  = zeros(length(CL_from0toS_new), 1);
        LW_fromStoA1_new = zeros(length(CL_from0toS_new), 1);
        LW_fromA1toC_new = zeros(length(CL_from0toS_new), 1);
        LW_fromCtoA2_new = zeros(length(CL_from0toS_new), 1);
        LW_fromA2toD_new = zeros(length(CL_from0toS_new), 1);
        LW_fromDto0_new  = zeros(length(CL_from0toS_new), 1);
        for i = 1:length(CL_from0toS_new)
            LW_from0toS_new(i)  = (0.5)*(V_from0toS(i)^2)* rho0 * S *(CL_from0toS_new(i))*(1e-1);  
            LW_fromStoA1_new(i) = (0.5)*(V_fromStoA1(i)^2)* rho0 * S *(CL_fromStoA1_new(i))*(1e-1);
            LW_fromA1toC_new(i) = (0.5)*(V_fromA1toC(i)^2)* rho0 * S *(CL_fromA1toC_new(i))*(1e-1);
            LW_fromCtoA2_new(i) = (0.5)*(V_fromCtoA2(i)^2)* rho0 * S *(CL_fromCtoA2_new(i))*(1e-1);
            LW_fromA2toD_new(i) = (0.5)*(V_fromA2toD(i)^2)* rho0 * S *(CL_fromA2toD_new(i))*(1e-1);
            LW_fromDto0_new(i)  = (0.5)*(V_fromDto0(i)^2)* rho0 * S *(CL_fromDto0_new(i))*(1e-1);
        end       
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS_new.value = CL_from0toS_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toS_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1_new.value = CL_fromStoA1_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromStoA1_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC_new.value = CL_fromA1toC_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA1toC_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoA2_new.value = CL_fromCtoA2_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoA2_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA2toD_new.value = CL_fromA2toD_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromA2toD_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.value = CL_fromDto0_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDto0_new.Attributes.unit = "Non dimensional";
        % ------------------------------------------------------------------------------------------------------------------------
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.value = LW_from0toS_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromStoA1_new.value = LW_fromStoA1_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromStoA1_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC_new.value = LW_fromA1toC_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA1toC_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoA2_new.value = LW_fromCtoA2_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromCtoA2_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA2toD_new.value = LW_fromA2toD_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromA2toD_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.value = LW_fromDto0_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromDto0_new.Attributes.unit = "daN";
        
        % TAIL LOADS DIAGRAM - POSITIVE SIDE   
        Wing_balancing_loads = figure(6); 
        hold on; grid on; grid minor; 
        
%         ylim([LHT_from0toS(1) LHT_fromDto0(1)])
%         xlim([0.0 V_fromDto0(1)])
        
        plot(V_from0toS,  LW_from0toS_new,  '-r', 'LineWidth', 1)
        plot(V_fromStoA1, LW_fromStoA1_new, '-r', 'LineWidth', 1)
        plot(V_fromA1toC, LW_fromA1toC_new, '-r', 'LineWidth', 1)
        plot(V_fromCtoA2, LW_fromCtoA2_new, '-r', 'LineWidth', 1)
        plot(V_fromA2toD, LW_fromA2toD_new, '-r', 'LineWidth', 1)
        plot(V_fromDto0,  LW_fromDto0_new,  '-r', 'LineWidth', 1)
        % UNIT LOAD DISTR.
        plot(V_unit_load_factor, LW_unit_load_factor,  '--k', 'LineWidth', 2)
        % ---------------------------------------------------------------------
        plot(V_from0toS(1),    LW_from0toS_new(1),    'k.', 'MarkerSize', 10)
        plot(V_from0toS(end),  LW_from0toS_new(end),  'k.', 'MarkerSize', 10)
        plot(VA,               LW_A_new,              'k.', 'MarkerSize', 10)
        plot(V_fromStoA1(end), LW_fromStoA1_new(end), 'k.', 'MarkerSize', 10)
        plot(V_fromA1toC(end), LW_fromA1toC_new(end), 'k.', 'MarkerSize', 10)
        plot(V_fromCtoA2(end), LW_fromCtoA2_new(end), 'k.', 'MarkerSize', 10)
        plot(V_fromA2toD(end), LW_fromA2toD_new(end), 'k.', 'MarkerSize', 10)
        plot(V_fromDto0(end),  LW_fromDto0_new(end),  'k.', 'MarkerSize', 10)
        % ---------------------------------------------------------------------
        text(V_fromStoA1(1),   LW_fromStoA1_new(1),   '  S',  'FontSize', 6)
        text(VA,               LW_A_new,              '  A',  'FontSize', 6)
        text(V_fromStoA1(end), LW_fromStoA1_new(end), '  A1', 'FontSize', 6)
        text(V_fromA1toC(end), LW_fromA1toC_new(end), '  C',  'FontSize', 6)
        text(V_fromCtoA2(end), LW_fromCtoA2_new(end), '  A2', 'FontSize', 6)
        text(V_fromA2toD(end), LW_fromA2toD_new(end), '  D',  'FontSize', 6)
        % text(40.75, -18, 'n = 1')
        % ---------------------------------------------------------------------
        xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
        ylabel("Full body lift - $L_{ht}$ (daN)", "Interpreter", "latex")
        title("Full body airloads per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct  
        
        % STORE INSIDE THE AIRCRAFT STRUCTURE VARIABLE

        % POINT S
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.value = q_from0toS(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.Attributes.unit = "Pa";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS_new.value = LW_from0toS_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTS.value = LHT_from0toS(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTS.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.value = alfa_from0toS(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_rad.value = deg2rad(alfa_from0toS(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alfaS_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CLHT_S.value = CLHT_from0toS(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CLHT_S.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CD_S.value = CD_from0toS(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CD_S.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CM_S.value = CMCG_from0toS(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CM_S.Attributes.unit = "Non dimensional";       
        % POINT A1
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.value = q_fromStoA1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.Attributes.unit = "Pa";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1_new.value = LW_fromStoA1_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.value = LHT_fromStoA1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.value = alfa_fromStoA1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_rad.value = deg2rad(alfa_fromStoA1(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CLHT_A1.value = CLHT_fromStoA1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CLHT_A1.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CD_A1.value = CD_fromStoA1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CD_A1.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CM_A1.value = CMCG_fromStoA1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CM_A1.Attributes.unit = "Non dimensional";
        % POINT C
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value = q_fromA1toC(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.Attributes.unit = "Pa"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.value = LW_fromA1toC_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.value = LHT_fromA1toC(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value = alfa_fromA1toC(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_rad.value = deg2rad(alfa_fromA1toC(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC_rad.Attributes.unit = "rad"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CLHT_C.value = CLHT_fromA1toC(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CLHT_C.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CD_C.value = CD_fromA1toC(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CD_C.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CM_C.value = CMCG_fromA1toC(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CM_C.Attributes.unit = "Non dimensional";        
        % POINT A2
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.qA2.value = q_fromCtoA2(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.qA2.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.LA2_new.value = LW_fromCtoA2_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.LA2_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.LHTA2.value = LHT_fromCtoA2(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.LHTA2.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.alfaA2.value = alfa_fromCtoA2(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.alfaA2.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.alfaA2_rad.value = deg2rad(alfa_fromCtoA2(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.alfaA2_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CLHT_A2.value = CLHT_fromCtoA2(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CLHT_A2.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CD_A2.value = CD_fromCtoA2(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CD_A2.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CM_A2.value = CMCG_fromCtoA2(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CM_A2.Attributes.unit = "Non dimensional";
        % POINT D
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value = q_fromA2toD(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.Attributes.unit = "Pa"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.value = LW_fromA2toD_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value = LHT_fromA2toD(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value = alfa_fromA2toD(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_rad.value = deg2rad(alfa_fromA2toD(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD_rad.Attributes.unit = "rad";  
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
            % =============================================================
            % FROM 0 TO S_INVERTED 
            CL_from0toSinv   = zeros(length(V_from0toSinv), 1);
            alfa_from0toSinv = zeros(length(V_from0toSinv), 1);
            CD_from0toSinv   = zeros(length(V_from0toSinv), 1);
            q_from0toSinv    = zeros(length(V_from0toSinv), 1);
            WBL_from0toSinv  = zeros(length(V_from0toSinv), 1);
            CMCL_from0toSinv = zeros(length(V_from0toSinv), 1);
            CMCD_from0toSinv = zeros(length(V_from0toSinv), 1);
            CMCT_from0toSinv = zeros(length(V_from0toSinv), 1);
            CMCG_from0toSinv = zeros(length(V_from0toSinv), 1); 
            CLHT_from0toSinv = zeros(length(V_from0toSinv), 1);
            LHT_from0toSinv  = zeros(length(V_from0toSinv), 1);
            for i = 1:length(V_from0toSinv)
                alfa_from0toSinv(i) = alfa_func(rho0, S, V_from0toSinv(i), WS, n_from0toSinv(i), CLalfa, alpha_zerol);
                CL_from0toSinv(i)   = CL_calc(obj1, n_from0toSinv(i), Mass, g, V_from0toSinv(i), rho0, S);
                if CL_from0toSinv(i) < CL_star
                    CL_from0toSinv(i) = CL_calc(obj1, n_from0toSinv(i), Mass, g, V_from0toSinv(i), rho0, S);
                elseif CL_from0toSinv(i) > CL_star
                    CL_from0toSinv(i) = CLmax_non_lin(alfa_from0toSinv(i));
                end
                CD_from0toSinv(i)   = cd_calc(obj1, CD0, CL_from0toSinv(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv.value = CL_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv.value = alfa_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_rad.value = deg2rad(alfa_from0toSinv);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toSinv.value = CD_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toSinv.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toSinv.value = q_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toSinv.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toSinv.value = WBL_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toSinv.Attributes.unit = "daN";
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
            CMCL_fromSinvtoG = zeros(length(V_fromSinvtoG), 1);
            CMCD_fromSinvtoG = zeros(length(V_fromSinvtoG), 1);
            CMCT_fromSinvtoG = zeros(length(V_fromSinvtoG), 1);
            CMCG_fromSinvtoG = zeros(length(V_fromSinvtoG), 1); 
            CLHT_fromSinvtoG = zeros(length(V_fromSinvtoG), 1);
            LHT_fromSinvtoG  = zeros(length(V_fromSinvtoG), 1);
            for i = 1:length(V_fromSinvtoG)
                alfa_fromSinvtoG(i) = alfa_func(rho0, S, V_fromSinvtoG(i), WS, n_fromSinvtoG(i), CLalfa, alpha_zerol);
                CL_fromSinvtoG(i)   = CL_calc(obj1, n_fromSinvtoG(i), Mass, g, V_fromSinvtoG(i), rho0, S);
                if CL_fromSinvtoG(i) < CL_star
                    CL_fromSinvtoG(i) = CL_calc(obj1, n_fromSinvtoG(i), Mass, g, V_fromSinvtoG(i), rho0, S);
                elseif CL_fromSinvtoG(i) > CL_star
                    CL_fromSinvtoG(i) = CLmax_non_lin(alfa_fromSinvtoG(i));
                end
                CD_fromSinvtoG(i)   = cd_calc(obj1, CD0, CL_fromSinvtoG(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG.value = CL_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG.value = alfa_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG_rad.value = deg2rad(alfa_fromSinvtoG);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromSinvtoG.value = CD_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromSinvtoG.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromSinvtoG.value = q_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromSinvtoG.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromSinvtoG.value = WBL_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromSinvtoG.Attributes.unit = "daN";
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
            CL_fromGtoG1   = zeros(length(V_fromGtoG1), 1);
            alfa_fromGtoG1 = zeros(length(V_fromGtoG1), 1);
            CD_fromGtoG1   = zeros(length(V_fromGtoG1), 1);
            q_fromGtoG1    = zeros(length(V_fromGtoG1), 1);
            WBL_fromGtoG1  = zeros(length(V_fromGtoG1), 1);
            CMCL_fromGtoG1 = zeros(length(V_fromGtoG1), 1);
            CMCD_fromGtoG1 = zeros(length(V_fromGtoG1), 1);
            CMCT_fromGtoG1 = zeros(length(V_fromGtoG1), 1);
            CMCG_fromGtoG1 = zeros(length(V_fromGtoG1), 1); 
            CLHT_fromGtoG1 = zeros(length(V_fromGtoG1), 1);
            LHT_fromGtoG1  = zeros(length(V_fromGtoG1), 1);
            for i = 1:length(V_fromGtoG1)
                alfa_fromGtoG1(i) = alfa_func(rho0, S, V_fromGtoG1(i), WS, n_fromGtoG1(i), CLalfa, alpha_zerol);
                CL_fromGtoG1(i)   = CL_calc(obj1, n_fromGtoG1(i), Mass, g, V_fromGtoG1(i), rho0, S);
                if CL_fromGtoG1(i) < CL_star
                    CL_fromGtoG1(i) = CL_calc(obj1, n_fromGtoG1(i), Mass, g, V_fromGtoG1(i), rho0, S);
                elseif CL_fromGtoG1(i) > CL_star
                    CL_fromGtoG1(i) = CLmax_non_lin(alfa_fromGtoG1(i));
                end
                CD_fromGtoG1(i)   = cd_calc(obj1, CD0, CL_fromGtoG1(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoG1.value = CL_fromGtoG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoG1.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoG1.value = alfa_fromGtoG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoG1.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoG1_rad.value = deg2rad(alfa_fromGtoG1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoG1_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromGtoG1.value = CD_fromGtoG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromGtoG1.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromGtoG1.value = q_fromGtoG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromGtoG1.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromGtoG1.value = WBL_fromGtoG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromGtoG1.Attributes.unit = "daN";
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
            CMCL_fromG1toF = zeros(length(V_fromG1toF), 1);
            CMCD_fromG1toF = zeros(length(V_fromG1toF), 1);
            CMCT_fromG1toF = zeros(length(V_fromG1toF), 1);
            CMCG_fromG1toF = zeros(length(V_fromG1toF), 1); 
            CLHT_fromG1toF = zeros(length(V_fromG1toF), 1);
            LHT_fromG1toF  = zeros(length(V_fromG1toF), 1);
            for i = 1:length(V_fromG1toF)
                alfa_fromG1toF(i) = alfa_func(rho0, S, V_fromG1toF(i), WS, n_fromG1toF(i), CLalfa, alpha_zerol);
                CL_fromG1toF(i)   = CL_calc(obj1, n_fromG1toF(i), Mass, g, V_fromG1toF(i), rho0, S);
                if CL_fromG1toF(i) < CL_star
                    CL_fromG1toF(i) = CL_calc(obj1, n_fromG1toF(i), Mass, g, V_fromG1toF(i), rho0, S);
                elseif CL_fromG1toF(i) > CL_star
                    CL_fromG1toF(i) = CLmax_non_lin(alfa_fromG1toF(i));
                end
                CD_fromG1toF(i)   = cd_calc(obj1, CD0, CL_fromG1toF(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF.value = CL_fromG1toF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF.value = alfa_fromG1toF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF_rad.value = deg2rad(alfa_fromG1toF);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromG1toF.value = CD_fromG1toF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromG1toF.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromG1toF.value = q_fromG1toF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromG1toF.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromG1toF.value = WBL_fromG1toF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromG1toF.Attributes.unit = "daN";
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
            CMCL_fromFtoG2 = zeros(length(V_fromFtoG2), 1);
            CMCD_fromFtoG2 = zeros(length(V_fromFtoG2), 1);
            CMCT_fromFtoG2 = zeros(length(V_fromFtoG2), 1);
            CMCG_fromFtoG2 = zeros(length(V_fromFtoG2), 1); 
            CLHT_fromFtoG2 = zeros(length(V_fromFtoG2), 1);
            LHT_fromFtoG2  = zeros(length(V_fromFtoG2), 1);
            for i = 1:length(V_fromFtoG2)
                alfa_fromFtoG2(i) = alfa_func(rho0, S, V_fromFtoG2(i), WS, n_fromFtoG2(i), CLalfa, alpha_zerol);
                CL_fromFtoG2(i)   = CL_calc(obj1, n_fromFtoG2(i), Mass, g, V_fromFtoG2(i), rho0, S);
                if CL_fromFtoG2(i) < CL_star
                    CL_fromFtoG2(i) = CL_calc(obj1, n_fromFtoG2(i), Mass, g, V_fromFtoG2(i), rho0, S);
                elseif CL_fromFtoG2(i) > CL_star
                    CL_fromFtoG2(i) = CLmax_non_lin(alfa_fromFtoG2(i));
                end
                CD_fromFtoG2(i)   = cd_calc(obj1, CD0, CL_fromFtoG2(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoG2.value = CL_fromFtoG2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoG2.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoG2.value = alfa_fromFtoG2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoG2.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoG2_rad.value = deg2rad(alfa_fromFtoG2);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoG2_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoG2.value = CD_fromFtoG2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoG2.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromFtoG2.value = q_fromFtoG2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromFtoG2.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoG2.value = WBL_fromFtoG2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoG2.Attributes.unit = "daN";
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
            CMCL_fromG2toE = zeros(length(V_fromG2toE), 1);
            CMCD_fromG2toE = zeros(length(V_fromG2toE), 1);
            CMCT_fromG2toE = zeros(length(V_fromG2toE), 1);
            CMCG_fromG2toE = zeros(length(V_fromG2toE), 1); 
            CLHT_fromG2toE = zeros(length(V_fromG2toE), 1);
            LHT_fromG2toE  = zeros(length(V_fromG2toE), 1);
            for i = 1:length(V_fromG2toE)
                alfa_fromG2toE(i) = alfa_func(rho0, S, V_fromG2toE(i), WS, n_fromG2toE(i), CLalfa, alpha_zerol);
                CL_fromG2toE(i)   = CL_calc(obj1, n_fromG2toE(i), Mass, g, V_fromG2toE(i), rho0, S);
                if CL_fromG2toE(i) < CL_star
                    CL_fromG2toE(i) = CL_calc(obj1, n_fromG2toE(i), Mass, g, V_fromG2toE(i), rho0, S);
                elseif CL_fromG2toE(i) > CL_star
                    CL_fromG2toE(i) = CLmax_non_lin(alfa_fromG2toE(i));
                end
                CD_fromG2toE(i)   = cd_calc(obj1, CD0, CL_fromG2toE(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG2toE.value = CL_fromG2toE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG2toE.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG2toE.value = alfa_fromG2toE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG2toE.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG2toE_rad.value = deg2rad(alfa_fromG2toE);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG2toE_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromG2toE.value = CD_fromG2toE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromG2toE.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromG2toE.value = q_fromG2toE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromG2toE.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromG2toE.value = WBL_fromG2toE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromG2toE.Attributes.unit = "daN";
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
            CMCL_fromEto0 = zeros(length(V_fromEto0), 1);
            CMCD_fromEto0 = zeros(length(V_fromEto0), 1);
            CMCT_fromEto0 = zeros(length(V_fromEto0), 1);
            CMCG_fromEto0 = zeros(length(V_fromEto0), 1); 
            CLHT_fromEto0 = zeros(length(V_fromEto0), 1);
            LHT_fromEto0  = zeros(length(V_fromEto0), 1);
            for i = 1:length(V_fromEto0)
                alfa_fromEto0(i) = alfa_func(rho0, S, V_fromEto0(i), WS, n_fromEto0(i), CLalfa, alpha_zerol);
                CL_fromEto0(i)   = CL_calc(obj1, n_fromEto0(i), Mass, g, V_fromEto0(i), rho0, S);
                if CL_fromEto0(i) < CL_star
                    CL_fromEto0(i) = CL_calc(obj1, n_fromEto0(i), Mass, g, V_fromEto0(i), rho0, S);
                elseif CL_fromEto0(i) > CL_star
                    CL_fromEto0(i) = CLmax_non_lin(alfa_fromEto0(i));
                end
                CD_fromEto0(i)   = cd_calc(obj1, CD0, CL_fromEto0(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0.value = CL_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0.value = alfa_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_rad.value = deg2rad(alfa_fromEto0);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromEto0.value = CD_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromEto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromEto0.value = q_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromEto0.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromEto0.value = WBL_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromEto0.Attributes.unit = "daN";
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
            % MAIN WING LOADS DIAGRAM        
            CL_from0toSinv_new = CL_from0toSinv  - CLHT_from0toSinv;        
            CL_fromSinvtoG_new = CL_fromSinvtoG - CLHT_fromSinvtoG;        
            CL_fromGtoG1_new   = CL_fromGtoG1    - CLHT_fromGtoG1;        
            CL_fromG1toF_new   = CL_fromG1toF    - CLHT_fromG1toF;        
            CL_fromFtoG2_new   = CL_fromFtoG2    - CLHT_fromFtoG2;        
            CL_fromG2toE_new   = CL_fromG2toE    - CLHT_fromG2toE;        
            CL_fromEto0_new    = CL_fromEto0     - CLHT_fromEto0;
            LW_from0toSinv_new = zeros(length(CL_from0toSinv_new), 1);
            LW_fromSinvtoG_new = zeros(length(CL_from0toSinv_new), 1);
            LW_fromGtoG1_new   = zeros(length(CL_from0toSinv_new), 1);
            LW_fromG1toF_new   = zeros(length(CL_from0toSinv_new), 1);
            LW_fromFtoG2_new   = zeros(length(CL_from0toSinv_new), 1);
            LW_fromG2toE_new   = zeros(length(CL_from0toSinv_new), 1);
            LW_fromEto0_new    = zeros(length(CL_from0toSinv_new), 1);
            for i = 1:length(CL_from0toSinv_new)
                LW_from0toSinv_new(i) = (0.5)*(V_from0toSinv(i)^2)* rho0 * S *(CL_from0toSinv_new(i))*(1e-1);  
                LW_fromSinvtoG_new(i) = (0.5)*(V_fromSinvtoG(i)^2)* rho0 * S *(CL_fromSinvtoG_new(i))*(1e-1);
                LW_fromGtoG1_new(i)   = (0.5)*(V_fromGtoG1(i)^2)* rho0 * S *(CL_fromGtoG1_new(i))*(1e-1);
                LW_fromG1toF_new(i)   = (0.5)*(V_fromG1toF(i)^2)* rho0 * S *(CL_fromG1toF_new(i))*(1e-1);
                LW_fromFtoG2_new(i)   = (0.5)*(V_fromFtoG2(i)^2)* rho0 * S *(CL_fromFtoG2_new(i))*(1e-1);
                LW_fromG2toE_new(i)   = (0.5)*(V_fromG2toE(i)^2)* rho0 * S *(CL_fromG2toE_new(i))*(1e-1);
                LW_fromEto0_new(i)    = (0.5)*(V_fromEto0(i)^2)* rho0 * S *(CL_fromEto0_new(i))*(1e-1);
            end       
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.value = CL_from0toSinv_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG_new.value = CL_fromSinvtoG_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoG1_new.value = CL_fromGtoG1_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoG1_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF_new.value = CL_fromG1toF_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoG2_new.value = CL_fromFtoG2_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoG2_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG2toE_new.value = CL_fromG2toE_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG2toE_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.value = CL_fromEto0_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.Attributes.unit = "Non dimensional";
            % ------------------------------------------------------------------------------------------------------------------------
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.value = LW_from0toSinv_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromSinvtoG_new.value = LW_fromSinvtoG_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromSinvtoG_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromGtoG1_new.value = LW_fromGtoG1_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromGtoG1_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG1toF_new.value = LW_fromG1toF_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG1toF_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoG2_new.value = LW_fromFtoG2_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoG2_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG2toE_new.value = LW_fromG2toE_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG2toE_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.value = LW_fromEto0_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.Attributes.unit = "daN"; 
            % =============================================================
            % WING LOADS DIAGRAM - NEGATIVE SIDE   
            figure(6)
            % ---------------------------------------------------------------------
            plot(V_from0toSinv, LW_from0toSinv_new, '-r', 'LineWidth', 1)
            plot(V_fromSinvtoG, LW_fromSinvtoG_new, '-r', 'LineWidth', 1)
            plot(V_fromGtoG1,   LW_fromGtoG1_new,   '-r', 'LineWidth', 1)
            plot(V_fromG1toF,   LW_fromG1toF_new,   '-r', 'LineWidth', 1)
            plot(V_fromFtoG2,   LW_fromFtoG2_new,   '-r', 'LineWidth', 1)
            plot(V_fromG2toE,   LW_fromG2toE_new,   '-r', 'LineWidth', 1)
            plot(V_fromEto0,    LW_fromEto0_new,    '-r', 'LineWidth', 1)
            % ---------------------------------------------------------------------
            plot(V_from0toSinv(1),   LW_from0toSinv_new(1),   'k.', 'MarkerSize', 10)
            plot(V_from0toSinv(end), LW_from0toSinv_new(end), 'k.', 'MarkerSize', 10)
            plot(V_fromSinvtoG(end), LW_fromSinvtoG_new(end), 'k.', 'MarkerSize', 10)
            plot(V_fromGtoG1(end),   LW_fromGtoG1_new(end),   'k.', 'MarkerSize', 10)
            plot(V_fromG1toF(end),   LW_fromG1toF_new(end),   'k.', 'MarkerSize', 10)
            plot(V_fromFtoG2(end),   LW_fromFtoG2_new(end),   'k.', 'MarkerSize', 10)
            plot(V_fromG2toE(end),   LW_fromG2toE_new(end),   'k.', 'MarkerSize', 10)
            % ---------------------------------------------------------------------
            text(V_from0toSinv(end),  LW_from0toSinv_new(end), ' S inv.', 'FontSize', 6)
            text(V_fromSinvtoG(end),  LW_fromSinvtoG_new(end), ' G',      'FontSize', 6)
            text(V_fromGtoG1(end),    LW_fromGtoG1_new(end),   ' G1',     'FontSize', 6)
            text(V_fromG1toF(end),    LW_fromG1toF_new(end),   ' F',      'FontSize', 6) 
            text(V_fromFtoG2(end),    LW_fromFtoG2_new(end),   ' G2',     'FontSize', 6)  
            text(V_fromG2toE(end),    LW_fromG2toE_new(end),   ' E',      'FontSize', 6) 
            
            % STORE INSIDE THE AIRCRAFT STRUCTURE VARIABLE

            % POINT S_INV
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.value = q_from0toSinv(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.Attributes.unit = "Pa";               
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv_new.value = LW_from0toSinv_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LHTS_inv.value = LHT_from0toSinv(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LHTS_inv.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv.value = alfa_from0toSinv(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_rad.value = deg2rad(alfa_from0toSinv(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_rad.Attributes.unit = "rad";  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CLHT_S_inv.value = CLHT_from0toSinv(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CLHT_S_inv.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CD_S_inv.value = CD_from0toSinv(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CD_S_inv.Attributes.unit = "Non dimensional";  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CM_S_inv.value = CMCG_from0toSinv(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CM_S_inv.Attributes.unit = "Non dimensional";    
            % POINT G
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value = q_fromSinvtoG(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.Attributes.unit = "Pa";             
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG_new.value = LW_fromSinvtoG_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.value = LHT_fromSinvtoG(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.value = alfa_fromSinvtoG(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_rad.value = deg2rad(alfa_fromSinvtoG(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_rad.Attributes.unit = "rad";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CLHT_G.value = CLHT_fromSinvtoG(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CLHT_G.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.value = CD_fromSinvtoG(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CM_G.value = CMCG_fromSinvtoG(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CM_G.Attributes.unit = "Non dimensional";
            % POINT G1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.qG1.value = q_fromGtoG1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.qG1.Attributes.unit = "Pa";            
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1_new.value = LW_fromGtoG1_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LHTG1.value = LHT_fromGtoG1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LHTG1.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1.value = alfa_fromGtoG1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_rad.value = deg2rad(alfa_fromGtoG1(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_rad.Attributes.unit = "rad";  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CLHT_G1.value = CLHT_fromGtoG1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CLHT_G1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CD_G1.value = CD_fromGtoG1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CD_G1.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CM_G1.value = CMCG_fromGtoG1(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CM_G1.Attributes.unit = "Non dimensional";         
            % POINT F
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value = q_fromG1toF(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.Attributes.unit = "Pa";             
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF_new.value = LW_fromG1toF_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.value = LHT_fromG1toF(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.value = alfa_fromG1toF(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_rad.value = deg2rad(alfa_fromG1toF(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_rad.Attributes.unit = "rad";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CLHT_F.value = CLHT_fromG1toF(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CLHT_F.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.value = CD_fromG1toF(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CM_F.value = CMCG_fromG1toF(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CM_F.Attributes.unit = "Non dimensional";
            % POINT G2
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.qG2.value = q_fromFtoG2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.qG2.Attributes.unit = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.LG2_new.value = LW_fromFtoG2_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.LG2_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.LHTG2.value = LHT_fromFtoG2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.LHTG2.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.alfaG2.value = alfa_fromFtoG2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.alfaG2.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.alfaG2_rad.value = deg2rad(alfa_fromFtoG2(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.alfaG2_rad.Attributes.unit = "rad"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CLHT_G2.value = CLHT_fromFtoG2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CLHT_G2.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CD_G2.value = CD_fromFtoG2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CD_G2.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CM_G2.value = CMCG_fromFtoG2(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CM_G2.Attributes.unit = "Non dimensional";
            % POINT E
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value = q_fromG2toE(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.Attributes.unit = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE_new.value = LW_fromG2toE_new(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.value = LHT_fromG2toE(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.Attributes.unit = "daN";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.value = alfa_fromG2toE(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_rad.value = deg2rad(alfa_fromG2toE(end));
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_rad.Attributes.unit = "rad";  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CLHT_E.value = CLHT_fromG2toE(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CLHT_E.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CD_E.value = CD_fromG2toE(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CD_E.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CM_E.value = CMCG_fromG2toE(end);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CM_E.Attributes.unit = "Non dimensional";           
            
        elseif abs(min(n_gust_cruise_neg)) < abs(nmin)             
            % =============================================================
            % FROM 0 TO S_INVERTED 
            CL_from0toSinv   = zeros(length(V_from0toSinv), 1);
            alfa_from0toSinv = zeros(length(V_from0toSinv), 1);
            CD_from0toSinv   = zeros(length(V_from0toSinv), 1);
            q_from0toSinv    = zeros(length(V_from0toSinv), 1);
            WBL_from0toSinv  = zeros(length(V_from0toSinv), 1);
            CMCL_from0toSinv = zeros(length(V_from0toSinv), 1);
            CMCD_from0toSinv = zeros(length(V_from0toSinv), 1);
            CMCT_from0toSinv = zeros(length(V_from0toSinv), 1);
            CMCG_from0toSinv = zeros(length(V_from0toSinv), 1); 
            CLHT_from0toSinv = zeros(length(V_from0toSinv), 1);
            LHT_from0toSinv  = zeros(length(V_from0toSinv), 1);
            for i = 1:length(V_from0toSinv)
                alfa_from0toSinv(i) = alfa_func(rho0, S, V_from0toSinv(i), WS, n_from0toSinv(i), CLalfa, alpha_zerol);
                CL_from0toSinv(i)   = CL_calc(obj1, n_from0toSinv(i), Mass, g, V_from0toSinv(i), rho0, S);
                if CL_from0toSinv(i) < CL_star
                    CL_from0toSinv(i) = CL_calc(obj1, n_from0toSinv(i), Mass, g, V_from0toSinv(i), rho0, S);
                elseif CL_from0toSinv(i) > CL_star
                    CL_from0toSinv(i) = CLmax_non_lin(alfa_from0toSinv(i));
                end
                CD_from0toSinv(i)   = cd_calc(obj1, CD0, CL_from0toSinv(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv.value = CL_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv.value = alfa_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_rad.value = deg2rad(alfa_from0toSinv);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toSinv.value = CD_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toSinv.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toSinv.value = q_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toSinv.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toSinv.value = WBL_from0toSinv;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toSinv.Attributes.unit = "daN";
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
            CMCL_fromSinvtoG = zeros(length(V_fromSinvtoG), 1);
            CMCD_fromSinvtoG = zeros(length(V_fromSinvtoG), 1);
            CMCT_fromSinvtoG = zeros(length(V_fromSinvtoG), 1);
            CMCG_fromSinvtoG = zeros(length(V_fromSinvtoG), 1); 
            CLHT_fromSinvtoG = zeros(length(V_fromSinvtoG), 1);
            LHT_fromSinvtoG  = zeros(length(V_fromSinvtoG), 1);
            for i = 1:length(V_fromSinvtoG)
                alfa_fromSinvtoG(i) = alfa_func(rho0, S, V_fromSinvtoG(i), WS, n_fromSinvtoG(i), CLalfa, alpha_zerol);
                CL_fromSinvtoG(i)   = CL_calc(obj1, n_fromSinvtoG(i), Mass, g, V_fromSinvtoG(i), rho0, S);
                if CL_fromSinvtoG(i) < CL_star
                    CL_fromSinvtoG(i) = CL_calc(obj1, n_fromSinvtoG(i), Mass, g, V_fromSinvtoG(i), rho0, S);
                elseif CL_fromSinvtoG(i) > CL_star
                    CL_fromSinvtoG(i) = CLmax_non_lin(alfa_fromSinvtoG(i));
                end
                CD_fromSinvtoG(i)   = cd_calc(obj1, CD0, CL_fromSinvtoG(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG.value = CL_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG.value = alfa_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG_rad.value = deg2rad(alfa_fromSinvtoG);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromSinvtoG.value = CD_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromSinvtoG.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromSinvtoG.value = q_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromSinvtoG.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromSinvtoG.value = WBL_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromSinvtoG.Attributes.unit = "daN";
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
            CMCL_fromGtoF = zeros(length(V_fromGtoF), 1);
            CMCD_fromGtoF = zeros(length(V_fromGtoF), 1);
            CMCT_fromGtoF = zeros(length(V_fromGtoF), 1);
            CMCG_fromGtoF = zeros(length(V_fromGtoF), 1); 
            CLHT_fromGtoF = zeros(length(V_fromGtoF), 1);
            LHT_fromGtoF  = zeros(length(V_fromGtoF), 1);
            for i = 1:length(V_fromGtoF)
                alfa_fromGtoF(i) = alfa_func(rho0, S, V_fromGtoF(i), WS, n_fromGtoF(i), CLalfa, alpha_zerol);
                CL_fromGtoF(i)   = CL_calc(obj1, n_fromGtoF(i), Mass, g, V_fromGtoF(i), rho0, S);
                if CL_fromGtoF(i) < CL_star
                    CL_fromGtoF(i) = CL_calc(obj1, n_fromGtoF(i), Mass, g, V_fromGtoF(i), rho0, S);
                elseif CL_fromGtoF(i) > CL_star
                    CL_fromGtoF(i) = CLmax_non_lin(alfa_fromGtoF(i));
                end
                CD_fromGtoF(i)   = cd_calc(obj1, CD0, CL_fromGtoF(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoF.value = CL_fromGtoF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoF.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoF.value = alfa_fromGtoF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoF.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoF_rad.value = deg2rad(alfa_fromGtoF);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromGtoF_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromGtoF.value = CD_fromGtoF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromGtoF.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromGtoF.value = q_fromGtoF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromGtoF.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromGtoF.value = WBL_fromGtoF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromGtoF.Attributes.unit = "daN";
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
            CMCL_fromFtoE = zeros(length(V_fromFtoE), 1);
            CMCD_fromFtoE = zeros(length(V_fromFtoE), 1);
            CMCT_fromFtoE = zeros(length(V_fromFtoE), 1);
            CMCG_fromFtoE = zeros(length(V_fromFtoE), 1); 
            CLHT_fromFtoE = zeros(length(V_fromFtoE), 1);
            LHT_fromFtoE  = zeros(length(V_fromFtoE), 1);
            for i = 1:length(V_fromFtoE)
                alfa_fromFtoE(i) = alfa_func(rho0, S, V_fromFtoE(i), WS, n_fromFtoE(i), CLalfa, alpha_zerol);
                CL_fromFtoE(i)   = CL_calc(obj1, n_fromFtoE(i), Mass, g, V_fromFtoE(i), rho0, S);
                if CL_fromFtoE(i) < CL_star
                    CL_fromFtoE(i) = CL_calc(obj1, n_fromFtoE(i), Mass, g, V_fromFtoE(i), rho0, S);
                elseif CL_fromFtoE(i) > CL_star
                    CL_fromFtoE(i) = CLmax_non_lin(alfa_fromFtoE(i));
                end
                CD_fromFtoE(i)   = cd_calc(obj1, CD0, CL_fromFtoE(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value = CL_fromFtoE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE.value = alfa_fromFtoE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE_rad.value = deg2rad(alfa_fromFtoE);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value = CD_fromFtoE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromFtoE.value = q_fromFtoE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromFtoE.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE.value = WBL_fromFtoE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE.Attributes.unit = "daN";
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
            CMCL_fromEto0 = zeros(length(V_fromEto0), 1);
            CMCD_fromEto0 = zeros(length(V_fromEto0), 1);
            CMCT_fromEto0 = zeros(length(V_fromEto0), 1);
            CMCG_fromEto0 = zeros(length(V_fromEto0), 1); 
            CLHT_fromEto0 = zeros(length(V_fromEto0), 1);
            LHT_fromEto0  = zeros(length(V_fromEto0), 1);
            for i = 1:length(V_fromEto0)
                alfa_fromEto0(i) = alfa_func(rho0, S, V_fromEto0(i), WS, n_fromEto0(i), CLalfa, alpha_zerol);
                CL_fromEto0(i)   = CL_calc(obj1, n_fromEto0(i), Mass, g, V_fromEto0(i), rho0, S);
                if CL_fromEto0(i) < CL_star
                    CL_fromEto0(i) = CL_calc(obj1, n_fromEto0(i), Mass, g, V_fromEto0(i), rho0, S);
                elseif CL_fromEto0(i) > CL_star
                    CL_fromEto0(i) = CLmax_non_lin(alfa_fromEto0(i));
                end
                CD_fromEto0(i)   = cd_calc(obj1, CD0, CL_fromEto0(i), AR, e, k1, k2);
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
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0.value = CL_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0.value = alfa_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0.Attributes.unit = "degrees";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_rad.value = deg2rad(alfa_fromEto0);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_rad.Attributes.unit = "rad";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromEto0.value = CD_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromEto0.Attributes.unit = "Non dimensional"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromEto0.value = q_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromEto0.Attributes.unit = "Pa";        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromEto0.value = WBL_fromEto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromEto0.Attributes.unit = "daN";
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
            % MAIN WING LOADS DIAGRAM        
            CL_from0toSinv_new = CL_from0toSinv - CLHT_from0toSinv;        
            CL_fromSinvtoG_new = CL_fromSinvtoG - CLHT_fromSinvtoG;        
            CL_fromGtoF_new    = CL_fromGtoF    - CLHT_fromGtoF;        
            CL_fromFtoE_new    = CL_fromFtoE    - CLHT_fromFtoE;       
            CL_fromEto0_new    = CL_fromEto0    - CLHT_fromEto0;
            LW_from0toSinv_new = zeros(length(CL_from0toSinv_new), 1);
            LW_fromSinvtoG_new = zeros(length(CL_from0toSinv_new), 1);
            LW_fromGtoF_new    = zeros(length(CL_from0toSinv_new), 1);
            LW_fromFtoE_new    = zeros(length(CL_from0toSinv_new), 1);
            LW_fromEto0_new    = zeros(length(CL_from0toSinv_new), 1);
            for i = 1:length(CL_from0toSinv_new)
                LW_from0toSinv_new(i)  = (0.5)*(V_from0toSinv(i)^2)* rho0 * S *(CL_from0toSinv_new(i))*(1e-1);  
                LW_fromSinvtoG_new(i) = (0.5)*(V_fromSinvtoG(i)^2)* rho0 * S *(CL_fromSinvtoG_new(i))*(1e-1);
                LW_fromGtoF_new(i)    = (0.5)*(V_fromGtoF(i)^2)* rho0 * S *(CL_fromGtoF_new(i))*(1e-1);
                LW_fromFtoE_new(i)     = (0.5)*(V_fromFtoE(i)^2)* rho0 * S *(CL_fromFtoE_new(i))*(1e-1);
                LW_fromEto0_new(i)     = (0.5)*(V_fromEto0(i)^2)* rho0 * S *(CL_fromEto0_new(i))*(1e-1);
            end       
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.value = CL_from0toSinv_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG_new.value = CL_fromSinvtoG_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoF_new.value = CL_fromGtoF_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromGtoF_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value = CL_fromFtoE_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.value = CL_fromEto0_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.Attributes.unit = "Non dimensional";
            % ------------------------------------------------------------------------------------------------------------------------
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.value = LW_from0toSinv_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromSinvtoG_new.value = LW_fromSinvtoG_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromSinvtoG_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromGtoF_new.value = LW_fromGtoF_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromGtoF_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoE_new.value = LW_fromFtoE_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoE_new.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.value = LW_fromEto0_new;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.Attributes.unit = "daN"; 
            % =============================================================
            % WING LOADS DIAGRAM - NEGATIVE SIDE   
            figure(6)
            % ---------------------------------------------------------------------
            plot(V_from0toSinv, LW_from0toSinv_new, '-r', 'LineWidth', 1)
            plot(V_fromSinvtoG, LW_fromSinvtoG_new, '-r', 'LineWidth', 1)
            plot(V_fromGtoF,    LW_fromGtoF_new,    '-r', 'LineWidth', 1)
            plot(V_fromFtoE,    LW_fromFtoE_new,    '-r', 'LineWidth', 1)
            plot(V_fromEto0,    LW_fromEto0_new,    '-r', 'LineWidth', 1)
            % ---------------------------------------------------------------------
            plot(V_from0toSinv(1),   LW_from0toSinv_new(1),   'k.', 'MarkerSize', 10)
            plot(V_from0toSinv(end), LW_from0toSinv_new(end), 'k.', 'MarkerSize', 10)
            plot(V_fromSinvtoG(end), LW_fromSinvtoG_new(end), 'k.', 'MarkerSize', 10)
            plot(V_fromGtoF(end),    LW_fromGtoF_new(end),    'k.', 'MarkerSize', 10)
            plot(V_fromFtoE(end),    LW_fromFtoE_new(end),    'k.', 'MarkerSize', 10)
            plot(V_fromEto0(end),    LW_fromEto0_new(end),    'k.', 'MarkerSize', 10)
            % ---------------------------------------------------------------------
            text(V_from0toSinv(end),  LW_from0toSinv_new(end), ' S inv.', 'FontSize', 6)
            text(V_fromSinvtoG(end), LW_fromSinvtoG_new(end),  ' G',      'FontSize', 6)
            text(V_fromGtoF(end),    LW_fromGtoF_new(end),     ' F',      'FontSize', 6)
            text(V_fromFtoE(end),    LW_fromFtoE_new(end),     ' E',      'FontSize', 6)  
            
        % STORE INSIDE THE AIRCRAFT STRUCTURE VARIABLE
        
        % POINT S_INV
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.value = q_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.Attributes.unit = "Pa"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv_new.value = LW_from0toSinv_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LHTS_inv.value = LHT_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LHTS_inv.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv.value = alfa_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_rad.value = deg2rad(alfa_from0toSinv(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_rad.Attributes.unit = "rad";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CLHT_S_inv.value = CLHT_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CLHT_S_inv.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CD_S_inv.value = CD_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CD_S_inv.Attributes.unit = "Non dimensional";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CM_S_inv.value = CMCG_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CM_S_inv.Attributes.unit = "Non dimensional";
        % POINT G
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value = q_fromSinvtoG(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG_new.value = LW_fromSinvtoG_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.value = LHT_fromSinvtoG(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.value = alfa_fromSinvtoG(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_rad.value = deg2rad(alfa_fromSinvtoG(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CLHT_G.value = CLHT_fromSinvtoG(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CLHT_G.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.value = CD_fromSinvtoG(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CM_G.value = CMCG_fromSinvtoG(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CM_G.Attributes.unit = "Non dimensional";
        % POINT F
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value = q_fromGtoF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.Attributes.unit = "Pa";         
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF_new.value = LW_fromGtoF_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.value = LHT_fromGtoF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.value = alfa_fromGtoF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_rad.value = deg2rad(alfa_fromGtoF(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CLHT_F.value = CLHT_fromGtoF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CLHT_F.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.value = CD_fromGtoF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CM_F.value = CMCG_fromGtoF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CM_F.Attributes.unit = "Non dimensional";
        % POINT E
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value = q_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.Attributes.unit = "Pa"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE_new.value = LW_fromFtoE_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.value = LHT_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.value = alfa_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_rad.value = deg2rad(alfa_fromFtoE(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_rad.Attributes.unit = "rad"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CLHT_E.value = CLHT_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CLHT_E.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CD_E.value = CD_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CD_E.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CM_E.value = CMCG_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CM_E.Attributes.unit = "Non dimensional";
        
        end
        
    % CASE 2: Real solutions of the intercept     
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
        % FROM 0 TO S_INVERTED 
        CL_from0toSinv   = zeros(length(V_from0toSinv), 1);
        alfa_from0toSinv = zeros(length(V_from0toSinv), 1);
        CD_from0toSinv   = zeros(length(V_from0toSinv), 1);
        q_from0toSinv    = zeros(length(V_from0toSinv), 1);
        WBL_from0toSinv  = zeros(length(V_from0toSinv), 1);
        CMCL_from0toSinv = zeros(length(V_from0toSinv), 1);
        CMCD_from0toSinv = zeros(length(V_from0toSinv), 1);
        CMCT_from0toSinv = zeros(length(V_from0toSinv), 1);
        CMCG_from0toSinv = zeros(length(V_from0toSinv), 1); 
        CLHT_from0toSinv = zeros(length(V_from0toSinv), 1);
        LHT_from0toSinv  = zeros(length(V_from0toSinv), 1);
        for i = 1:length(V_from0toSinv)
            alfa_from0toSinv(i) = alfa_func(rho0, S, V_from0toSinv(i), WS, n_from0toSinv(i), CLalfa, alpha_zerol);
            CL_from0toSinv(i)   = CL_calc(obj1, n_from0toSinv(i), Mass, g, V_from0toSinv(i), rho0, S);
            if CL_from0toSinv(i) < CL_star
                CL_from0toSinv(i) = CL_calc(obj1, n_from0toSinv(i), Mass, g, V_from0toSinv(i), rho0, S);
            elseif CL_from0toSinv(i) > CL_star
                CL_from0toSinv(i) = CLmax_non_lin(alfa_from0toSinv(i));
            end
            CD_from0toSinv(i)   = cd_calc(obj1, CD0, CL_from0toSinv(i), AR, e, k1, k2);
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
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv.value = CL_from0toSinv;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv.value = alfa_from0toSinv;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_rad.value = deg2rad(alfa_from0toSinv);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_from0toSinv_rad.Attributes.unit = "rad";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toSinv.value = CD_from0toSinv;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_from0toSinv.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toSinv.value = q_from0toSinv;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_from0toSinv.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toSinv.value = WBL_from0toSinv;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_from0toSinv.Attributes.unit = "daN";
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
        alfa_G = alfa_func(rho0, S, VG, WS, nG, CLalfa, alpha_zerol);
        CL_G   = CL_calc(obj1, abs(nG), Mass, g, VG, rho0, S);
            if CL_G < CL_star - 0.03
                CL_G = CL_calc(obj1, abs(nG), Mass, g, VG, rho0, S);
            elseif CL_G > CL_star + 0.03
                CL_G = CLmax_non_lin(alfa_G);
            end
        CD_G   = cd_calc(obj1, CD0, CL_G, AR, e, k1, k2);
        q_G    = 0.5*rho0*(VG)^2;
        WBL_G  = q_G*S*CL_G*1e-1; 
        CMCL_G = CLWB_contrib(obj1, CL_G, deg2rad(alfa_G), XAC, XCG, bCG, MAC);
        CMCD_G = CDWB_contrib(obj1, CL_G, deg2rad(alfa_G), XAC, XCG, bCG, MAC);
        CMCT_G = CT_contr(obj1, CD_G, Thrust_axes, MAC);
        CMCG_G = CM_aboutcg(obj1, CM0, CM_landing_gear, CM_slope, CL_G);
        % HORIZONTAL TAIL LIFT COEFFICIENT
        CLHT_G = CL_Tail(obj1, CMCL_G, CMCD_G, CMCT_G, CMCG_G, l_ht, MAC, XAC, XCG, deg2rad(alfa_G)); 
        % HORIZONTAL TAIL LIFT
        LHT_G = (0.5)*(VG^2)*(S)*(rho0)*(CLHT_G)*(1e-1);
        % NEW LIFT COEFFICIENT
        CL_G_new = CL_G - CLHT_G;
        % NEW LIFT COEFFICIENT
        LW_G_new = q_G*S*CL_G_new*1e-1;
        
        % STORE ALL THE DATA INSIDE THE STRUCT VARIABLE
        % POINT G
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value = q_G;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.Attributes.unit = "Pa";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG_new.value = LW_G_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.value = LHT_G;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.value = alfa_G;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_rad.value = deg2rad(alfa_G);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CLHT_G.value = CLHT_G;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CLHT_G.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.value = CD_G;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CM_G.value = CMCG_G;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CM_G.Attributes.unit = "Non dimensional";          
        % =================================================================  
        % FROM Sinv TO G1  
        CL_fromSinvtoG1   = zeros(length(V_fromSinvtoG1), 1);
        alfa_fromSinvtoG1 = zeros(length(V_fromSinvtoG1), 1);
        CD_fromSinvtoG1   = zeros(length(V_fromSinvtoG1), 1);
        q_fromSinvtoG1    = zeros(length(V_fromSinvtoG1), 1);
        WBL_fromSinvtoG1  = zeros(length(V_fromSinvtoG1), 1);
        CMCL_fromSinvtoG1 = zeros(length(V_fromSinvtoG1), 1);
        CMCD_fromSinvtoG1 = zeros(length(V_fromSinvtoG1), 1);
        CMCT_fromSinvtoG1 = zeros(length(V_fromSinvtoG1), 1);
        CMCG_fromSinvtoG1 = zeros(length(V_fromSinvtoG1), 1); 
        CLHT_fromSinvtoG1 = zeros(length(V_fromSinvtoG1), 1);
        LHT_fromSinvtoG1  = zeros(length(V_fromSinvtoG1), 1);
        for i = 1:length(V_fromSinvtoG1)
            alfa_fromSinvtoG1(i) = alfa_func(rho0, S, V_fromSinvtoG1(i), WS, n_fromSinvtoG1(i), CLalfa, alpha_zerol);
            CL_fromSinvtoG1(i)   = CL_calc(obj1, n_fromSinvtoG1(i), Mass, g, V_fromSinvtoG1(i), rho0, S);
            if CL_fromSinvtoG1(i) < CL_star
                CL_fromSinvtoG1(i) = CL_calc(obj1, n_fromSinvtoG1(i), Mass, g, V_fromSinvtoG1(i), rho0, S);
            elseif CL_fromSinvtoG1(i) > CL_star
                CL_fromSinvtoG1(i) = CLmax_non_lin(alfa_fromSinvtoG1(i));
            end
            CD_fromSinvtoG1(i)   = cd_calc(obj1, CD0, CL_fromSinvtoG1(i), AR, e, k1, k2);
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
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG1.value = CL_fromSinvtoG1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG1.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG1.value = alfa_fromSinvtoG1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG1.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG1_rad.value = deg2rad(alfa_fromSinvtoG1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromSinvtoG1_rad.Attributes.unit = "rad";        
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
        CL_fromG1toF   = zeros(length(V_fromG1toF), 1);
        alfa_fromG1toF = zeros(length(V_fromG1toF), 1);
        CD_fromG1toF   = zeros(length(V_fromG1toF), 1);
        q_fromG1toF    = zeros(length(V_fromG1toF), 1);
        WBL_fromG1toF  = zeros(length(V_fromG1toF), 1);
        CMCL_fromG1toF = zeros(length(V_fromG1toF), 1);
        CMCD_fromG1toF = zeros(length(V_fromG1toF), 1);
        CMCT_fromG1toF = zeros(length(V_fromG1toF), 1);
        CMCG_fromG1toF = zeros(length(V_fromG1toF), 1); 
        CLHT_fromG1toF = zeros(length(V_fromG1toF), 1);
        LHT_fromG1toF  = zeros(length(V_fromG1toF), 1);
        for i = 1:length(V_fromG1toF)
            alfa_fromG1toF(i) = alfa_func(rho0, S, V_fromG1toF(i), WS, n_fromG1toF(i), CLalfa, alpha_zerol);
            CL_fromG1toF(i)   = CL_calc(obj1, n_fromG1toF(i), Mass, g, V_fromG1toF(i), rho0, S);
            if CL_fromG1toF(i) < CL_star
                CL_fromG1toF(i) = CL_calc(obj1, n_fromG1toF(i), Mass, g, V_fromG1toF(i), rho0, S);
            elseif CL_fromG1toF(i) > CL_star
                CL_fromG1toF(i) = CLmax_non_lin(alfa_fromG1toF(i));
            end
            CD_fromG1toF(i)   = cd_calc(obj1, CD0, CL_fromG1toF(i), AR, e, k1, k2);
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
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF.value = CL_fromG1toF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF.value = alfa_fromG1toF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF_rad.value = deg2rad(alfa_fromG1toF);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromG1toF_rad.Attributes.unit = "rad";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromG1toF.value = CD_fromG1toF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromG1toF.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromG1toF.value = q_fromG1toF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromG1toF.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromG1toF.value = WBL_fromG1toF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromG1toF.Attributes.unit = "daN";
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
        CL_fromFtoE   = zeros(length(V_fromFtoE), 1);
        alfa_fromFtoE = zeros(length(V_fromFtoE), 1);
        CD_fromFtoE   = zeros(length(V_fromFtoE), 1);
        q_fromFtoE    = zeros(length(V_fromFtoE), 1);
        WBL_fromFtoE  = zeros(length(V_fromFtoE), 1);
        CMCL_fromFtoE = zeros(length(V_fromFtoE), 1);
        CMCD_fromFtoE = zeros(length(V_fromFtoE), 1);
        CMCT_fromFtoE = zeros(length(V_fromFtoE), 1);
        CMCG_fromFtoE = zeros(length(V_fromFtoE), 1); 
        CLHT_fromFtoE = zeros(length(V_fromFtoE), 1);
        LHT_fromFtoE  = zeros(length(V_fromFtoE), 1);
        for i = 1:length(V_fromFtoE)
            alfa_fromFtoE(i) = alfa_func(rho0, S, V_fromFtoE(i), WS, n_fromFtoE(i), CLalfa, alpha_zerol);
            CL_fromFtoE(i)   = CL_calc(obj1, n_fromFtoE(i), Mass, g, V_fromFtoE(i), rho0, S);
            if CL_fromFtoE(i) < CL_star
                CL_fromFtoE(i) = CL_calc(obj1, n_fromFtoE(i), Mass, g, V_fromFtoE(i), rho0, S);
            elseif CL_fromFtoE(i) > CL_star
                CL_fromFtoE(i) = CLmax_non_lin(alfa_fromFtoE(i));
            end
            CD_fromFtoE(i)   = cd_calc(obj1, CD0, CL_fromFtoE(i), AR, e, k1, k2);
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
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value = CL_fromFtoE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE.value = alfa_fromFtoE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE_rad.value = deg2rad(alfa_fromFtoE);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromFtoE_rad.Attributes.unit = "rad";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value = CD_fromFtoE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromFtoE.value = q_fromFtoE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromFtoE.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE.value = WBL_fromFtoE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE.Attributes.unit = "daN";
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
        CL_fromEto0   = zeros(length(V_fromEto0), 1);
        alfa_fromEto0 = zeros(length(V_fromEto0), 1);
        CD_fromEto0   = zeros(length(V_fromEto0), 1);
        q_fromEto0    = zeros(length(V_fromEto0), 1);
        WBL_fromEto0  = zeros(length(V_fromEto0), 1);
        CMCL_fromEto0 = zeros(length(V_fromEto0), 1);
        CMCD_fromEto0 = zeros(length(V_fromEto0), 1);
        CMCT_fromEto0 = zeros(length(V_fromEto0), 1);
        CMCG_fromEto0 = zeros(length(V_fromEto0), 1); 
        CLHT_fromEto0 = zeros(length(V_fromEto0), 1);
        LHT_fromEto0  = zeros(length(V_fromEto0), 1);
        for i = 1:length(V_fromEto0)
            alfa_fromEto0(i) = alfa_func(rho0, S, V_fromEto0(i), WS, n_fromEto0(i), CLalfa, alpha_zerol);
            CL_fromEto0(i)   = CL_calc(obj1, n_fromEto0(i), Mass, g, V_fromEto0(i), rho0, S);
            if CL_fromEto0(i) < CL_star
                CL_fromEto0(i) = CL_calc(obj1, n_fromEto0(i), Mass, g, V_fromEto0(i), rho0, S);
            elseif CL_fromEto0(i) > CL_star
                CL_fromEto0(i) = CLmax_non_lin(alfa_fromEto0(i));
            end
            CD_fromEto0(i)   = cd_calc(obj1, CD0, CL_fromEto0(i), AR, e, k1, k2);
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
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0.value = CL_fromEto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0.value = alfa_fromEto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_rad.value = deg2rad(alfa_fromEto0);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alfa_fromEto0_rad.Attributes.unit = "rad";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromEto0.value = CD_fromEto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromEto0.Attributes.unit = "Non dimensional"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromEto0.value = q_fromEto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_fromEto0.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromEto0.value = WBL_fromEto0;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromEto0.Attributes.unit = "daN";
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
        
        % MAIN WING LOADS DIAGRAM        
        CL_from0toSinv_new  = CL_from0toSinv  - CLHT_from0toSinv;        
        CL_fromSinvtoG1_new = CL_fromSinvtoG1 - CLHT_fromSinvtoG1;        
        CL_fromG1toF_new    = CL_fromG1toF - CLHT_fromG1toF;        
        CL_fromFtoE_new     = CL_fromFtoE - CLHT_fromFtoE;        
        CL_fromEto0_new     = CL_fromEto0 - CLHT_fromEto0;
        LW_from0toSinv_new  = zeros(length(CL_from0toSinv_new), 1);
        LW_fromSinvtoG1_new = zeros(length(CL_from0toSinv_new), 1);
        LW_fromG1toF_new    = zeros(length(CL_from0toSinv_new), 1);
        LW_fromFtoE_new     = zeros(length(CL_from0toSinv_new), 1);
        LW_fromEto0_new     = zeros(length(CL_from0toSinv_new), 1);
        for i = 1:length(CL_from0toSinv_new)
            LW_from0toSinv_new(i)  = (0.5)*(V_from0toSinv(i)^2)* rho0 * S *(CL_from0toSinv_new(i))*(1e-1);  
            LW_fromSinvtoG1_new(i) = (0.5)*(V_fromSinvtoG1(i)^2)* rho0 * S *(CL_fromSinvtoG1_new(i))*(1e-1);
            LW_fromG1toF_new(i)    = (0.5)*(V_fromG1toF(i)^2)* rho0 * S *(CL_fromG1toF_new(i))*(1e-1);
            LW_fromFtoE_new(i)     = (0.5)*(V_fromFtoE(i)^2)* rho0 * S *(CL_fromFtoE_new(i))*(1e-1);
            LW_fromEto0_new(i)     = (0.5)*(V_fromEto0(i)^2)* rho0 * S *(CL_fromEto0_new(i))*(1e-1);
        end       
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.value = CL_from0toSinv_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_from0toSinv_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG1_new.value = CL_fromSinvtoG1_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromSinvtoG1_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF_new.value = CL_fromG1toF_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromG1toF_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value = CL_fromFtoE_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.value = CL_fromEto0_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromEto0_new.Attributes.unit = "Non dimensional";
        % ------------------------------------------------------------------------------------------------------------------------
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.value = LW_from0toSinv_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toSinv_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromSinvtoG1_new.value = LW_fromSinvtoG1_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromSinvtoG1_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG1toF_new.value = LW_fromG1toF_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromG1toF_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoE_new.value = LW_fromFtoE_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromFtoE_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.value = LW_fromEto0_new;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_fromEto0_new.Attributes.unit = "daN";
        
        % WING LOADS DIAGRAM - NEGATIVE SIDE   
        figure(6)
        % ---------------------------------------------------------------------
        plot(V_from0toSinv,  LW_from0toSinv_new,  '-r', 'LineWidth', 1)
        plot(V_fromSinvtoG1, LW_fromSinvtoG1_new, '-r', 'LineWidth', 1)
        plot(V_fromG1toF,    LW_fromG1toF_new,    '-r', 'LineWidth', 1)
        plot(V_fromFtoE,     LW_fromFtoE_new,     '-r', 'LineWidth', 1)
        plot(V_fromEto0,     LW_fromEto0_new,     '-r', 'LineWidth', 1)
        % ---------------------------------------------------------------------
        plot(V_from0toSinv(1),    LW_from0toSinv_new(1),    'k.', 'MarkerSize', 10)
        plot(V_from0toSinv(end),  LW_from0toSinv_new(end),  'k.', 'MarkerSize', 10)
        plot(VG,                  LW_G_new,                 'k.', 'MarkerSize', 10)
        plot(V_fromSinvtoG1(end), LW_fromSinvtoG1_new(end), 'k.', 'MarkerSize', 10)
        plot(V_fromG1toF(end),    LW_fromG1toF_new(end),    'k.', 'MarkerSize', 10)
        plot(V_fromFtoE(end),     LW_fromFtoE_new(end),     'k.', 'MarkerSize', 10)
        plot(V_fromEto0(end),     LW_fromEto0_new(end),     'k.', 'MarkerSize', 10)
        % ---------------------------------------------------------------------
        text(V_from0toSinv(end),  LW_from0toSinv_new(end),  ' S inv.', 'FontSize', 6)
        text(VG,                  LW_G_new,                 ' G',      'FontSize', 6)
        text(V_fromSinvtoG1(end), LW_fromSinvtoG1_new(end), ' G1',     'FontSize', 6)
        text(V_fromG1toF(end),    LW_fromG1toF_new(end),    ' F',      'FontSize', 6)
        text(V_fromFtoE(end),     LW_fromFtoE_new(end),     ' E',      'FontSize', 6)  
        
        % STORE INSIDE THE AIRCRAFT STRUCTURE VARIABLE
        
        % POINT S_INV
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.value = q_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.Attributes.unit = "Pa";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv_new.value = LW_from0toSinv_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LHTS_inv.value = LHT_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LHTS_inv.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv.value = alfa_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_rad.value = deg2rad(alfa_from0toSinv(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.alfaS_inv_rad.Attributes.unit = "rad"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CLHT_S_inv.value = CLHT_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CLHT_S_inv.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CD_S_inv.value = CD_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CD_S_inv.Attributes.unit = "Non dimensional";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CM_S_inv.value = CMCG_from0toSinv(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CM_S_inv.Attributes.unit = "Non dimensional";     
        % POINT G1
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.qG1.value = q_fromSinvtoG1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.qG1.Attributes.unit = "Pa";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1_new.value = LW_fromSinvtoG1_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LHTG1.value = LHT_fromSinvtoG1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LHTG1.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1.value = alfa_fromSinvtoG1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_rad.value = deg2rad(alfa_fromSinvtoG1(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.alfaG1_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CLHT_G1.value = CLHT_fromSinvtoG1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CLHT_G1.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CD_G1.value = CD_fromSinvtoG1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CD_G1.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CM_G1.value = CMCG_fromSinvtoG1(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CM_G1.Attributes.unit = "Non dimensional";
        % POINT F
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value = q_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.Attributes.unit = "Pa";         
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF_new.value = LW_fromG1toF_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.value = LHT_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.value = alfa_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_rad.value = deg2rad(alfa_fromG1toF(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CLHT_F.value = CLHT_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CLHT_F.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.value = CD_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CM_F.value = CMCG_fromG1toF(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CM_F.Attributes.unit = "Non dimensional";
        % POINT E
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value = q_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.Attributes.unit = "Pa";           
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE_new.value = LW_fromFtoE_new(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE_new.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.value = LHT_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.Attributes.unit = "daN";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.value = alfa_fromFtoE(end);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.Attributes.unit = "degrees";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_rad.value = deg2rad(alfa_fromFtoE(end));
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE_rad.Attributes.unit = "rad";
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
% -------------------------------------------------------------------------
% %% CL CALCULATIONS - POSITIVE LOAD FACTORS
% % ------------------------------------------------------------------------- 
% % Calculation of the CL for the wing - body configuration. It must be
% % noticed that the calculations relative to lift coefficient in this
% % particular range of airspeed V and load factor n are going to give a
% % constant as a results, which is the maximum lift coefficient following
% % the lift curve of the aircrat. The previously defined CLMAX is a little
% % bit greater than the one obtained from V - n diagram flight condition to
% % take into account the extra lift produced by the horizontal empennage.
% 
% % *** From point S to point A ***
% CL_positivestall = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
%     V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value(i);
%     n = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value(i);
%     
%     % A complete documentation of the function CL_calc(...) used here is
%     % inside the class file aero_model.m, which can be found inside the
%     % 'utilities' folder of this library.
%     CL_positivestall(i) = CL_calc(obj1, n, ...
%                     Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
%                     Aircraft.Constants.g.value, ...
%                     V, ...
%                     Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value, ...
%                     Aircraft.Geometry.Wing.S.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value = CL_positivestall;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.Attributes.unit = "Non dimensional";
% 
% % *** From point C to point D ***
% % Now it is necessary to construct the straight line from point C to point
% % fg of the Final Envelope. Then, from fg, we go to Point D.
% 
% % FROM C TO FG
% V_fromCtofg  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.value, ...
%                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_speed.value,   ...
%                         length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value))';
% n_fromCtofg  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.value, ...
%                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_load_factor.value,   ...
%                         length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value))';
% CL_fromCtofg = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
%     CL_fromCtofg(i) = CL_calc(obj1, n_fromCtofg(i), ...
%                     Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
%                     Aircraft.Constants.g.value, ...
%                     V_fromCtofg(i), ...
%                     Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value, ...
%                     Aircraft.Geometry.Wing.S.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value = CL_fromCtofg;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.Attributes.unit = "Non dimensional";
% 
% % FROM FG TO D 
% V_fromfgtoD  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_speed.value, ...
%                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value,   ...
%                         length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value))';
% n_fromfgtoD  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_load_factor.value, ...
%                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value,   ...
%                         length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value))';
% CL_fromfgtoD = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
%     CL_fromfgtoD(i) = CL_calc(obj1, n_fromfgtoD(i), ...
%                     Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
%                     Aircraft.Constants.g.value, ...
%                     V_fromfgtoD(i), ...
%                     Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value, ...
%                     Aircraft.Geometry.Wing.S.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value = CL_fromfgtoD;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.Attributes.unit = "Non dimensional";
% 
% % V_fromCtoD  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.value, ...
% %                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value,   ...
% %                       length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value));
% % n_fromCtoD  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.value, ...
% %                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value,   ...
% %                       length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value));
% % CL_fromCtoD = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
% % for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
% %     CL_fromCtoD(i) = CL_calc(obj1, n_fromCtoD(i), ...
% %                     Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
% %                     Aircraft.Constants.g.value, ...
% %                     V_fromCtoD(i), ...
% %                     Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value, ...
% %                     Aircraft.Geometry.Wing.S.value);
% % end
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoD.value = CL_fromCtoD;
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoD.Attributes.unit = "Non dimensional";
% 
% %% CL CALCULATIONS - NEGATIVE LOAD FACTORS
% % ------------------------------------------------------------------------- 
% % Calculation of the CL for the wing - body configuration
% % *** From point S to point F ***
% CL_negativestall = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value)
%     V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value(i);
%     n = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.value(i);
%     CL_negativestall(i) = CL_calc(obj1, n, ...
%                     Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
%                     Aircraft.Constants.g.value, ...
%                     V, ...
%                     Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value, ...
%                     Aircraft.Geometry.Wing.S.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value = CL_negativestall;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.Attributes.unit = "Non dimensional";
% 
% % *** From point F to point E ***
% % Now it is necessary to construct the straight line from point C to point
% % D of the Final Envelope. 
% V_fromFtoE  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.value, ...
%                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value,   ...
%                       length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value));
% n_fromFtoE  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.value, ...
%                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value,   ...
%                       length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value));
% CL_fromFtoE = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
%     CL_fromFtoE(i) = CL_calc(obj1, n_fromFtoE(i), ...
%                     Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
%                     Aircraft.Constants.g.value, ...
%                     V_fromFtoE(i), ...
%                     Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value, ...
%                     Aircraft.Geometry.Wing.S.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value = CL_fromFtoE;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.Attributes.unit = "Non dimensional";
% 
% %% CL CALCULATION - FROM POINT D TO POINT E
% V_fromDtoE  = repmat(Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value, ...
%                     [length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value),1]);
% n_fromDtoE  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value, ...
%                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value,   ...
%                       length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value));
% CL_fromDtoE = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
%     CL_fromDtoE(i) = CL_calc(obj1, n_fromDtoE(i), ...
%                     Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
%                     Aircraft.Constants.g.value, ...
%                     V_fromDtoE(i), ...
%                     Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value, ...
%                     Aircraft.Geometry.Wing.S.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value = CL_fromDtoE;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.Attributes.unit = "Non dimensional";
% 
% %% CD CALCULATION - POSITIVE LOAD FACTOR
% % *** From point S to point C ***
% CD_positivestall = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
%     
%     % A complete documentation of the function cd_calc(...) used here is
%     % inside the class file aero_model.m, which can be found inside the
%     % 'utilities' folder of this library.   
%     CD_positivestall(i) = cd_calc(obj1, Aircraft.Certification.Aerodynamic_data.CD0.value, ...
%                              Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value(i), ...
%                              Aircraft.Geometry.Wing.AR.value, ...
%                              Aircraft.Certification.Aerodynamic_data.e.value, ...
%                              Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value, ...
%                              Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value = CD_positivestall;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.Attributes.unit = "Non dimensional";
% 
% % *** From point C to point D ***
% 
% % FROM C TO FG
% CD_fromCtofg = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
%     CD_fromCtofg(i) = cd_calc(obj1, Aircraft.Certification.Aerodynamic_data.CD0.value, ...
%                              Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value(i), ...
%                              Aircraft.Geometry.Wing.AR.value, ...
%                              Aircraft.Certification.Aerodynamic_data.e.value, ...
%                              Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value, ...
%                              Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtofg.value = CD_fromCtofg;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtofg.Attributes.unit = "Non dimensional";
% 
% 
% % FROM FG TO D
% CD_fromfgtoD = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
%     CD_fromfgtoD(i) = cd_calc(obj1, Aircraft.Certification.Aerodynamic_data.CD0.value, ...
%                              Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value(i), ...
%                              Aircraft.Geometry.Wing.AR.value, ...
%                              Aircraft.Certification.Aerodynamic_data.e.value, ...
%                              Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value, ...
%                              Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromfgtoD.value = CD_fromfgtoD;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromfgtoD.Attributes.unit = "Non dimensional";
% 
% %% CD CALCULATION - NEGATIVE LOAD FACTOR
% % *** From point S to point F ***
% CD_negativestall = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value)
%     CD_negativestall(i) = cd_calc(obj1, Aircraft.Certification.Aerodynamic_data.CD0.value, ...
%                              Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value(i), ...
%                              Aircraft.Geometry.Wing.AR.value, ...
%                              Aircraft.Certification.Aerodynamic_data.e.value, ...
%                              Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value, ...
%                              Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_negativestall.value = CD_negativestall;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_negativestall.Attributes.unit = "Non dimensional";
% 
% % *** From point F to E ***
% CD_fromFtoE = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value)
%     CD_fromFtoE(i) = cd_calc(obj1, Aircraft.Certification.Aerodynamic_data.CD0.value, ...
%                              Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value(i), ...
%                              Aircraft.Geometry.Wing.AR.value, ...
%                              Aircraft.Certification.Aerodynamic_data.e.value, ...
%                              Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value, ...
%                              Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value = CD_fromFtoE;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.Attributes.unit = "Non dimensional";
% 
% %% CD CALCULATION - FROM POINT D TO POINT E
% 
% CD_fromDtoE = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value)
%     CD_fromDtoE(i) = cd_calc(obj1, Aircraft.Certification.Aerodynamic_data.CD0.value, ...
%                              Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value(i), ...
%                              Aircraft.Geometry.Wing.AR.value, ...
%                              Aircraft.Certification.Aerodynamic_data.e.value, ...
%                              Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value, ...
%                              Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDtoE.value = CD_fromDtoE;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDtoE.Attributes.unit = "Non dimensional";
% 
% %% ALPHA CALCULATION - POSITIVE LOAD FACTOR
% 
% % Interpolation coefficient from actual aerodynamic data
% CL_supp    = str2num(Aircraft.Certification.Aerodynamic_data.CL.value);
% alpha_supp = str2num(Aircraft.Certification.Aerodynamic_data.alpha.value);
% x          = 0.03*ones(length(CL_supp), 1);
% p = polyfit(CL_supp + x, alpha_supp, 2);
% 
% % Normal force curve slope in 1/deg
% % Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value = 0.0913;
% % Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.Attributes.unit = "1/degree";
% 
% % Alpha zero lift
% Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value = alpha_calc_lin(obj1, ...
%                                                                                0.0, ...
%                                                                                Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                                                                Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
% Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.Attributes.unit = "degree";
% 
% % *** From point S to point C ***
% % alpha_positivestall = zeros(length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value), 1);
% % for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value)
%     
%     % A complete documentation of the function alpha_calc(...) used here is
%     % inside the class file aero_model.m, which can be found inside the
%     % 'utilities' folder of this library.      
%     % alpha_positivestall(i) = alpha_calc_lin(obj1, ...
%     %                                    Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_positivestall.value(i), ...
%     %                                    Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%     %                                    Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
%     % Aircraft.Certification.Aerodynamic_data.CL.value, Aircraft.Certification.Aerodynamic_data.alpha.value,
%     alpha_positivestall = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% % end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall.value = alpha_positivestall;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall.Attributes.unit = "Degrees";
% 
% % *** From point C to point D ***
% % alpha_fromCtoD = zeros(length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value), 1);
% % for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value)
% %     alpha_fromCtoD(i) = alpha_calc_lin(obj1, ...
% %                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromCtoD.value(i), ...
% %                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
% 
% % FROM C TO fg
%       alpha_fromCtofg = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% % end
% 
% % Alpha from C to fg
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg.value = alpha_fromCtofg;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg.Attributes.unit = "Degrees";
% 
% % FROM fg TO D
%       alpha_fromfgtoD = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% % end
% 
% % Alpha from fg to D
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD.value = alpha_fromfgtoD;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD.Attributes.unit = "Degrees";
% 
% %% ALPHA CALCULATION - NEGATIVE LOAD FACTOR
% 
% % *** From point S to point F *** 
% % alpha_negativestall = zeros(length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Negative_stall_speed.value), 1);
% % for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Negative_stall_speed.value)
% %     alpha_negativestall(i) = alpha_calc_lin(obj1, ...
% %                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_negativestall.value(i), ...
% %                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
%       alpha_negativestall = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% 
% % end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall.value = alpha_negativestall - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall.Attributes.unit = "Degrees";
% 
% % *** From point F to point E ***
% % alpha_fromFtoE = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromFtoE.value), 1);
% % for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromFtoE.value)
% %     alpha_fromFtoE(i) = alpha_calc_lin(obj1, ...
% %                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromFtoE.value(i), ...
% %                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
%       alpha_fromFtoE = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% % end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE.value = alpha_fromFtoE - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE.Attributes.unit = "Degrees";
% 
% %% ALPHA CALCULATION - FROM POINT D TO POINT E
% % alpha_fromDtoE = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromDtoE.value), 1);
% % for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromDtoE.value)
% %     alpha_fromDtoE(i) = alpha_calc_lin(obj1, ...
% %                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromDtoE.value(i), ...
% %                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
%       alpha_fromDtoE = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% % end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE.value = alpha_fromDtoE - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE.Attributes.unit = "Degrees";
% 
% %% ALPHA CONVERSION IN RADIANS
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall_rad.value           = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall_rad.Attributes.unit = "Radians";
% % *************************************************************************
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_rad.value                = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_rad.Attributes.unit      = "Radians";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_rad.value                = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_rad.Attributes.unit      = "Radians";
% % -------------------------------------------------------------------------
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall_rad.value           = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall_rad.Attributes.unit = "Radians";
% % *************************************************************************
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_rad.value                = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_rad.Attributes.unit      = "Radians";
% % ---------------------------------------------------------------------------
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_rad.value                = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_rad.Attributes.unit      = "Radians";
% 
% %% WING BODY LIFT CALCULATION
% % In this section of the code the following formula is applied to evaluate
% % the wing-body lift with the formula 
% %        
% %   L = 0.5*rho*(V^2)*Sw*CL 
% %
% % using the previously calculated lift coefficients, coming from the 
% % various curves of the final envelope diagram.
% 
% %   Positive stall speed
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_positivestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value)   
%     V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_positivestall.value(i) = (0.5)*(V^2)* ...
%         (Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_positivestall.Attributes.unit = "daN";
% 
% % *** From point C to point D ***
% % FROM C TO FG
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtofg.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value)   
%     V = V_fromCtofg(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtofg.value(i) = (0.5)*(V^2)* ...
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtofg.Attributes.unit = "daN";
% 
% % FROM FG TO D
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromfgtoD.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value)   
%     V = V_fromfgtoD(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromfgtoD.value(i) = (0.5)*(V^2)* ...
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromfgtoD.Attributes.unit = "daN";
% % *************************************************************************
% 
% % Negative stall speed
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_negativestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value)   
%     V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_negativestall.value(i) = (0.5)*(V^2)* ...
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_negativestall.Attributes.unit = "daN";
% 
% % *** From point F to point E ***
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value)   
%     V = V_fromFtoE(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE.value(i) = (0.5)*(V^2)* ...
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE.Attributes.unit = "daN";
% % *************************************************************************
% 
% % *** From point D to point E ***
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value)   
%     V = V_fromDtoE(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDtoE.value(i) = (0.5)*(V^2)* ... 
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDtoE.Attributes.unit = "daN";
% 
% %% PITCHING MOMENT CALCULATIONS - CL WING BODY CONTRIBUTIONS
% % Main lifting surface lift pitching moment contribution - POSITIVE LOAD FACTOR
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_positivestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value)
%     
%     % A complete documentation of the function CLWB_contrib(...) used here 
%     % is inside the class file aero_model.m, which can be found inside the
%     % 'utilities' folder of this library.      
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_positivestall.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value(i), ...
%                    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall_rad.value(i), ...
%                    Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
%                    Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                    Aircraft.Certification.Aerodynamic_data.bcg.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_positivestall.Attributes.unit = "Non dimensional";
% 
% % *** From point C to point D ***
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromCtofg.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromCtofg.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value(i), ...
%                    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_rad.value(i), ...
%                    Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
%                    Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                    Aircraft.Certification.Aerodynamic_data.bcg.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromCtofg.Attributes.unit = "Non dimensional";
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromfgtoD.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromfgtoD.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value(i), ...
%                    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_rad.value(i), ...
%                    Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
%                    Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                    Aircraft.Certification.Aerodynamic_data.bcg.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromfgtoD.Attributes.unit = "Non dimensional";
% 
% % Main lifting surface lift pitching moment contribution - NEGATIVE LOAD FACTOR
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_negativestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_negativestall.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value(i), ...
%                    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall_rad.value(i), ...
%                    Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
%                    Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                    Aircraft.Certification.Aerodynamic_data.bcg.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_negativestall.Attributes.unit = "Non dimensional";
% 
% % *** From point F to point E ***
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromFtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromFtoE.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value(i), ...
%                    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_rad.value(i), ...
%                    Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
%                    Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                    Aircraft.Certification.Aerodynamic_data.bcg.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromFtoE.Attributes.unit = "Non dimensional";
% 
% % Main lifting surface lift moment contribution - FROM POINT D TO POINT E 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromDtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromDtoE.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value(i), ...
%                    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_rad.value(i), ...
%                    Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
%                    Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                    Aircraft.Certification.Aerodynamic_data.bcg.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromDtoE.Attributes.unit = "Non dimensional";
% 
% %% PITCHING MOMENT CALCULATION - CD WING BODY CONTRIBUTIONS
% % Main lifting surface drag pitching moment contribution - POSITIVE LOAD FACTORS 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_positivestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value)
%     
%     % A complete documentation of the function CLWB_contrib(...) used here 
%     % is inside the class file aero_model.m, which can be found inside the
%     % 'utilities' folder of this library.      
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_positivestall.value(i) = CDWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value(i), ...
%                    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall_rad.value(i), ...
%                    Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
%                    Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                    Aircraft.Certification.Aerodynamic_data.bcg.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_positivestall.Attributes.unit = "Non dimensional";
% 
% % *** From point C to point D ***
% 
% % FROM C TO FG
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromCtofg.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtofg.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtofg.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromCtofg.value(i) = CDWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtofg.value(i), ...
%                    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_rad.value(i), ...
%                    Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
%                    Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                    Aircraft.Certification.Aerodynamic_data.bcg.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromCtofg.Attributes.unit = "Non dimensional";
% 
% % FROM FG TO D
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromfgtoD.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromfgtoD.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromfgtoD.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromfgtoD.value(i) = CDWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromfgtoD.value(i), ...
%                    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_rad.value(i), ...
%                    Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
%                    Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                    Aircraft.Certification.Aerodynamic_data.bcg.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromfgtoD.Attributes.unit = "Non dimensional";
% 
% % Main lifting surface drag pitching moment contribution - NEGATIVE LOAD FACTOR
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_negativestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_negativestall.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_negativestall.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_negativestall.value(i) = CDWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_negativestall.value(i), ...
%                    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall_rad.value(i), ...
%                    Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
%                    Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                    Aircraft.Certification.Aerodynamic_data.bcg.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_negativestall.Attributes.unit = "Non dimensional";
% 
% % *** From point F to point E ***
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromFtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromFtoE.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value(i), ...
%                    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_rad.value(i), ...
%                    Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
%                    Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                    Aircraft.Certification.Aerodynamic_data.bcg.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromFtoE.Attributes.unit = "Non dimensional";
% 
% % Main lifting surface drag moment contribution - FROM POINT D TO POINT E 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromDtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDtoE.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDtoE.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromDtoE.value(i) = CDWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value(i), ...
%                    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_rad.value(i), ...
%                    Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
%                    Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                    Aircraft.Certification.Aerodynamic_data.bcg.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromDtoE.Attributes.unit = "Non dimensional";
% 
% %% PITCHING MOMENT CALCULATION - CT CONTRIBUTIONS
% % Non - baricentral thrust - POSITIVE LOAD FACTORS DRAG COEFFICIENTS
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_positivestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value)
%     
%     % A complete documentation of the function CT_contrib(...) used here 
%     % is inside the class file aero_model.m, which can be found inside the
%     % 'utilities' folder of this library.     
%     
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_positivestall.value(i) = CT_contr(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value(i), ...
%                    Aircraft.Geometry.Engine.Primary.Thrust_axes.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_positivestall.Attributes.unit = "Non dimensional";
% 
% % *** From point C to point D ***
% 
% % FROM C TO FG
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromCtofg.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtofg.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtofg.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromCtofg.value(i) = CT_contr(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtofg.value(i), ...
%                    Aircraft.Geometry.Engine.Primary.Thrust_axes.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromCtofg.Attributes.unit = "Non dimensional";
% 
% % FROM FG TO D
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromfgtoD.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromfgtoD.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromfgtoD.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromfgtoD.value(i) = CT_contr(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromfgtoD.value(i), ...
%                    Aircraft.Geometry.Engine.Primary.Thrust_axes.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromfgtoD.Attributes.unit = "Non dimensional";
% 
% % Non - baricentral thrust - NEGATIVE LOAD FACTORS DRAG COEFFICIENTS
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_negativestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_negativestall.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_negativestall.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_negativestall.value(i) = CT_contr(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_negativestall.value(i), ...
%                    Aircraft.Geometry.Engine.Primary.Thrust_axes.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_negativestall.Attributes.unit = "Non dimensional";
% 
% % *** From point F to point E ***
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromFtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromFtoE.value(i) = CT_contr(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value(i), ...
%                    Aircraft.Geometry.Engine.Primary.Thrust_axes.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromFtoE.Attributes.unit = "Non dimensional";
% 
% % Non - baricental thrust - FROM POINT D TO POINT E 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromDtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDtoE.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDtoE.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromDtoE.value(i) = CT_contr(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDtoE.value(i), ...
%                    Aircraft.Geometry.Engine.Primary.Thrust_axes.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromDtoE.Attributes.unit = "Non dimensional";
% 
% %% PITCHING MOMENT COEFFICIENT OF THE WING BODY CONFIGURATION 
% % Wing body pitching moment - POSITIVE LOAD FACTORS
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_positivestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value)
%     
%     % A complete documentation of the function CM_aboutcg(...) used here 
%     % is inside the class file aero_model.m, which can be found inside the
%     % 'utilities' folder of this library.     
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_positivestall.value(i) = CM_aboutcg(obj1, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CM0.value, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CMCL.value, ...
%                                                                                             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value(i));
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_positivestall.Attributes.unit = "Non dimensional";
% 
% % *** From point C to point D ***
% 
% % FROM C TO FG
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromCtofg.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromCtofg.value(i) = CM_aboutcg(obj1, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CM0.value, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CMCL.value, ...
%                                                                                             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value(i));
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromCtofg.Attributes.unit = "Non dimensional";
% 
% % FROM FG TO D
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromfgtoD.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromfgtoD.value(i) = CM_aboutcg(obj1, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CM0.value, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CMCL.value, ...
%                                                                                             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value(i));
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromfgtoD.Attributes.unit = "Non dimensional";
% 
% % Wing body pitching moment - NEGATIVE LOAD FACTORS
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_negativestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_negativestall.value(i) = CM_aboutcg(obj1, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CM0.value, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CMCL.value, ...
%                                                                                             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value(i));
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_negativestall.Attributes.unit = "Non dimensional";
% 
% % *** From point F to point E ***
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromFtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromFtoE.value(i) = CM_aboutcg(obj1, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CM0.value, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CMCL.value, ...
%                                                                                             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value(i));
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromFtoE.Attributes.unit = "Non dimensional";
% 
% % *** From point D to point E ***
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromDtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromDtoE.value(i) = CM_aboutcg(obj1, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CM0.value, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CMCL.value, ...
%                                                                                             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value(i));
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromDtoE.Attributes.unit = "Non dimensional";
% 
% %% HORIZONTAL TAIL LIFT COEFFICIENT CALCULATION
% % POSITIVE LOAD FACTORS
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_positivestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall_rad.value) , 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall_rad.value)
%     
%     % A complete documentation of the function CL_Tail(...) used here 
%     % is inside the class file aero_model.m, which can be found inside the
%     % 'utilities' folder of this library.        
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_positivestall.value(i) = CL_Tail(obj1, ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_positivestall.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_positivestall.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_positivestall.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_positivestall.value(i), ...
%                                                         Aircraft.Geometry.Horizontal.l.value, ...
%                                                         Aircraft.Geometry.Wing.mac.value, ...
%                                                         Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ...
%                                                         Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall_rad.value(i));
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_positivestall.Attributes.unit = "Non dimensional";
%                                                    
% % *** From point C to point D ***
% 
% % FROM C TO FG
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromCtofg.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_rad.value) , 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_rad.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromCtofg.value(i) = CL_Tail(obj1, ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromCtofg.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromCtofg.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromCtofg.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromCtofg.value(i), ...
%                                                         Aircraft.Geometry.Horizontal.l.value, ...
%                                                         Aircraft.Geometry.Wing.mac.value, ...
%                                                         Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ...
%                                                         Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_rad.value(i));
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromCtofg.Attributes.unit = "Non dimensional";
% 
% % FROM FG TO D
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromfgtoD.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_rad.value) , 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_rad.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromfgtoD.value(i) = CL_Tail(obj1, ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromfgtoD.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromfgtoD.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromfgtoD.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromfgtoD.value(i), ...
%                                                         Aircraft.Geometry.Horizontal.l.value, ...
%                                                         Aircraft.Geometry.Wing.mac.value, ...
%                                                         Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ...
%                                                         Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_rad.value(i));
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromfgtoD.Attributes.unit = "Non dimensional";
% 
% % NEGATIVE LOAD FACTORS
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_negativestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall_rad.value) , 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall_rad.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_negativestall.value(i) = CL_Tail(obj1, ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_negativestall.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_negativestall.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_negativestall.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_negativestall.value(i), ...
%                                                         Aircraft.Geometry.Horizontal.l.value, ...
%                                                         Aircraft.Geometry.Wing.mac.value, ...
%                                                         Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ...
%                                                         Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall_rad.value(i));
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_negativestall.Attributes.unit = "Non dimensional";
% 
% % *** From point F to point E ***
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromFtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_rad.value) , 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_rad.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromFtoE.value(i) = CL_Tail(obj1, ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromFtoE.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromFtoE.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromFtoE.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromFtoE.value(i), ...
%                                                         Aircraft.Geometry.Horizontal.l.value, ...
%                                                         Aircraft.Geometry.Wing.mac.value, ...
%                                                         Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ...
%                                                         Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_rad.value(i));
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromFtoE.Attributes.unit = "Non dimensional";
% 
% % *** From point D to point E ***
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromDtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_rad.value) , 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_rad.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromDtoE.value(i) = CL_Tail(obj1, ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromDtoE.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromDtoE.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromDtoE.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromDtoE.value(i), ...
%                                                         Aircraft.Geometry.Horizontal.l.value, ...
%                                                         Aircraft.Geometry.Wing.mac.value, ...
%                                                         Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ...
%                                                         Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_rad.value(i));
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromDtoE.Attributes.unit = "Non dimensional";
% 
% %% HORIZONTAL TAIL DIMENSIONAL AIRLOADS
% % Positive stall speed
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_positivestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_positivestall.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_positivestall.value)   
%     V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_positivestall.value(i) = (0.5)*(V^2)* ...
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_positivestall.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_positivestall.Attributes.unit = "daN";
% 
% % *** From point C to point D ***
% % FROM C TO FG
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromCtofg.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromCtofg.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromCtofg.value)   
%     V = V_fromCtofg(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromCtofg.value(i) = (0.5)*(V^2)* ... 
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromCtofg.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromCtofg.Attributes.unit = "daN";
% % FROM FG TO D
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromfgtoD.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromfgtoD.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromfgtoD.value)   
%     V = V_fromfgtoD(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromfgtoD.value(i) = (0.5)*(V^2)* ... 
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromfgtoD.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromfgtoD.Attributes.unit = "daN";
% % *************************************************************************
% 
% % Negative stall speed
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_negativestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_negativestall.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_negativestall.value)   
%     V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_negativestall.value(i) = (0.5)*(V^2)* ...
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_negativestall.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_negativestall.Attributes.unit = "daN";
% 
% % *** From point F to point E ***
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromFtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromFtoE.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value)   
%     V = V_fromFtoE(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromFtoE.value(i) = (0.5)*(V^2)* ... 
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromFtoE.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromFtoE.Attributes.unit = "daN";
% % *************************************************************************
% 
% % *** From point D to point E ***
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromDtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromDtoE.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromDtoE.value)   
%     V = V_fromDtoE(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromDtoE.value(i) = (0.5)*(V^2)* ...
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromDtoE.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromDtoE.Attributes.unit = "daN";
% 
%% TAIL BALANCING LOADS - N = 1.0 
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
% VS = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value(1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value(1) = VS;
% VD = V_fromfgtoD(end);
% for i = 2:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value)
%         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value(i-1) + (VD - VS)*(1/length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value));
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.Attributes.unit = "m/s";
% 
% % DYNAMIC PRESSURE
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.value = 0.5*Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value.^2);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.Attributes.unit = "Pa";
% 
% % WING BODY LIFT COEFFICIENT
% Dyn_press_force = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.value*Aircraft.Geometry.Wing.S.value;
% Aircraft_weight = ones(length(Dyn_press_force),1)*Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value = zeros(length(Dyn_press_force), 1);
% for i = 1:length(Dyn_press_force)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value(i) = ((Aircraft_weight(i))*(1/Dyn_press_force(i)));
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.Attributes.unit = "Non dimensional";
% 
% % DRAG LIFT COEFFICIENT
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value(i) = cd_calc(obj1, Aircraft.Certification.Aerodynamic_data.CD0.value, ...
%                              Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value(i), ...
%                              Aircraft.Geometry.Wing.AR.value, ...
%                              Aircraft.Certification.Aerodynamic_data.e.value, ...
%                              Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value, ...
%                              Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.Attributes.unit = "Non dimensional";
% 
% % ALPHA CALCULATION 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromDtoE.value), 1);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.value = alpha_calc(obj1, ...
%                                                    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value, ...
%                                                    Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                                    Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                                    Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                                    p);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.Attributes.unit = "degrees";
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.Attributes.unit = "radians";
% 
% % WING BODY LIFT 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value)   
%     V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.value(i) = (0.5)*(V^2)* ...
%         (Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.Attributes.unit = "daN";
% 
% obj1 = aero_model;
% % EVALUATION OF TAIL LOADS CONTRIBUTION
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value),1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_unit_load_factor.value)
% %     y  = Aircraft.Certification.Aerodynamic_data.XAC_nondim.value - Aircraft.Certification.Aerodynamic_data.XCG_nondim.value;
% %     x  = Aircraft.Certification.Aerodynamic_data.bcg.value/Aircraft.Geometry.Wing.mac.value;
% %     t1 = cos(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value(i));
% %     t2 = sin(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value(i));
% %     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_unit_load_factor.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value(i)*t1*y - Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value(i)*2*(x);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_unit_load_factor.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value(i), ...
%                    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value(i), ...
%                    Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
%                    Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                    Aircraft.Certification.Aerodynamic_data.bcg.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_unit_load_factor.Attributes.unit = "Non dimensional";
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromFtoE.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value(i), ...
%                    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value(i), ...
%                    Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
%                    Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                    Aircraft.Certification.Aerodynamic_data.bcg.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_unit_load_factor.Attributes.unit = "Non dimensional";
% 
% % *** From point C to point D ***
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_unit_load_factor.value(i) = CT_contr(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value(i), ...
%                    Aircraft.Geometry.Engine.Primary.Thrust_axes.value, ...
%                    Aircraft.Geometry.Wing.mac.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_unit_load_factor.Attributes.unit = "Non dimensional";
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_unit_load_factor.value(i) = CM_aboutcg(obj1, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CM0.value, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value, ...
%                                                                                             Aircraft.Certification.Aerodynamic_data.CMCL.value, ...
%                                                                                             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value(i));
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_unit_load_factor.Attributes.unit = "Non dimensional";
% 
% % TAIL LIFT COEFFICIENT 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value) , 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value)
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_unit_load_factor.value(i) = CL_Tail(obj1, ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_unit_load_factor.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_unit_load_factor.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_unit_load_factor.value(i), ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_unit_load_factor.value(i), ...
%                                                         Aircraft.Geometry.Horizontal.l.value, ...
%                                                         Aircraft.Geometry.Wing.mac.value, ...
%                                                         Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ...
%                                                         Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
%                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value(i));
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_unit_load_factor.Attributes.unit = "Non dimensional";
% 
% % *** TAIL LIFT AT N = 1.0 ***
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_unit_load_factor.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_unit_load_factor.value)   
%     V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_unit_load_factor.value(i) = (0.5)*(V^2)* ... 
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_unit_load_factor.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_unit_load_factor.Attributes.unit = "daN";
% 
% % *** WING LIFT AT N = 1.0 ***
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Wing_Lift_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_unit_load_factor.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_unit_load_factor.value) 
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Wing_Lift_unit_load_factor.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.value(i) - Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_unit_load_factor.value(i);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Wing_Lift_unit_load_factor.Attributes.unit = "daN";
% 
% %% BALANCING LOADS DIAGRAM 
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
% % -------------------------------------------------------------------------
% % To plot the airloads acting on the horizontal surface (which are
% % necessary to evaluate the whole aircraft Lift, which must be used to size
% % the wing structure) the following function will be used: 
% %
% % HTailAirloadsDiagram = Balancing_loads(HT_Lift_posstall, VSpos, ...
% %                                        HT_Lift_negstall, VSneg, ... 
% %                                        HT_Lift_fromCtoD, V_fromCtoD, ...
% %                                        HT_Lift_fromDtoE, V_fromDtoE, ...
% %                                        HT_Lift_fromFtoE, V_fromFtoE, ...
% %                                        HTail_Lift_unit_load_factor.value, V_unit_load_factor, ... 
% %                                        Reg, Aircraft_name)
% %
% % A complete documentation of this function is included inside the file
% % csvla.m
% 
% disp(" ")
% disp(" ++++ FIGURE 7 - HORIZONTAL EMPENNAGE AIRLOADS ++++ ");
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTailAirloadsDiagram.value = Balancing_loads(obj, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_positivestall.value, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_negativestall.value, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value, ... 
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromCtofg.value, ...
%                 V_fromCtofg, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromfgtoD.value, ...
%                 V_fromfgtoD, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromDtoE.value, ...
%                 V_fromDtoE, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromFtoE.value, ...
%                 V_fromFtoE, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_unit_load_factor.value, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.value, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_Design_manoeuvring_speed_VG.value, ...
%                 Aircraft.Certification.Regulation.value, ... 
%                 Aircraft.Certification.Aircraft_Name.value);
% pause(1);
% % Saving figures inside correct folder
% fprintf('Saving Balancingloads.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile Balancingloads.pdf Output
% movefile Balancingloads.png Output
% 
% %% MAIN WING LOADS DIAGRAM 
% %   In this section we have to take into account the lift coefficient
% %   produced by the horizontal tail. The global lift coefficient will be
% %   the algebric sum of the Wing - Body lift coefficient and the horizontal
% %   tail lift coefficient.
% 
% % LIFT COEFFICIENT CALCULATIONS
% % POSITIVE STALL
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_posstall_new.value = ...
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value - Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_positivestall.value;
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_posstall_new.Attributes.unit = "Non dimensional";
% 
% % FROM C TO D
% % FROM C TO FG
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg_new.value = ...
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value - Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromCtofg.value;
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg_new.Attributes.unit = "Non dimensional";
% 
% % FROM FG TO D
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD_new.value = ...
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value - Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromfgtoD.value;
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD_new.Attributes.unit = "Non dimensional";
% 
% % NEGATIVE STALL
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negstall_new.value = ...
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value - Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_negativestall.value;
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negstall_new.Attributes.unit = "Non dimensional";
% 
% % FROM F TO E 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value = ...
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value - Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromFtoE.value;
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.Attributes.unit = "Non dimensional";
% 
% % FROM D TO E
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE_new.value = ...
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value - Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromDtoE.value;
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE_new.Attributes.unit = "Non dimensional";
% 
% %% LIFT CALCULATIONS WITH NEW VALUES OF THE LIFT COEFFICIENT
% % Positive stall speed
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_posstall_new.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_posstall_new.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_posstall_new.value)   
%     V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_posstall_new.value(i) = (0.5)*(V^2)* ...
%         (Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_posstall_new.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_posstall_new.Attributes.unit = "daN";
% 
% % *** From point C to point D ***
% 
% % FROM C TO FG
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtofg_new.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg_new.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg_new.value)   
%     V = V_fromCtofg(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtofg_new.value(i) = (0.5)*(V^2)* ...
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg_new.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtofg_new.Attributes.unit = "daN";
% 
% % FROM FG TO Ds
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromfgtoD_new.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD_new.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD_new.value)   
%     V = V_fromfgtoD(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromfgtoD_new.value(i) = (0.5)*(V^2)* ...
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD_new.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromfgtoD_new.Attributes.unit = "daN";
% % *************************************************************************
% 
% % Negative stall speed
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_negstall_new.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negstall_new.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negstall_new.value)   
%     V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_negstall_new.value(i) = (0.5)*(V^2)* ...
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negstall_new.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_negstall_new.Attributes.unit = "daN";
% 
% % *** From point F to point E ***
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE_new.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value)   
%     V = V_fromFtoE(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE_new.value(i) = (0.5)*(V^2)* ...
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE_new.Attributes.unit = "daN";
% % *************************************************************************
% 
% % *** From point D to point E ***
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDtoE_new.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE_new.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE_new.value)   
%     V = V_fromDtoE(i);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDtoE_new.value(i) = (0.5)*(V^2)* ... 
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE_new.value(i))*(1E-1);   
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDtoE_new.Attributes.unit = "daN";
% 
% %% MAIN WING AIRLOADS DIAGRAM
% % To plot the airloads acting on the main lifting surface (evaluated as the
% % sum of the wing body Lift plus the horizontal tail lift at equilibrium)
% % the following function will be used:
% %
% % Mainwing_loads(obj, WING_Lift_posstall, VSpos, ...
% %                     WING_Lift_negstall, VSneg, ... 
% %                     WING_Lift_fromCtoD, V_fromCtoD, ...
% %                     WING_Lift_fromDtoE, V_fromDtoE, ...
% %                     WING_Lift_fromFtoE, V_fromFtoE, ...
% %                     Reg, Aircraft_name)
% %
% % A complete documentation of this function is included inside the file
% % csvla.m
% 
% disp(" ")
% disp(" ++++ FIGURE 8 - MAIN WING AIRLOADS ++++ ");                                
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WingAirloadsDiagram.value = Mainwing_loads(obj, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_posstall_new.value, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_negstall_new.value, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value, ... 
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtofg_new.value, ...
%                 V_fromCtofg, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromfgtoD_new.value, ...
%                 V_fromfgtoD, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDtoE_new.value, ...
%                 V_fromDtoE, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE_new.value, ...
%                 V_fromFtoE, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Wing_Lift_unit_load_factor.value, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.value, ...
%                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_Design_manoeuvring_speed_VG.value, ...
%                 Aircraft.Certification.Regulation.value, ... 
%                 Aircraft.Certification.Aircraft_Name.value);
% % Saving figures inside correct folder
% fprintf('Saving Wingairloads.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile Wingairloads.pdf Output
% movefile Wingairloads.png Output
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
% 
% %% ALPHA CALCULATION WITH NEW LIFT COEFFICIENT - POSITIVE LOAD FACTOR
% 
% % *** From point S to point C ***
% % alpha_posstall_new = zeros(length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value), 1);
% % for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value)
% %     alpha_posstall_new(i) = alpha_calc_lin(obj1, ...
% %                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_posstall_new.value(i), ...
% %                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
%       alpha_posstall_new = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_posstall_new.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% 
% % end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_posstall_new.value = alpha_posstall_new;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_posstall_new.Attributes.unit = "Degrees";
% 
% % *** From point C to point D ***
% % alpha_fromCtoD_new = zeros(length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value), 1);
% % for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value)
% %     alpha_fromCtoD_new(i) = alpha_calc_lin(obj1, ...
% %                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromCtoD_new.value(i), ...
% %                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
% 
% % FROM C TO FG
%       alpha_fromCtofg_new = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg_new.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_new.value = alpha_fromCtofg_new;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_new.Attributes.unit = "Degrees";
% 
% % FROM FG TO D 
%       alpha_fromfgtoD_new = alpha_calc(obj1, ...
%                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD_new.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_new.value = alpha_fromfgtoD_new;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_new.Attributes.unit = "Degrees";
% 
% %% ALPHA CALCULATION WITH NEW LIFT COEFFICIENT - NEGATIVE LOAD FACTOR
% 
% % *** From point S to point F *** 
% % alpha_negstall_new = zeros(length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Negative_stall_speed.value), 1);
% % for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Negative_stall_speed.value)
% %     alpha_negstall_new(i) = alpha_calc_lin(obj1, ...
% %                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_negstall_new.value(i), ...
% %                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
%       alpha_negstall_new = alpha_calc(obj1, ...
%                                        abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negstall_new.value), ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
%                                    
% % end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negstall_new.value = alpha_negstall_new - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negstall_new.Attributes.unit = "Degrees";
% 
% % *** From point F to point E ***
% % alpha_fromFtoE_new = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromFtoE_new.value), 1);
% % for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromFtoE_new.value)
% %     alpha_fromFtoE_new(i) = alpha_calc_lin(obj1, ...
% %                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromFtoE_new.value(i), ...
% %                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
%       alpha_fromFtoE_new = alpha_calc(obj1, ...
%                                        abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value), ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% 
% % end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_new.value = alpha_fromFtoE_new - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_new.Attributes.unit = "Degrees";
% 
% %% ALPHA CALCULATION WITH NEW LIFT COEFFICIENT - FROM POINT D TO POINT E
% % alpha_fromDtoE_new = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromDtoE_new.value), 1);
% % for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromDtoE_new.value)
% %     alpha_fromDtoE_new(i) = alpha_calc_lin(obj1, ...
% %                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromDtoE_new.value(i), ...
% %                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
% %                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
%       alpha_fromDtoE_new = alpha_calc(obj1, ...
%                                        abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value), ...
%                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%                                        p);
% 
% % end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_new.value = alpha_fromDtoE_new - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_new.Attributes.unit = "Degrees";
% 
% %% ALPHA_NEW CONVERSION IN RADIANS
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_posstall_new_rad.value           = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_posstall_new.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_posstall_new_rad.Attributes.unit = "Radians";
% % *************************************************************************
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_new_rad.value           = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_new.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_new_rad.Attributes.unit = "Radians";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_new_rad.value           = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_new.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_new_rad.Attributes.unit = "Radians";
% % -------------------------------------------------------------------------
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negstall_new_rad.value           = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negstall_new.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negstall_new_rad.Attributes.unit = "Radians";
% % *************************************************************************
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_new_rad.value           = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_new.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_new_rad.Attributes.unit = "Radians";
% % ---------------------------------------------------------------------------
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_new_rad.value           = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_new.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_new_rad.Attributes.unit = "Radians";
% 
% %% WING GEOMETRY FOR OPEN VSP
% 
% % Define all the geometrical panels for the Open VSP calculations; in
% % general the wing will be subdivided in a certain number of panels (three
% % in our actual case) with different values of the chords and different
% % span.
% 
% % Uncomment this section for a general, trapezoidal wing
% % Aircraft.Geometry.Wing.Kinks.Croot1.value = Aircraft.Geometry.Wing.croot.value;
% % Aircraft.Geometry.Wing.Kinks.Croot1.Attributes.unit = "m";
% % Aircraft.Geometry.Wing.Kinks.Ctip3.value = Aircraft.Geometry.Wing.ctip.value;
% % Aircraft.Geometry.Wing.Kinks.Ctip3.Attributes.unit = "m";
% % Aircraft.Geometry.Wing.Kinks.Croot2.value = NaN;
% % Aircraft.Geometry.Wing.Kinks.Croot2.Attributes.unit = "m";
% % Aircraft.Geometry.Wing.Kinks.Ctip1.value = NaN;
% % Aircraft.Geometry.Wing.Kinks.Ctip1.Attributes.unit = "m";
% % Aircraft.Geometry.Wing.Kinks.Croot3.value = NaN;
% % Aircraft.Geometry.Wing.Kinks.Croot3.Attributes.unit = "m";
% % Aircraft.Geometry.Wing.Kinks.Ctip2.value = NaN;
% % Aircraft.Geometry.Wing.Kinks.Ctip2.Attributes.unit = "m";
% % Aircraft.Geometry.Wing.Kinks.b1.value = Aircraft.Geometry.Wing.b*(1/3);
% % Aircraft.Geometry.Wing.Kinks.b1.Attributes.unit = "m"; 
% % Aircraft.Geometry.Wing.Kinks.b2.value = Aircraft.Geometry.Wing.b*(1/3);
% % Aircraft.Geometry.Wing.Kinks.b2.Attributes.unit = "m"; 
% % Aircraft.Geometry.Wing.Kinks.b3.value = Aircraft.Geometry.Wing.b*(1/3);
% % Aircraft.Geometry.Wing.Kinks.b3.Attributes.unit = "m"; 
% 
% % Chord and span for the panels in our actual case 
% Aircraft.Geometry.Kinks.Croot1.value = Aircraft.Geometry.Wing.croot.value;
% Aircraft.Geometry.Kinks.Croot1.Attributes.unit = "m";
% Aircraft.Geometry.Kinks.Croot2.value = Aircraft.Geometry.Wing.croot.value;
% Aircraft.Geometry.Kinks.Croot2.Attributes.unit = "m";
% Aircraft.Geometry.Kinks.Ctip1.value = Aircraft.Geometry.Kinks.Croot2.value;
% Aircraft.Geometry.Kinks.Ctip1.Attributes.unit = "m";
% Aircraft.Geometry.Kinks.Croot3.value = Aircraft.Geometry.Wing.croot.value;
% Aircraft.Geometry.Kinks.Croot3.Attributes.unit = "m";
% Aircraft.Geometry.Kinks.Ctip2.value = Aircraft.Geometry.Kinks.Croot3.value;
% Aircraft.Geometry.Kinks.Ctip2.Attributes.unit = "m";
% Aircraft.Geometry.Kinks.Ctip3.value = Aircraft.Geometry.Wing.ctip.value;
% Aircraft.Geometry.Kinks.Ctip3.Attributes.unit = "m";
% Aircraft.Geometry.Kinks.b1.value = Aircraft.Geometry.Wing.b.value*(1/3);
% Aircraft.Geometry.Kinks.b1.Attributes.unit = "m"; 
% Aircraft.Geometry.Kinks.b2.value = Aircraft.Geometry.Wing.b.value*(1/3);
% Aircraft.Geometry.Kinks.b2.Attributes.unit = "m"; 
% Aircraft.Geometry.Kinks.b3.value = Aircraft.Geometry.Wing.b.value*(1/3);
% Aircraft.Geometry.Kinks.b3.Attributes.unit = "m"; 
