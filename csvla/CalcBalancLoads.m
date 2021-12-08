% ==== USEFUL FUNCTION DEFINED LOCALLY ====
% -------------------------------------------------------------------------
% CLMAX FUNCTION
CLmax_func = @(rho, S, V, WS, n) (2 / rho) * (1 / V^2) * (WS) * n;
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
obj1 = aero_model; 

%% LIFT CHARACTERISTIC
% Before to start the balancing loads analysis, it is required to plot lift
% coefficient values calculated with a non-linear formulation. Then, linear
% and non-linear calculation are compared with raw data available. The
% selected range is -2.0<CL<2.0. To find a complete documentation of the
% function 'CL_Non_linear_model(...)' search inside aero_model.m file. 

% =========================================================================
% DRAG POLYNOMIAL COEFFICIENTS DEFINITION 
% =========================================================================
Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value = 0.079;                       % Coefficient inside an expression for the CD in polynomial form
Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value = 0.365;                       % Coefficient inside an expressione for the CD in polynomial form
Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.Attributes.unit = "Non dimensional"; 
k1 = Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value;
k2 = Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value
% =========================================================================
% =========================================================================

% NUMBER OF ELEMENTS
numb = 1e3;

% Alpha_star and Alpha_max 
a = -0.021837866;
b = 0.436064773;
c = -0.56312855;
Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a.value = a;
Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b.value = b;
Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c.value = c;

% Alpha_star and Alpha_max 
% a = Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a.value;
% b = Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b.value;
% c = Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c.value;

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CL_WB_model = @(alpha) a*alpha.^2 + b*alpha + c; 
alpha_plus  = @(CL) (-b + sqrt(b^2 - 4*a*(c - CL)))/(2*a);
alpha_meno  = @(CL) (-b - sqrt(b^2 - 4*a*(c - CL)))/(2*a);
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
alpha_star = alpha_plus(Aircraft.Certification.Aerodynamic_data.CL_star.value);
alpha_max  = alpha_meno(Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value);

% Lift coefficient curve
Aircraft.Certification.Aerodynamic_data.AOA_aux.value = linspace(-8.0, 25*alpha_star*1e-3 - alpha_star, 1000)';
Aircraft.Certification.Aerodynamic_data.AOA_aux.Attributes.unit = "degree";
Aircraft.Certification.Aerodynamic_data.AOA_aux1.value = linspace(alpha_star + 50*alpha_star*1e-3, 13.0 , 1000)';
Aircraft.Certification.Aerodynamic_data.AOA_aux1.Attributes.unit = "degrees";

% % Non linear lift curve
% Aircraft.Certification.Aerodynamic_data.CL_NonLinear.value = CL_Non_linear_model(obj1, ...
%     Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a.value, ...
%     Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b.value, ...
%     Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c.value);
% Aircraft.Certification.Aerodynamic_data.CL_NonLinear.Attributes.unit = "Non dimensional";

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
% Lift model comparison
% disp(" ++++ FIGURE 5 - SELECTED INTERPOLATION CURVES FOR AERO DATA ++++ ");
% Aircraft.Certification.Aerodynamic_data.CL_comparison_diagram.value = Lift_comparison(Aircraft.Certification.Aerodynamic_data.AOA_aux.value, ... 
%                   Aircraft.Certification.Aerodynamic_data.CL_NonLinear.value, ...
%                   Aircraft.Certification.Aerodynamic_data.CL_Linear.value, ...
%                   Aircraft.Certification.Aerodynamic_data.CL.value, ...
%                   Aircraft.Certification.Aerodynamic_data.alpha.value);
% pause(1);
% exportgraphics(Aircraft.Certification.Aerodynamic_data.CL_comparison_diagram.value, 'LiftComparison.pdf', 'ContentType', 'vector');
% % Saving figures inside correct folder
% fprintf('Saving LiftComparison.pdf in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile LiftComparison.pdf Output
              
% % Lift curve - Full model
% Aircraft.Certification.Aerodynamic_data.CL_Full_model.value = CL_WB_CompleteCurve(obj1, ...
%     Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%     Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%     Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a.value, ...
%     Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b.value, ...
%     Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c.value, ...
%     Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%     Aircraft.Certification.Aerodynamic_data.AOA_aux.value);
% Aircraft.Certification.Aerodynamic_data.CL_Full_model.Attributes.unit = "Non dimensional";
% 
% % Lift curve - Full model - Inverted flight
% Aircraft.Certification.Aerodynamic_data.CL_Full_model_invertedflight.value =  -(1/Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value) + Aircraft.Certification.Aerodynamic_data.CL_Full_model.value;
% 
% % Complete lift curve diagram 
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
              
%% CL CALCULATIONS - POSITIVE LOAD FACTORS
% ------------------------------------------------------------------------- 
% Calculation of the CL for the wing - body configuration. It must be
% noticed that the calculations relative to lift coefficient in this
% particular range of airspeed V and load factor n are going to give a
% constant as a results, which is the maximum lift coefficient following
% the lift curve of the aircrat. The previously defined CLMAX is a little
% bit greater than the one obtained from V - n diagram flight condition to
% take into account the extra lift produced by the horizontal empennage.

% *** From point S to point A ***
CL_positivestall = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
    V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value(i);
    n = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value(i);
    
    % A complete documentation of the function CL_calc(...) used here is
    % inside the class file aero_model.m, which can be found inside the
    % 'utilities' folder of this library.
    CL_positivestall(i) = CL_calc(obj1, n, ...
                    Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
                    Aircraft.Constants.g.value, ...
                    V, ...
                    Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value, ...
                    Aircraft.Geometry.Wing.S.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value = CL_positivestall;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.Attributes.unit = "Non dimensional";

% *** From point C to point D ***
% Now it is necessary to construct the straight line from point C to point
% fg of the Final Envelope. Then, from fg, we go to Point D.

% FROM C TO FG
V_fromCtofg  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.value, ...
                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_speed.value,   ...
                        length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value))';
n_fromCtofg  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.value, ...
                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_load_factor.value,   ...
                        length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value))';
CL_fromCtofg = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
    CL_fromCtofg(i) = CL_calc(obj1, n_fromCtofg(i), ...
                    Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
                    Aircraft.Constants.g.value, ...
                    V_fromCtofg(i), ...
                    Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value, ...
                    Aircraft.Geometry.Wing.S.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value = CL_fromCtofg;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.Attributes.unit = "Non dimensional";

% FROM FG TO D 
V_fromfgtoD  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_speed.value, ...
                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value,   ...
                        length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value))';
n_fromfgtoD  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_load_factor.value, ...
                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value,   ...
                        length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value))';
CL_fromfgtoD = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
    CL_fromfgtoD(i) = CL_calc(obj1, n_fromfgtoD(i), ...
                    Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
                    Aircraft.Constants.g.value, ...
                    V_fromfgtoD(i), ...
                    Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value, ...
                    Aircraft.Geometry.Wing.S.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value = CL_fromfgtoD;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.Attributes.unit = "Non dimensional";

% V_fromCtoD  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.value, ...
%                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value,   ...
%                       length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value));
% n_fromCtoD  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.value, ...
%                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value,   ...
%                       length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value));
% CL_fromCtoD = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
%     CL_fromCtoD(i) = CL_calc(obj1, n_fromCtoD(i), ...
%                     Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
%                     Aircraft.Constants.g.value, ...
%                     V_fromCtoD(i), ...
%                     Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value, ...
%                     Aircraft.Geometry.Wing.S.value);
% end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoD.value = CL_fromCtoD;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtoD.Attributes.unit = "Non dimensional";

%% CL CALCULATIONS - NEGATIVE LOAD FACTORS
% ------------------------------------------------------------------------- 
% Calculation of the CL for the wing - body configuration
% *** From point S to point F ***
CL_negativestall = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value)
    V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value(i);
    n = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.value(i);
    CL_negativestall(i) = CL_calc(obj1, n, ...
                    Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
                    Aircraft.Constants.g.value, ...
                    V, ...
                    Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value, ...
                    Aircraft.Geometry.Wing.S.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value = CL_negativestall;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.Attributes.unit = "Non dimensional";

% *** From point F to point E ***
% Now it is necessary to construct the straight line from point C to point
% D of the Final Envelope. 
V_fromFtoE  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.value, ...
                      Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value,   ...
                      length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value));
n_fromFtoE  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.value, ...
                      Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value,   ...
                      length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value));
CL_fromFtoE = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
    CL_fromFtoE(i) = CL_calc(obj1, n_fromFtoE(i), ...
                    Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
                    Aircraft.Constants.g.value, ...
                    V_fromFtoE(i), ...
                    Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value, ...
                    Aircraft.Geometry.Wing.S.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value = CL_fromFtoE;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.Attributes.unit = "Non dimensional";

%% CL CALCULATION - FROM POINT D TO POINT E
V_fromDtoE  = repmat(Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value, ...
                    [length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value),1]);
n_fromDtoE  = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value, ...
                      Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value,   ...
                      length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value));
CL_fromDtoE = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
    CL_fromDtoE(i) = CL_calc(obj1, n_fromDtoE(i), ...
                    Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
                    Aircraft.Constants.g.value, ...
                    V_fromDtoE(i), ...
                    Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value, ...
                    Aircraft.Geometry.Wing.S.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value = CL_fromDtoE;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.Attributes.unit = "Non dimensional";

%% CD CALCULATION - POSITIVE LOAD FACTOR
% *** From point S to point C ***
CD_positivestall = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
    
    % A complete documentation of the function cd_calc(...) used here is
    % inside the class file aero_model.m, which can be found inside the
    % 'utilities' folder of this library.   
    CD_positivestall(i) = cd_calc(obj1, Aircraft.Certification.Aerodynamic_data.CD0.value, ...
                             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value(i), ...
                             Aircraft.Geometry.Wing.AR.value, ...
                             Aircraft.Certification.Aerodynamic_data.e.value, ...
                             Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value, ...
                             Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value = CD_positivestall;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.Attributes.unit = "Non dimensional";

% *** From point C to point D ***

% FROM C TO FG
CD_fromCtofg = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
    CD_fromCtofg(i) = cd_calc(obj1, Aircraft.Certification.Aerodynamic_data.CD0.value, ...
                             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value(i), ...
                             Aircraft.Geometry.Wing.AR.value, ...
                             Aircraft.Certification.Aerodynamic_data.e.value, ...
                             Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value, ...
                             Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtofg.value = CD_fromCtofg;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtofg.Attributes.unit = "Non dimensional";


% FROM FG TO D
CD_fromfgtoD = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value)
    CD_fromfgtoD(i) = cd_calc(obj1, Aircraft.Certification.Aerodynamic_data.CD0.value, ...
                             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value(i), ...
                             Aircraft.Geometry.Wing.AR.value, ...
                             Aircraft.Certification.Aerodynamic_data.e.value, ...
                             Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value, ...
                             Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromfgtoD.value = CD_fromfgtoD;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromfgtoD.Attributes.unit = "Non dimensional";

%% CD CALCULATION - NEGATIVE LOAD FACTOR
% *** From point S to point F ***
CD_negativestall = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value)
    CD_negativestall(i) = cd_calc(obj1, Aircraft.Certification.Aerodynamic_data.CD0.value, ...
                             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value(i), ...
                             Aircraft.Geometry.Wing.AR.value, ...
                             Aircraft.Certification.Aerodynamic_data.e.value, ...
                             Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value, ...
                             Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_negativestall.value = CD_negativestall;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_negativestall.Attributes.unit = "Non dimensional";

% *** From point F to E ***
CD_fromFtoE = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value)
    CD_fromFtoE(i) = cd_calc(obj1, Aircraft.Certification.Aerodynamic_data.CD0.value, ...
                             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value(i), ...
                             Aircraft.Geometry.Wing.AR.value, ...
                             Aircraft.Certification.Aerodynamic_data.e.value, ...
                             Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value, ...
                             Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value = CD_fromFtoE;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.Attributes.unit = "Non dimensional";

%% CD CALCULATION - FROM POINT D TO POINT E

CD_fromDtoE = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value)
    CD_fromDtoE(i) = cd_calc(obj1, Aircraft.Certification.Aerodynamic_data.CD0.value, ...
                             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value(i), ...
                             Aircraft.Geometry.Wing.AR.value, ...
                             Aircraft.Certification.Aerodynamic_data.e.value, ...
                             Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value, ...
                             Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDtoE.value = CD_fromDtoE;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDtoE.Attributes.unit = "Non dimensional";

%% ALPHA CALCULATION - POSITIVE LOAD FACTOR

% Interpolation coefficient from actual aerodynamic data
CL_supp    = str2num(Aircraft.Certification.Aerodynamic_data.CL.value);
alpha_supp = str2num(Aircraft.Certification.Aerodynamic_data.alpha.value);
x          = 0.03*ones(length(CL_supp), 1);
p = polyfit(CL_supp + x, alpha_supp, 2);

% Normal force curve slope in 1/deg
% Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value = 0.0913;
% Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.Attributes.unit = "1/degree";

% Alpha zero lift
Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value = alpha_calc_lin(obj1, ...
                                                                               0.0, ...
                                                                               Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                                                               Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.Attributes.unit = "degree";

% *** From point S to point C ***
% alpha_positivestall = zeros(length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value)
    
    % A complete documentation of the function alpha_calc(...) used here is
    % inside the class file aero_model.m, which can be found inside the
    % 'utilities' folder of this library.      
    % alpha_positivestall(i) = alpha_calc_lin(obj1, ...
    %                                    Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_positivestall.value(i), ...
    %                                    Aircraft.Certification.Aerodynamic_data.CL0.value, ...
    %                                    Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
    % Aircraft.Certification.Aerodynamic_data.CL.value, Aircraft.Certification.Aerodynamic_data.alpha.value,
    alpha_positivestall = alpha_calc(obj1, ...
                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p);
% end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall.value = alpha_positivestall;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall.Attributes.unit = "Degrees";

% *** From point C to point D ***
% alpha_fromCtoD = zeros(length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value)
%     alpha_fromCtoD(i) = alpha_calc_lin(obj1, ...
%                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromCtoD.value(i), ...
%                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);

% FROM C TO fg
      alpha_fromCtofg = alpha_calc(obj1, ...
                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p);
% end

% Alpha from C to fg
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg.value = alpha_fromCtofg;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg.Attributes.unit = "Degrees";

% FROM fg TO D
      alpha_fromfgtoD = alpha_calc(obj1, ...
                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p);
% end

% Alpha from fg to D
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD.value = alpha_fromfgtoD;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD.Attributes.unit = "Degrees";

%% ALPHA CALCULATION - NEGATIVE LOAD FACTOR

% *** From point S to point F *** 
% alpha_negativestall = zeros(length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Negative_stall_speed.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Negative_stall_speed.value)
%     alpha_negativestall(i) = alpha_calc_lin(obj1, ...
%                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_negativestall.value(i), ...
%                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
      alpha_negativestall = alpha_calc(obj1, ...
                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p);

% end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall.value = alpha_negativestall - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall.Attributes.unit = "Degrees";

% *** From point F to point E ***
% alpha_fromFtoE = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromFtoE.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromFtoE.value)
%     alpha_fromFtoE(i) = alpha_calc_lin(obj1, ...
%                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromFtoE.value(i), ...
%                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
      alpha_fromFtoE = alpha_calc(obj1, ...
                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p);
% end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE.value = alpha_fromFtoE - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE.Attributes.unit = "Degrees";

%% ALPHA CALCULATION - FROM POINT D TO POINT E
% alpha_fromDtoE = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromDtoE.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromDtoE.value)
%     alpha_fromDtoE(i) = alpha_calc_lin(obj1, ...
%                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromDtoE.value(i), ...
%                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
      alpha_fromDtoE = alpha_calc(obj1, ...
                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p);
% end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE.value = alpha_fromDtoE - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE.Attributes.unit = "Degrees";

%% ALPHA CONVERSION IN RADIANS
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall_rad.value           = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall_rad.Attributes.unit = "Radians";
% *************************************************************************
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_rad.value                = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_rad.Attributes.unit      = "Radians";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_rad.value                = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_rad.Attributes.unit      = "Radians";
% -------------------------------------------------------------------------
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall_rad.value           = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall_rad.Attributes.unit = "Radians";
% *************************************************************************
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_rad.value                = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_rad.Attributes.unit      = "Radians";
% ---------------------------------------------------------------------------
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_rad.value                = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_rad.Attributes.unit      = "Radians";

%% WING BODY LIFT CALCULATION
% In this section of the code the following formula is applied to evaluate
% the wing-body lift with the formula 
%        
%   L = 0.5*rho*(V^2)*Sw*CL 
%
% using the previously calculated lift coefficients, coming from the 
% various curves of the final envelope diagram.

%   Positive stall speed
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_positivestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value)   
    V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_positivestall.value(i) = (0.5)*(V^2)* ...
        (Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_positivestall.Attributes.unit = "daN";

% *** From point C to point D ***
% FROM C TO FG
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtofg.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value)   
    V = V_fromCtofg(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtofg.value(i) = (0.5)*(V^2)* ...
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtofg.Attributes.unit = "daN";

% FROM FG TO D
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromfgtoD.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value)   
    V = V_fromfgtoD(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromfgtoD.value(i) = (0.5)*(V^2)* ...
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromfgtoD.Attributes.unit = "daN";
% *************************************************************************

% Negative stall speed
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_negativestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value)   
    V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_negativestall.value(i) = (0.5)*(V^2)* ...
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_negativestall.Attributes.unit = "daN";

% *** From point F to point E ***
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value)   
    V = V_fromFtoE(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE.value(i) = (0.5)*(V^2)* ...
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE.Attributes.unit = "daN";
% *************************************************************************

% *** From point D to point E ***
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value)   
    V = V_fromDtoE(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDtoE.value(i) = (0.5)*(V^2)* ... 
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDtoE.Attributes.unit = "daN";

%% PITCHING MOMENT CALCULATIONS - CL WING BODY CONTRIBUTIONS
% Main lifting surface lift pitching moment contribution - POSITIVE LOAD FACTOR
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_positivestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value)
    
    % A complete documentation of the function CLWB_contrib(...) used here 
    % is inside the class file aero_model.m, which can be found inside the
    % 'utilities' folder of this library.      
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_positivestall.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value(i), ...
                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall_rad.value(i), ...
                   Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
                   Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                   Aircraft.Certification.Aerodynamic_data.bcg.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_positivestall.Attributes.unit = "Non dimensional";

% *** From point C to point D ***
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromCtofg.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromCtofg.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value(i), ...
                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_rad.value(i), ...
                   Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
                   Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                   Aircraft.Certification.Aerodynamic_data.bcg.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromCtofg.Attributes.unit = "Non dimensional";

Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromfgtoD.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromfgtoD.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value(i), ...
                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_rad.value(i), ...
                   Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
                   Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                   Aircraft.Certification.Aerodynamic_data.bcg.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromfgtoD.Attributes.unit = "Non dimensional";

% Main lifting surface lift pitching moment contribution - NEGATIVE LOAD FACTOR
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_negativestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_negativestall.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value(i), ...
                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall_rad.value(i), ...
                   Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
                   Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                   Aircraft.Certification.Aerodynamic_data.bcg.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_negativestall.Attributes.unit = "Non dimensional";

% *** From point F to point E ***
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromFtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromFtoE.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value(i), ...
                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_rad.value(i), ...
                   Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
                   Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                   Aircraft.Certification.Aerodynamic_data.bcg.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromFtoE.Attributes.unit = "Non dimensional";

% Main lifting surface lift moment contribution - FROM POINT D TO POINT E 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromDtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromDtoE.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value(i), ...
                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_rad.value(i), ...
                   Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
                   Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                   Aircraft.Certification.Aerodynamic_data.bcg.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromDtoE.Attributes.unit = "Non dimensional";

%% PITCHING MOMENT CALCULATION - CD WING BODY CONTRIBUTIONS
% Main lifting surface drag pitching moment contribution - POSITIVE LOAD FACTORS 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_positivestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value)
    
    % A complete documentation of the function CLWB_contrib(...) used here 
    % is inside the class file aero_model.m, which can be found inside the
    % 'utilities' folder of this library.      
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_positivestall.value(i) = CDWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value(i), ...
                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall_rad.value(i), ...
                   Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
                   Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                   Aircraft.Certification.Aerodynamic_data.bcg.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_positivestall.Attributes.unit = "Non dimensional";

% *** From point C to point D ***

% FROM C TO FG
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromCtofg.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtofg.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtofg.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromCtofg.value(i) = CDWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtofg.value(i), ...
                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_rad.value(i), ...
                   Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
                   Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                   Aircraft.Certification.Aerodynamic_data.bcg.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromCtofg.Attributes.unit = "Non dimensional";

% FROM FG TO D
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromfgtoD.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromfgtoD.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromfgtoD.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromfgtoD.value(i) = CDWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromfgtoD.value(i), ...
                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_rad.value(i), ...
                   Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
                   Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                   Aircraft.Certification.Aerodynamic_data.bcg.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromfgtoD.Attributes.unit = "Non dimensional";

% Main lifting surface drag pitching moment contribution - NEGATIVE LOAD FACTOR
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_negativestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_negativestall.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_negativestall.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_negativestall.value(i) = CDWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_negativestall.value(i), ...
                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall_rad.value(i), ...
                   Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
                   Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                   Aircraft.Certification.Aerodynamic_data.bcg.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_negativestall.Attributes.unit = "Non dimensional";

% *** From point F to point E ***
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromFtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromFtoE.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value(i), ...
                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_rad.value(i), ...
                   Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
                   Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                   Aircraft.Certification.Aerodynamic_data.bcg.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromFtoE.Attributes.unit = "Non dimensional";

% Main lifting surface drag moment contribution - FROM POINT D TO POINT E 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromDtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDtoE.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDtoE.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromDtoE.value(i) = CDWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value(i), ...
                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_rad.value(i), ...
                   Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
                   Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                   Aircraft.Certification.Aerodynamic_data.bcg.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromDtoE.Attributes.unit = "Non dimensional";

%% PITCHING MOMENT CALCULATION - CT CONTRIBUTIONS
% Non - baricentral thrust - POSITIVE LOAD FACTORS DRAG COEFFICIENTS
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_positivestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value)
    
    % A complete documentation of the function CT_contrib(...) used here 
    % is inside the class file aero_model.m, which can be found inside the
    % 'utilities' folder of this library.     
    
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_positivestall.value(i) = CT_contr(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value(i), ...
                   Aircraft.Geometry.Engine.Primary.Thrust_axes.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_positivestall.Attributes.unit = "Non dimensional";

% *** From point C to point D ***

% FROM C TO FG
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromCtofg.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtofg.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtofg.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromCtofg.value(i) = CT_contr(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromCtofg.value(i), ...
                   Aircraft.Geometry.Engine.Primary.Thrust_axes.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromCtofg.Attributes.unit = "Non dimensional";

% FROM FG TO D
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromfgtoD.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromfgtoD.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromfgtoD.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromfgtoD.value(i) = CT_contr(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromfgtoD.value(i), ...
                   Aircraft.Geometry.Engine.Primary.Thrust_axes.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromfgtoD.Attributes.unit = "Non dimensional";

% Non - baricentral thrust - NEGATIVE LOAD FACTORS DRAG COEFFICIENTS
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_negativestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_negativestall.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_negativestall.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_negativestall.value(i) = CT_contr(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_negativestall.value(i), ...
                   Aircraft.Geometry.Engine.Primary.Thrust_axes.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_negativestall.Attributes.unit = "Non dimensional";

% *** From point F to point E ***
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromFtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromFtoE.value(i) = CT_contr(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value(i), ...
                   Aircraft.Geometry.Engine.Primary.Thrust_axes.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromFtoE.Attributes.unit = "Non dimensional";

% Non - baricental thrust - FROM POINT D TO POINT E 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromDtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDtoE.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDtoE.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromDtoE.value(i) = CT_contr(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDtoE.value(i), ...
                   Aircraft.Geometry.Engine.Primary.Thrust_axes.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromDtoE.Attributes.unit = "Non dimensional";

%% PITCHING MOMENT COEFFICIENT OF THE WING BODY CONFIGURATION 
% Wing body pitching moment - POSITIVE LOAD FACTORS
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_positivestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value)
    
    % A complete documentation of the function CM_aboutcg(...) used here 
    % is inside the class file aero_model.m, which can be found inside the
    % 'utilities' folder of this library.     
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_positivestall.value(i) = CM_aboutcg(obj1, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CM0.value, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CMCL.value, ...
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value(i));
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_positivestall.Attributes.unit = "Non dimensional";

% *** From point C to point D ***

% FROM C TO FG
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromCtofg.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromCtofg.value(i) = CM_aboutcg(obj1, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CM0.value, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CMCL.value, ...
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value(i));
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromCtofg.Attributes.unit = "Non dimensional";

% FROM FG TO D
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromfgtoD.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromfgtoD.value(i) = CM_aboutcg(obj1, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CM0.value, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CMCL.value, ...
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value(i));
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromfgtoD.Attributes.unit = "Non dimensional";

% Wing body pitching moment - NEGATIVE LOAD FACTORS
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_negativestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_negativestall.value(i) = CM_aboutcg(obj1, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CM0.value, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CMCL.value, ...
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value(i));
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_negativestall.Attributes.unit = "Non dimensional";

% *** From point F to point E ***
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromFtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromFtoE.value(i) = CM_aboutcg(obj1, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CM0.value, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CMCL.value, ...
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value(i));
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromFtoE.Attributes.unit = "Non dimensional";

% *** From point D to point E ***
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromDtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromDtoE.value(i) = CM_aboutcg(obj1, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CM0.value, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CMCL.value, ...
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value(i));
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromDtoE.Attributes.unit = "Non dimensional";

%% HORIZONTAL TAIL LIFT COEFFICIENT CALCULATION
% POSITIVE LOAD FACTORS
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_positivestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall_rad.value) , 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall_rad.value)
    
    % A complete documentation of the function CL_Tail(...) used here 
    % is inside the class file aero_model.m, which can be found inside the
    % 'utilities' folder of this library.        
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_positivestall.value(i) = CL_Tail(obj1, ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_positivestall.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_positivestall.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_positivestall.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_positivestall.value(i), ...
                                                        Aircraft.Geometry.Horizontal.l.value, ...
                                                        Aircraft.Geometry.Wing.mac.value, ...
                                                        Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ...
                                                        Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_positivestall_rad.value(i));
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_positivestall.Attributes.unit = "Non dimensional";
                                                   
% *** From point C to point D ***

% FROM C TO FG
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromCtofg.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_rad.value) , 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_rad.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromCtofg.value(i) = CL_Tail(obj1, ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromCtofg.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromCtofg.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromCtofg.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromCtofg.value(i), ...
                                                        Aircraft.Geometry.Horizontal.l.value, ...
                                                        Aircraft.Geometry.Wing.mac.value, ...
                                                        Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ...
                                                        Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_rad.value(i));
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromCtofg.Attributes.unit = "Non dimensional";

% FROM FG TO D
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromfgtoD.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_rad.value) , 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_rad.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromfgtoD.value(i) = CL_Tail(obj1, ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromfgtoD.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromfgtoD.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromfgtoD.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromfgtoD.value(i), ...
                                                        Aircraft.Geometry.Horizontal.l.value, ...
                                                        Aircraft.Geometry.Wing.mac.value, ...
                                                        Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ...
                                                        Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_rad.value(i));
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromfgtoD.Attributes.unit = "Non dimensional";

% NEGATIVE LOAD FACTORS
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_negativestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall_rad.value) , 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall_rad.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_negativestall.value(i) = CL_Tail(obj1, ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_negativestall.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_negativestall.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_negativestall.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_negativestall.value(i), ...
                                                        Aircraft.Geometry.Horizontal.l.value, ...
                                                        Aircraft.Geometry.Wing.mac.value, ...
                                                        Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ...
                                                        Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negativestall_rad.value(i));
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_negativestall.Attributes.unit = "Non dimensional";

% *** From point F to point E ***
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromFtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_rad.value) , 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_rad.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromFtoE.value(i) = CL_Tail(obj1, ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromFtoE.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromFtoE.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromFtoE.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromFtoE.value(i), ...
                                                        Aircraft.Geometry.Horizontal.l.value, ...
                                                        Aircraft.Geometry.Wing.mac.value, ...
                                                        Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ...
                                                        Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_rad.value(i));
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromFtoE.Attributes.unit = "Non dimensional";

% *** From point D to point E ***
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromDtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_rad.value) , 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_rad.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromDtoE.value(i) = CL_Tail(obj1, ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_fromDtoE.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromDtoE.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_fromDtoE.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_fromDtoE.value(i), ...
                                                        Aircraft.Geometry.Horizontal.l.value, ...
                                                        Aircraft.Geometry.Wing.mac.value, ...
                                                        Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ...
                                                        Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_rad.value(i));
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromDtoE.Attributes.unit = "Non dimensional";

%% HORIZONTAL TAIL DIMENSIONAL AIRLOADS
% Positive stall speed
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_positivestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_positivestall.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_positivestall.value)   
    V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_positivestall.value(i) = (0.5)*(V^2)* ...
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_positivestall.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_positivestall.Attributes.unit = "daN";

% *** From point C to point D ***
% FROM C TO FG
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromCtofg.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromCtofg.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromCtofg.value)   
    V = V_fromCtofg(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromCtofg.value(i) = (0.5)*(V^2)* ... 
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromCtofg.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromCtofg.Attributes.unit = "daN";
% FROM FG TO D
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromfgtoD.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromfgtoD.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromfgtoD.value)   
    V = V_fromfgtoD(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromfgtoD.value(i) = (0.5)*(V^2)* ... 
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromfgtoD.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromfgtoD.Attributes.unit = "daN";
% *************************************************************************

% Negative stall speed
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_negativestall.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_negativestall.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_negativestall.value)   
    V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_negativestall.value(i) = (0.5)*(V^2)* ...
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_negativestall.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_negativestall.Attributes.unit = "daN";

% *** From point F to point E ***
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromFtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromFtoE.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value)   
    V = V_fromFtoE(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromFtoE.value(i) = (0.5)*(V^2)* ... 
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromFtoE.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromFtoE.Attributes.unit = "daN";
% *************************************************************************

% *** From point D to point E ***
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromDtoE.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromDtoE.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromDtoE.value)   
    V = V_fromDtoE(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromDtoE.value(i) = (0.5)*(V^2)* ...
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromDtoE.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromDtoE.Attributes.unit = "daN";

%% TAIL BALANCING LOADS - N = 1.0 

Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value), 1);
VS = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value(1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value(1) = VS;
VD = V_fromfgtoD(end);
for i = 2:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value)
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value(i-1) + (VD - VS)*(1/length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value));
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.Attributes.unit = "m/s";

% DYNAMIC PRESSURE
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.value = 0.5*Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value.^2);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.Attributes.unit = "Pa";

% WING BODY LIFT COEFFICIENT
Dyn_press_force = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.q_unit_load_factor.value*Aircraft.Geometry.Wing.S.value;
Aircraft_weight = ones(length(Dyn_press_force),1)*Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value = zeros(length(Dyn_press_force), 1);
for i = 1:length(Dyn_press_force)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value(i) = ((Aircraft_weight(i))*(1/Dyn_press_force(i)));
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.Attributes.unit = "Non dimensional";

% DRAG LIFT COEFFICIENT
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value(i) = cd_calc(obj1, Aircraft.Certification.Aerodynamic_data.CD0.value, ...
                             Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value(i), ...
                             Aircraft.Geometry.Wing.AR.value, ...
                             Aircraft.Certification.Aerodynamic_data.e.value, ...
                             Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value, ...
                             Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.Attributes.unit = "Non dimensional";

% ALPHA CALCULATION 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromDtoE.value), 1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.value = alpha_calc(obj1, ...
                                                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value, ...
                                                   Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                                   Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                                   Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                                   p);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.Attributes.unit = "degrees";

Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.Attributes.unit = "radians";

% WING BODY LIFT 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value)   
    V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.value(i) = (0.5)*(V^2)* ...
        (Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.Attributes.unit = "daN";

obj1 = aero_model;
% EVALUATION OF TAIL LOADS CONTRIBUTION
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value),1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_unit_load_factor.value)
%     y  = Aircraft.Certification.Aerodynamic_data.XAC_nondim.value - Aircraft.Certification.Aerodynamic_data.XCG_nondim.value;
%     x  = Aircraft.Certification.Aerodynamic_data.bcg.value/Aircraft.Geometry.Wing.mac.value;
%     t1 = cos(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value(i));
%     t2 = sin(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value(i));
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_unit_load_factor.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value(i)*t1*y - Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value(i)*2*(x);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_unit_load_factor.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value(i), ...
                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value(i), ...
                   Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
                   Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                   Aircraft.Certification.Aerodynamic_data.bcg.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_unit_load_factor.Attributes.unit = "Non dimensional";

Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_fromFtoE.value(i) = CLWB_contrib(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value(i), ...
                   Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value(i), ...
                   Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ... 
                   Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                   Aircraft.Certification.Aerodynamic_data.bcg.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_unit_load_factor.Attributes.unit = "Non dimensional";

% *** From point C to point D ***
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_unit_load_factor.value(i) = CT_contr(obj1, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_unit_load_factor.value(i), ...
                   Aircraft.Geometry.Engine.Primary.Thrust_axes.value, ...
                   Aircraft.Geometry.Wing.mac.value);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_unit_load_factor.Attributes.unit = "Non dimensional";

Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_unit_load_factor.value(i) = CM_aboutcg(obj1, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CM0.value, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.CMCL.value, ...
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLWB_unit_load_factor.value(i));
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_unit_load_factor.Attributes.unit = "Non dimensional";

% TAIL LIFT COEFFICIENT 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value) , 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_unit_load_factor.value(i) = CL_Tail(obj1, ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLwb_pitch_contrib_unit_load_factor.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CDwb_pitch_contrib_unit_load_factor.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CTwb_pitch_contrib_unit_load_factor.value(i), ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CMwb_pitch_contrib_unit_load_factor.value(i), ...
                                                        Aircraft.Geometry.Horizontal.l.value, ...
                                                        Aircraft.Geometry.Wing.mac.value, ...
                                                        Aircraft.Certification.Aerodynamic_data.XAC_nondim.value, ...
                                                        Aircraft.Certification.Aerodynamic_data.XCG_nondim.value, ...
                                                        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_unit_load_factor_rad.value(i));
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_unit_load_factor.Attributes.unit = "Non dimensional";

% *** TAIL LIFT AT N = 1.0 ***
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_unit_load_factor.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_unit_load_factor.value)   
    V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_unit_load_factor.value(i) = (0.5)*(V^2)* ... 
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_unit_load_factor.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_unit_load_factor.Attributes.unit = "daN";

% *** WING LIFT AT N = 1.0 ***
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Wing_Lift_unit_load_factor.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_unit_load_factor.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_unit_load_factor.value) 
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Wing_Lift_unit_load_factor.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_unit_load_factor.value(i) - Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_unit_load_factor.value(i);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Wing_Lift_unit_load_factor.Attributes.unit = "daN";

%% BALANCING LOADS DIAGRAM 
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
% -------------------------------------------------------------------------
% To plot the airloads acting on the horizontal surface (which are
% necessary to evaluate the whole aircraft Lift, which must be used to size
% the wing structure) the following function will be used: 
%
% HTailAirloadsDiagram = Balancing_loads(HT_Lift_posstall, VSpos, ...
%                                        HT_Lift_negstall, VSneg, ... 
%                                        HT_Lift_fromCtoD, V_fromCtoD, ...
%                                        HT_Lift_fromDtoE, V_fromDtoE, ...
%                                        HT_Lift_fromFtoE, V_fromFtoE, ...
%                                        HTail_Lift_unit_load_factor.value, V_unit_load_factor, ... 
%                                        Reg, Aircraft_name)
%
% A complete documentation of this function is included inside the file
% csvla.m

disp(" ")
disp(" ++++ FIGURE 7 - HORIZONTAL EMPENNAGE AIRLOADS ++++ ");
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTailAirloadsDiagram.value = Balancing_loads(obj, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_positivestall.value, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_negativestall.value, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value, ... 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromCtofg.value, ...
                V_fromCtofg, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromfgtoD.value, ...
                V_fromfgtoD, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromDtoE.value, ...
                V_fromDtoE, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromFtoE.value, ...
                V_fromFtoE, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_unit_load_factor.value, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.value, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_Design_manoeuvring_speed_VG.value, ...
                Aircraft.Certification.Regulation.value, ... 
                Aircraft.Certification.Aircraft_Name.value);
pause(1);
% Saving figures inside correct folder
fprintf('Saving Balancingloads.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile Balancingloads.pdf Output
movefile Balancingloads.png Output

%% MAIN WING LOADS DIAGRAM 
%   In this section we have to take into account the lift coefficient
%   produced by the horizontal tail. The global lift coefficient will be
%   the algebric sum of the Wing - Body lift coefficient and the horizontal
%   tail lift coefficient.

% LIFT COEFFICIENT CALCULATIONS
% POSITIVE STALL
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_posstall_new.value = ...
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_positivestall.value - Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_positivestall.value;

Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_posstall_new.Attributes.unit = "Non dimensional";

% FROM C TO D
% FROM C TO FG
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg_new.value = ...
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg.value - Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromCtofg.value;

Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg_new.Attributes.unit = "Non dimensional";

% FROM FG TO D
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD_new.value = ...
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD.value - Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromfgtoD.value;

Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD_new.Attributes.unit = "Non dimensional";

% NEGATIVE STALL
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negstall_new.value = ...
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negativestall.value - Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_negativestall.value;

Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negstall_new.Attributes.unit = "Non dimensional";

% FROM F TO E 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value = ...
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE.value - Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromFtoE.value;

Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.Attributes.unit = "Non dimensional";

% FROM D TO E
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE_new.value = ...
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value - Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_fromDtoE.value;

Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE_new.Attributes.unit = "Non dimensional";

%% LIFT CALCULATIONS WITH NEW VALUES OF THE LIFT COEFFICIENT
% Positive stall speed
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_posstall_new.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_posstall_new.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_posstall_new.value)   
    V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_posstall_new.value(i) = (0.5)*(V^2)* ...
        (Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_posstall_new.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_posstall_new.Attributes.unit = "daN";

% *** From point C to point D ***

% FROM C TO FG
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtofg_new.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg_new.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg_new.value)   
    V = V_fromCtofg(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtofg_new.value(i) = (0.5)*(V^2)* ...
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg_new.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtofg_new.Attributes.unit = "daN";

% FROM FG TO Ds
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromfgtoD_new.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD_new.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD_new.value)   
    V = V_fromfgtoD(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromfgtoD_new.value(i) = (0.5)*(V^2)* ...
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD_new.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromfgtoD_new.Attributes.unit = "daN";
% *************************************************************************

% Negative stall speed
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_negstall_new.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negstall_new.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negstall_new.value)   
    V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_negstall_new.value(i) = (0.5)*(V^2)* ...
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negstall_new.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_negstall_new.Attributes.unit = "daN";

% *** From point F to point E ***
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE_new.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value)   
    V = V_fromFtoE(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE_new.value(i) = (0.5)*(V^2)* ...
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE_new.Attributes.unit = "daN";
% *************************************************************************

% *** From point D to point E ***
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDtoE_new.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE_new.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE_new.value)   
    V = V_fromDtoE(i);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDtoE_new.value(i) = (0.5)*(V^2)* ... 
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE_new.value(i))*(1E-1);   
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDtoE_new.Attributes.unit = "daN";

%% MAIN WING AIRLOADS DIAGRAM
% To plot the airloads acting on the main lifting surface (evaluated as the
% sum of the wing body Lift plus the horizontal tail lift at equilibrium)
% the following function will be used:
%
% Mainwing_loads(obj, WING_Lift_posstall, VSpos, ...
%                     WING_Lift_negstall, VSneg, ... 
%                     WING_Lift_fromCtoD, V_fromCtoD, ...
%                     WING_Lift_fromDtoE, V_fromDtoE, ...
%                     WING_Lift_fromFtoE, V_fromFtoE, ...
%                     Reg, Aircraft_name)
%
% A complete documentation of this function is included inside the file
% csvla.m

disp(" ")
disp(" ++++ FIGURE 8 - MAIN WING AIRLOADS ++++ ");                                
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WingAirloadsDiagram.value = Mainwing_loads(obj, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_posstall_new.value, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_negstall_new.value, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value, ... 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromCtofg_new.value, ...
                V_fromCtofg, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromfgtoD_new.value, ...
                V_fromfgtoD, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDtoE_new.value, ...
                V_fromDtoE, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE_new.value, ...
                V_fromFtoE, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Wing_Lift_unit_load_factor.value, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.value, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_Design_manoeuvring_speed_VG.value, ...
                Aircraft.Certification.Regulation.value, ... 
                Aircraft.Certification.Aircraft_Name.value);
% Saving figures inside correct folder
fprintf('Saving Wingairloads.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile Wingairloads.pdf Output
movefile Wingairloads.png Output

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

%% ALPHA CALCULATION WITH NEW LIFT COEFFICIENT - POSITIVE LOAD FACTOR

% *** From point S to point C ***
% alpha_posstall_new = zeros(length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value)
%     alpha_posstall_new(i) = alpha_calc_lin(obj1, ...
%                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_posstall_new.value(i), ...
%                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
      alpha_posstall_new = alpha_calc(obj1, ...
                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_posstall_new.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p);

% end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_posstall_new.value = alpha_posstall_new;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_posstall_new.Attributes.unit = "Degrees";

% *** From point C to point D ***
% alpha_fromCtoD_new = zeros(length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value)
%     alpha_fromCtoD_new(i) = alpha_calc_lin(obj1, ...
%                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromCtoD_new.value(i), ...
%                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);

% FROM C TO FG
      alpha_fromCtofg_new = alpha_calc(obj1, ...
                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromCtofg_new.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_new.value = alpha_fromCtofg_new;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_new.Attributes.unit = "Degrees";

% FROM FG TO D 
      alpha_fromfgtoD_new = alpha_calc(obj1, ...
                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromfgtoD_new.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_new.value = alpha_fromfgtoD_new;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_new.Attributes.unit = "Degrees";

%% ALPHA CALCULATION WITH NEW LIFT COEFFICIENT - NEGATIVE LOAD FACTOR

% *** From point S to point F *** 
% alpha_negstall_new = zeros(length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Negative_stall_speed.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Negative_stall_speed.value)
%     alpha_negstall_new(i) = alpha_calc_lin(obj1, ...
%                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_negstall_new.value(i), ...
%                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
      alpha_negstall_new = alpha_calc(obj1, ...
                                       abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_negstall_new.value), ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p);
                                   
% end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negstall_new.value = alpha_negstall_new - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negstall_new.Attributes.unit = "Degrees";

% *** From point F to point E ***
% alpha_fromFtoE_new = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromFtoE_new.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromFtoE_new.value)
%     alpha_fromFtoE_new(i) = alpha_calc_lin(obj1, ...
%                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromFtoE_new.value(i), ...
%                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
      alpha_fromFtoE_new = alpha_calc(obj1, ...
                                       abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value), ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p);

% end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_new.value = alpha_fromFtoE_new - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_new.Attributes.unit = "Degrees";

%% ALPHA CALCULATION WITH NEW LIFT COEFFICIENT - FROM POINT D TO POINT E
% alpha_fromDtoE_new = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromDtoE_new.value), 1);
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromDtoE_new.value)
%     alpha_fromDtoE_new(i) = alpha_calc_lin(obj1, ...
%                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Flightloads.Balancingloads.CL_fromDtoE_new.value(i), ...
%                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
      alpha_fromDtoE_new = alpha_calc(obj1, ...
                                       abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE.value), ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p);

% end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_new.value = alpha_fromDtoE_new - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_new.Attributes.unit = "Degrees";

%% ALPHA_NEW CONVERSION IN RADIANS
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_posstall_new_rad.value           = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_posstall_new.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_posstall_new_rad.Attributes.unit = "Radians";
% *************************************************************************
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_new_rad.value           = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_new.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromCtofg_new_rad.Attributes.unit = "Radians";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_new_rad.value           = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_new.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromfgtoD_new_rad.Attributes.unit = "Radians";
% -------------------------------------------------------------------------
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negstall_new_rad.value           = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negstall_new.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_negstall_new_rad.Attributes.unit = "Radians";
% *************************************************************************
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_new_rad.value           = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_new.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_new_rad.Attributes.unit = "Radians";
% ---------------------------------------------------------------------------
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_new_rad.value           = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_new.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_new_rad.Attributes.unit = "Radians";

%% WING GEOMETRY FOR OPEN VSP

% Define all the geometrical panels for the Open VSP calculations; in
% general the wing will be subdivided in a certain number of panels (three
% in our actual case) with different values of the chords and different
% span.

% Uncomment this section for a general, trapezoidal wing
% Aircraft.Geometry.Wing.Kinks.Croot1.value = Aircraft.Geometry.Wing.croot.value;
% Aircraft.Geometry.Wing.Kinks.Croot1.Attributes.unit = "m";
% Aircraft.Geometry.Wing.Kinks.Ctip3.value = Aircraft.Geometry.Wing.ctip.value;
% Aircraft.Geometry.Wing.Kinks.Ctip3.Attributes.unit = "m";
% Aircraft.Geometry.Wing.Kinks.Croot2.value = NaN;
% Aircraft.Geometry.Wing.Kinks.Croot2.Attributes.unit = "m";
% Aircraft.Geometry.Wing.Kinks.Ctip1.value = NaN;
% Aircraft.Geometry.Wing.Kinks.Ctip1.Attributes.unit = "m";
% Aircraft.Geometry.Wing.Kinks.Croot3.value = NaN;
% Aircraft.Geometry.Wing.Kinks.Croot3.Attributes.unit = "m";
% Aircraft.Geometry.Wing.Kinks.Ctip2.value = NaN;
% Aircraft.Geometry.Wing.Kinks.Ctip2.Attributes.unit = "m";
% Aircraft.Geometry.Wing.Kinks.b1.value = Aircraft.Geometry.Wing.b*(1/3);
% Aircraft.Geometry.Wing.Kinks.b1.Attributes.unit = "m"; 
% Aircraft.Geometry.Wing.Kinks.b2.value = Aircraft.Geometry.Wing.b*(1/3);
% Aircraft.Geometry.Wing.Kinks.b2.Attributes.unit = "m"; 
% Aircraft.Geometry.Wing.Kinks.b3.value = Aircraft.Geometry.Wing.b*(1/3);
% Aircraft.Geometry.Wing.Kinks.b3.Attributes.unit = "m"; 

% Chord and span for the panels in our actual case 
Aircraft.Geometry.Kinks.Croot1.value = Aircraft.Geometry.Wing.croot.value;
Aircraft.Geometry.Kinks.Croot1.Attributes.unit = "m";
Aircraft.Geometry.Kinks.Croot2.value = Aircraft.Geometry.Wing.croot.value;
Aircraft.Geometry.Kinks.Croot2.Attributes.unit = "m";
Aircraft.Geometry.Kinks.Ctip1.value = Aircraft.Geometry.Kinks.Croot2.value;
Aircraft.Geometry.Kinks.Ctip1.Attributes.unit = "m";
Aircraft.Geometry.Kinks.Croot3.value = Aircraft.Geometry.Wing.croot.value;
Aircraft.Geometry.Kinks.Croot3.Attributes.unit = "m";
Aircraft.Geometry.Kinks.Ctip2.value = Aircraft.Geometry.Kinks.Croot3.value;
Aircraft.Geometry.Kinks.Ctip2.Attributes.unit = "m";
Aircraft.Geometry.Kinks.Ctip3.value = Aircraft.Geometry.Wing.ctip.value;
Aircraft.Geometry.Kinks.Ctip3.Attributes.unit = "m";
Aircraft.Geometry.Kinks.b1.value = Aircraft.Geometry.Wing.b.value*(1/3);
Aircraft.Geometry.Kinks.b1.Attributes.unit = "m"; 
Aircraft.Geometry.Kinks.b2.value = Aircraft.Geometry.Wing.b.value*(1/3);
Aircraft.Geometry.Kinks.b2.Attributes.unit = "m"; 
Aircraft.Geometry.Kinks.b3.value = Aircraft.Geometry.Wing.b.value*(1/3);
Aircraft.Geometry.Kinks.b3.Attributes.unit = "m"; 
