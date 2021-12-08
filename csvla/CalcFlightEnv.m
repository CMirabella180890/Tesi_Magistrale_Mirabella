close all; 
%% STARTING A TEST CASE FOR FINAL ENVELOPE DIAGRAM

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

%% STANDARD ATMOSPHERE

h0 = Aircraft.Certification.ISA_Condition.Sea_level.Altitude.value;
[T0, a0, p0, rho0] = atmosisa(h0);
fprintf("--------------------------------------");
fprintf('\n');
fprintf("--------------------------------------");
fprintf('\n');
fprintf("### Standard atmosphere - Sea Level ###");
fprintf('\n');
fprintf("--------------------------------------");
fprintf('\n');
fprintf('Temperature [K]: ');
fprintf('%5.4f%', T0);
fprintf('\n');
fprintf('Speed of sound [m/s]: ');
fprintf('%5.4f%', a0);
fprintf('\n');
fprintf('Pressure [Pa]: ');
fprintf('%5.4f%', p0);
fprintf('\n'); 
fprintf('Density [kg/m^3]: ');
fprintf('%5.4f%', rho0);
fprintf('\n');
fprintf("--------------------------------------"); 
Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value = rho0;
Aircraft.Certification.ISA_Condition.Sea_Level.rho0.Attributes.unit = 'kg/m^3';
Aircraft.Certification.ISA_Condition.Sea_Level.Altitude.value = h0;
Aircraft.Certification.ISA_Condition.Sea_Level.Altitude.Attributes.unit = 'm';
Aircraft.Certification.ISA_Condition.Sea_Level.T0.value = T0; 
Aircraft.Certification.ISA_Condition.Sea_Level.T0.Attributes.unit = 'K'; 
Aircraft.Certification.ISA_Condition.Sea_Level.p0.value = p0; 
Aircraft.Certification.ISA_Condition.Sea_Level.p0.Attributes.unit = 'Pa';
Aircraft.Certification.ISA_Condition.Sea_Level.Speed_of_sound0.value = a0; 
Aircraft.Certification.ISA_Condition.Sea_Level.Speed_of_sound0.Attributes.unit = 'm/s';
% -------------------------------------------------------------------------

h_operative = Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.value;
[T_operative, a_operative, p_operative, rho_operative] = atmosisa(h_operative);
fprintf("--------------------------------------");
fprintf('\n');
fprintf("--------------------------------------");
fprintf('\n');
fprintf("### Standard atmosphere - Operative Ceiling ###");
fprintf('\n');
fprintf("--------------------------------------");
fprintf('\n');
fprintf('Temperature [K]: ');
fprintf('%5.4f%', T_operative);
fprintf('\n');
fprintf('Speed of sound [m/s]: ');
fprintf('%5.4f%', a_operative);
fprintf('\n');
fprintf('Pressure [Pa]: ');
fprintf('%5.4f%', p_operative);
fprintf('\n'); 
fprintf('Density [kg/m^3]: ');
fprintf('%5.4f%', rho_operative);
fprintf('\n');
fprintf("--------------------------------------"); 
Aircraft.Certification.ISA_Condition.Operative_ceiling.rho0.value = rho_operative;
Aircraft.Certification.ISA_Condition.Operative_ceiling.rho0.Attributes.unit = 'kg/m^3';
Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.value = h_operative;
Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.Attributes.unit = 'm';
Aircraft.Certification.ISA_Condition.Operative_ceiling.T0.value = T_operative; 
Aircraft.Certification.ISA_Condition.Operative_ceiling.T0.Attributes.unit = 'K'; 
Aircraft.Certification.ISA_Condition.Operative_ceiling.p0.value = p_operative; 
Aircraft.Certification.ISA_Condition.Operative_ceiling.p0.Attributes.unit = 'Pa';
Aircraft.Certification.ISA_Condition.Operative_ceiling.Speed_of_sound0.value = a_operative; 
Aircraft.Certification.ISA_Condition.Operative_ceiling.Speed_of_sound0.Attributes.unit = 'm/s';
% -------------------------------------------------------------------------

h_theoretical = Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.value;
[T_theoretical, a_theoretical, p_theoretical, rho_theoretical] = atmosisa(h_theoretical);
fprintf("--------------------------------------");
fprintf('\n');
fprintf("--------------------------------------");
fprintf('\n');
fprintf("### Standard atmosphere - Theoretical Ceiling ###");
fprintf('\n');
fprintf("--------------------------------------");
fprintf('\n');
fprintf('Temperature [K]: ');
fprintf('%5.4f%', T_theoretical);
fprintf('\n');
fprintf('Speed of sound [m/s]: ');
fprintf('%5.4f%', a_theoretical);
fprintf('\n');
fprintf('Pressure [Pa]: ');
fprintf('%5.4f%', p_theoretical);
fprintf('\n'); 
fprintf('Density [kg/m^3]: ');
fprintf('%5.4f%', rho_theoretical);
fprintf('\n');
fprintf("--------------------------------------"); 
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.rho0.value = rho_theoretical;
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.rho0.Attributes.unit = 'kg/m^3';
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.value = h_theoretical;
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.Attributes.unit = 'm';
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.T0.value = T_theoretical; 
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.T0.Attributes.unit = 'K'; 
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.p0.value = p_theoretical; 
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.p0.Attributes.unit = 'Pa';
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Speed_of_sound0.value = a_theoretical; 
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Speed_of_sound0.Attributes.unit = 'm/s';

%% NUMBER OF ELEMENTS
numb = 1e3;

% REGULATION APPLIED
Reg = Aircraft.Certification.Regulation.value;

%% WING LOADING DEFINITION
% x = calcn(obj, nmax) - from csvla.m
% This function defines a vector with load factor values between two pre-
% scribed limits. Check the class file csvla.m to have a complete
% documentation.
g    = Aircraft.Constants.g.value;
S1   = 0.0;
S    = Aircraft.Geometry.Wing.S.value + S1;
x    = 90;
Mass = Aircraft.Weight.I_Level.W_maxTakeOff.value+x;
Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value = (Mass*g)/(S);
Aircraft.Certification.Performance.I_Level.Wing_loading_SI.Attributes.unit = "Pa";

% LOCAL VARIABLE WITH WING LOADING
WS = Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value;

% LOCAL VARIABLE WITH OTHER USEFUL QUANTITIES
rho0 = Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value;

% MEAN AERODYNAMIC CHORD 
MAC = Aircraft.Geometry.Wing.mac.value;

% LIFT/NORMAL FORCE SLOPE COEFFICIENT
CLalfa = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;

% LOCALLY DEFINED GUST SPEED 
Ude_cruise = 15.24;
Ude_dive   = 7.62;

% LOCAL VARIABLE WITH MAX LIFT COEFFICIENT
CLMAX_clean          = Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value;
CLMAX_clean_inverted = Aircraft.Certification.Aerodynamic_data.Min_Lift_Coefficient.value;

% LOAD FACTOR VALUES
% Positive load factor values
nmax = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value = calcn(obj, nmax);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.Attributes.unit = "g's";
npos = Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value;
% Negative load factor values
nmin = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.value = calcn(obj, nmin);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.Attributes.unit = "g's";
nneg = Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.value;

%% STALL SPEED CALCULATION
% x = calcvs(obj, rho, WingLoading, MaxLiftCoeff, PositiveLoadFactors)
% This function defines a vector with stall airspeed for the chosen
% aircraft, within the precribed range of load factors. Check the
% class file csvla.m to have a complete documentation.
% Positive stall speed 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_VS.value = calcvs(obj, rho0, ...   % Standard atmosphere density
                                                                                 WS, ...          % Wing Loading in SI units 
                                                                                 CLMAX_clean, ... % Maximum Lift coefficient
                                                                                 npos);           % A vector of load factors
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_VS.Attributes.unit = "m/s"; 
VSpos = Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_VS.value;
% Negative stall speed 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_VS.value = calcvs(obj, rho0, ...             % Standard atmosphere density
                                                                                  WS, ...                   % Wing Loading in SI units 
                                                                                  CLMAX_clean_inverted, ... % Maximum Lift coefficient
                                                                                  nneg);                    % A vector of load factors
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_VS.Attributes.unit = "m/s"; 
VSneg = Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_VS.value;

%% CALCULATION OF THE CRUISE SPEED 
% x = calcvc(obj, WingLoading, MaxContinuousPowerSpeedVH)
% This function identifies (following CS-VLA airworthiness reg.)
% maximum cruise speed (Point C) for flight envelope calculations. 
% To have a complete documentation check the class file csvla.m
% VH design speed for max continous power: this airspeed is not available
% but must be known. From CS - VLA Airworthiness rules
Aircraft.Certification.Regulation.SubpartC.Flightloads.Cruise_Speed_VC.value = calcvc(obj, WS); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Cruise_Speed_VC.Attributes.unit = "m/s";
Aircraft.Certification.Regulation.SubpartC.Flightloads.nC.value = nmax; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.nC.Attributes.unit = "g";
VC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Cruise_Speed_VC.value;
nC = nmax;

%% CALCULATION OF THE DIVE SPEED 
% x = calcvd(obj, MinDesignCruiseSpeed, CruiseSpeedVC)
% This function identifies (following CS-VLA airworthiness reg.)
% the maximum dive speed (Point D) for flight envelope
% calculations. To have a complete documentation check the class
% file csvla.m
Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value = calcvd(obj, Aircraft.Certification.Regulation.SubpartC.Flightloads.Min_Design_Cruise_Speed.value, ... % Min design cruise speed 
                                                                                         VC);                                                                                      % Cruise speed from previous calculations
Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.Attributes.unit = "m/s";
Aircraft.Certification.Regulation.SubpartC.Flightloads.nD.value = nmax; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.nD.Attributes.unit = "g";
VD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value;
nD = nmax;

% INVERTED FLIGHT DIVE SPEED 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VE.value = VD; % Speed at points E and D are equal 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VE.Attributes.unit = "m/s";
Aircraft.Certification.Regulation.SubpartC.Flightloads.nE.value = nmin; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.nE.Attributes.unit = "g";
VE = VD;
nE = nmin;

%% Point S definition 
VS = Vstall(WS, rho0, CLMAX_clean, 1.0); 
nS = 1.0;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Stall_speed_VS.value = VS;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Stall_speed_VS.Attributes.unit = "m/s";

% FLIGHT ENVELOPE STARTING SECTION 
n_from1toS = linspace(0.0, nS, numb);
V_from1toS = VS*ones(numb, 1);

%% Point S_inverted definition 
nS_inv = -1.0;
VS_inv = Vstall(WS, rho0, abs(CLMAX_clean_inverted), abs(nS_inv)); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_Stall_speed_VS.value = VS_inv;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_Stall_speed_VS.Attributes.unit = "m/s";

% FLIGHT ENVELOPE STARTING SECTION 
n_from1toS_inv = linspace(0.0, nS_inv, numb);
V_from1toS_inv = VS_inv*ones(numb, 1);

%% Point A definition
% Assign speed at Point A (Maneuver point) equals to the maximum
% permissible positive load factor stall speed. 

% POSITIVE STALL SPEED AT POINT VA
VA = Vstall(WS, rho0, CLMAX_clean, nmax);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.value = VA;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.Attributes.unit = "m/s";
Aircraft.Certification.Regulation.SubpartC.Flightloads.nA.value = nmax; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.nA.Attributes.unit = "g";
nA = nmax;

% POSITIVE STALL SPEED - FLIGHT ENVELOPE
n_fromStoA = linspace(nS, nA, numb)';
V_fromStoA = Vstall(WS, rho0, CLMAX_clean, n_fromStoA);

%% Point G definition 
% NEGATIVE STALL SPEED AT POINT VG 
VG = Vstall(WS, rho0, abs(CLMAX_clean_inverted), abs(nmin));
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_Design_manoeuvring_speed_VG.value = VG;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_Design_manoeuvring_speed_VG.Attributes.unit = "m/s";
Aircraft.Certification.Regulation.SubpartC.Flightloads.nG.value = nmin; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.nG.Attributes.unit = "g";
nG = nmin;

% NEGATIVE STALL SPEED - FLIGHT ENVELOPE
n_fromStoG = linspace(nS_inv, nG, numb)';
V_fromStoG = Vstall(WS, rho0, abs(CLMAX_clean_inverted), abs(n_fromStoG));

%% Point C definition 
V_fromAtoC = linspace(VA, VC, numb)'; 
n_fromAtoC = nmax*ones(numb, 1); 
nC         = nmax;

%% Point D definition 
V_fromCtoD = linspace(VC, VD, numb)';
n_fromCtoD = nmax*ones(numb, 1); 

%% Point E definition 
V_fromGtoE = linspace(VG, VE, numb)';
n_fromGtoE = nmin*ones(numb, 1);

%% Flight envelope limit 
V_fromDtoE = VD*ones(numb, 1);
n_fromDtoE = linspace(nD, nE, numb)';

%% FLIGHT ENVELOPE

%% INPUT TRACKING - VN DIAGRAM 
% A possible way to track inputs for the various data will be provided
% inside the .txt file used as a log for the program.

disp(" ")
disp(" ++++ INPUT TO V - N DIAGRAM ++++");
% Input to the flight envelope 
Data1 = [  VS, ...              % Straight flight stall airspeed 
           VS_inv, ...          % Inverted flight stall airspeed
           VA, ...              % Design manoeuvring speed
           VG];                 % VG = VD on the negative side of V - n diagram
disp(" ++++++++++ DATA USED TO PLOT V - N DIAGRAM ++++++++++ ")
format = ' %6.6f          %6.6f          %6.6f          %6.6f\n';
label  = ' VS+                VS-                 VA                VG\n';
fprintf(label);
fprintf(format, Data1.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

disp(" ")
% Input to the flight envelope
Data1 = [  nmax, ...            % Max positive value of load factors
           nmin, ...            % Min (negative) value of load factors
           VD, ...              % Max dive speed from V - n diagram
           VE];                 % VG = VD on the negative side of V - n diagram
disp(" ++++++++++ DATA USED TO PLOT V - N DIAGRAM ++++++++++ ")
format = ' %6.6f          %6.6f          %6.6f          %6.6f\n';
label  = ' nmax                nmin                 VD                VE\n';
fprintf(label);
fprintf(format, Data1.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

% AIRCRAFT STRUCT VARIABLE FILLING
% POINT S
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS.VS.value = VS; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS.VS.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS.nS.value = nS; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS.nS.Attributes.unit = "g's"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS.CL_S.value = CLmax_func(rho0, S, VS, WS, nS);
CL_S = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS.CL_S.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS.CL_S.Attribute.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS.qS.value = 0.5*rho0*VS^2;
qS = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS.qS.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS.qS.Attributes.unit = "Pa"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS.LS.value = CL_S*qS*S;

% POINT A
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointA.VA.value = VA; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointA.VA.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointA.nA.value = nA; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointA.nA.Attributes.unit = "g's"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointA.CL_A.value = CLmax_func(rho0, S, VA, WS, nA);
CL_A = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointA.CL_A.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointA.CL_A.Attribute.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointA.qA.value = 0.5*rho0*VA^2;
qA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointA.qA.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointA.qA.Attributes.unit = "Pa"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointA.LA.value = CL_A*qA*S;

% POINT C
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointC.VC.value = VC; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointC.VC.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointC.nC.value = nC; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointC.nC.Attributes.unit = "g's"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointC.CL_C.value = CLmax_func(rho0, S, VC, WS, nC);
CL_C = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointC.CL_C.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointC.CL_C.Attribute.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointC.qC.value = 0.5*rho0*VC^2;
qC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointC.qC.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointC.qC.Attributes.unit = "Pa"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointC.LC.value = CL_C*qC*S;

% POINT D
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointD.VD.value = VD; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointD.VD.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointD.nD.value = nD; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointD.nD.Attributes.unit = "g's"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointD.CL_D.value = CLmax_func(rho0, S, VD, WS, nD);
CL_D = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointD.CL_D.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointD.CL_D.Attribute.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointD.qD.value = 0.5*rho0*VD^2;
qD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointD.qD.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointD.qD.Attributes.unit = "Pa"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointD.LD.value = CL_D*qD*S;

% POINT S_inv
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS_inverted.VS_inverted.value = VS_inv; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS_inverted.VS_inverted.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS_inverted.nS_inverted.value = nS_inv; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS_inverted.nS_inverted.Attributes.unit = "g's"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS_inverted.CL_S_inverted.value = abs(CLMAX_clean_inverted);
CL_S_inverted = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS_inverted.CL_S_inverted.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS_inverted.CL_S_inverted.Attribute.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS_inverted.qS_inverted.value = 0.5*rho0*VS_inv^2;
qS_inverted = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS_inverted.qS_inverted.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS_inverted.qS_inverted.Attributes.unit = "Pa"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS_inverted.LD.value = CL_S_inverted*qS_inverted*S;

% POINT F
VF = VC;
nF = nmin;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointF.VF.value = VF; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointF.VF.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointF.nF.value = nF; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointF.nF.Attributes.unit = "g's"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointF.CL_F.value = CLmax_func(rho0, S, VF, WS, abs(nF));
CL_F = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointF.CL_F.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointF.CL_F.Attribute.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointF.qF.value = 0.5*rho0*VC^2;
qF = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointF.qF.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointF.qF.Attributes.unit = "Pa"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointF.LF.value = CL_F*qF*S;

% POINT G
nG = nmin;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointG.VG.value = VG; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointG.VG.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointG.nG.value = nG; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointG.nG.Attributes.unit = "g's"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointG.CL_G.value = CLmax_func(rho0, S, VG, WS, abs(nG));
CL_G = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointG.CL_G.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointG.CL_G.Attribute.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointG.qG.value = 0.5*rho0*VG^2;
qG = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointG.qG.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointG.qG.Attributes.unit = "Pa"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointG.LG.value = CL_G*qG*S;

% POINT E
nE = nmin;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointE.VE.value = VE; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointE.VE.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointE.nE.value = nE; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointE.nE.Attributes.unit = "g's"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointE.CL_E.value = CLmax_func(rho0, S, VE, WS, abs(nE));
CL_E = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointE.CL_E.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointE.CL_E.Attribute.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointE.qE.value = 0.5*rho0*VE^2;
qE = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointE.qE.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointE.qE.Attributes.unit = "Pa"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointE.LE.value = CL_E*qE*S;

flight_envelope = figure; 
hold on; grid on; grid minor;
ylim([nmin-1.0 nmax+1.0])
xlim([0 VD+10])
plot(VSpos, npos, ':r', 'LineWidth', 0.25)
plot(VSneg, nneg, ':r', 'LineWidth', 0.25)
plot(V_from1toS, n_from1toS, '-r', 'LineWidth', 1)
plot(V_from1toS_inv, n_from1toS_inv, '-r', 'LineWidth', 1)
plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth', 1)
plot(V_fromAtoC, n_fromAtoC, '-r', 'LineWidth', 1)
plot(V_fromCtoD, n_fromCtoD, '-r', 'LineWidth', 1)
plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth', 1)
plot(V_fromDtoE, n_fromDtoE, '-r', 'LineWidth', 1)
plot(V_fromGtoE, n_fromGtoE, '-r', 'LineWidth', 1)
plot(V_fromStoG, n_fromStoG, '-r', 'LineWidth', 1)
plot(VA, nA, 'k.', 'MarkerSize', 12)
plot(VC, nC, 'k.', 'MarkerSize', 12)
plot(VD, nD, 'k.', 'MarkerSize', 12)
plot(VE, nE, 'k.', 'MarkerSize', 12)
plot(VF, nF, 'k.', 'MarkerSize', 12)
plot(VG, nG, 'k.', 'MarkerSize', 12)
plot(VS, nS, 'k.', 'MarkerSize', 12)
plot(VS_inv, nS_inv, 'k.', 'MarkerSize', 12)
text(VA, nA, 'Point A', 'FontSize', 6)
text(VC, nC, 'Point C', 'FontSize', 6)
text(VD, nD, 'Point D', 'FontSize', 6)
text(VE, nE, 'Point E', 'FontSize', 6)
text(VF, nF, 'Point F', 'FontSize', 6)
text(VG, nG, 'Point G', 'FontSize', 6)
text(VS, nS, 'Point S', 'FontSize', 6)
text(VS_inv, nS_inv, 'Point S invert.', 'FontSize', 6)
xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
title("V~-~n diagram per ", Reg, "Interpreter", "latex")

% EXPORT FIGURE
exportgraphics(flight_envelope, 'Vndiagram.pdf', 'ContentType', 'vector')
exportgraphics(flight_envelope, 'Vndiagram.png', 'ContentType', 'vector')

% STORE FIGURE
Aircraft.Certification.Regulation.SubpartC.Flightloads.V_n_diagram.value = flight_envelope;

% Saving figures inside correct folder
dir = pwd;
fprintf("--------------------------------------");
fprintf('\n');
fprintf('### Saving outpus inside correct Folder ###');
fprintf('\n');
SaveFolder = strcat(dir,'\Output');
fprintf('Saving Vndiagram.pdf in: ');
fprintf('\n');      
fprintf('%s\n', SaveFolder);
% SaveFolder
% Moving file inside correct folder
movefile Vndiagram.pdf Output
movefile Vndiagram.png Output
% -------------------------------------------------------------------------

%% GUST ENVELOPE 
% Vectors with airspeed values
% Two vectors with airspeed values within following ranges:
% 1st ---> [0, VC] where VC = MaxCruiseSpeed
% 2nd ---> [0, VD] where VD = MaxDiveSpeed
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value = linspace(0.0, VC, numb)'; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.Attributes.unit = 'm/s';
V_gust_cruise = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value;

Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_dive.value = linspace(0.0, VD, numb)'; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_dive.Attributes.unit = 'm/s';
V_gust_dive = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_dive.value;
% -------------------------------------------------------------------------

%% x = calcmug(obj, Wingloading, MAC, NormalForceCurveSlope, g)
% This function calculates the MASS RATIO for the selected airplane
% following the CS-VLA airworthiness prescriptions. To have a
% complete documentation check the class file csvla.m
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Mass_ratio.value = calcmug(obj, ...
                                                                                WS, ...      % Wing loading in SI units
                                                                                MAC, ...     % Mean Aerodynamic Chord in meters
                                                                                CLalfa, ...  % Normal force curve slope (practically equal to lift curve slope)
                                                                                rho_operative, ...     % Air density
                                                                                g);          % Gravity acceleration g
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Mass_ratio.Attributes.unit = 'Non dim.';
mu_g = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Mass_ratio.value;

%% x = calckg(obj, MassRatio)
% This function calculates the GUST ALLEVIATION FACTOR for the
% selected airplane and flight conditions following the CS-VLA
% airworthiness prescriprions. To have a complete documentation
% check the class fil csvla.m
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor.value = calckg(obj, mu_g);  
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor.Attributes.unit = 'Non dim.'; 
KG = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor.value;
% -------------------------------------------------------------------------

%% x = calcngust(rho0, NormalForceCurveSlope, GustAlleviationFact, GustSpeedCruiseVect, WingLoading, CruiseSpeed, DiveSpeed, FlagToCalc)
% This function is able to calculates in any possible case a vector
% which contains gust load factors value, following CS-VLA
% airworthiness prescription. To have a complete documentation
% check the class file csvla.m

Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.value = calcngust(obj, rho_operative, ... % Standard atmosphere density
                                                                                            CLalfa, ...               % Normal force curve slope [1/rad]
                                                                                            KG, ...                   % Gust alleviation factor KG
                                                                                            Ude_cruise, ...           % Gust speed at cruise VC
                                                                                            WS, ...                   % Wing loading in SI units
                                                                                            VC, ...                   % Cruise speed from the V - n diagram 
                                                                                            VD, ...                   % Dive speed from the V - n diagram
                                                                                            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.Attributes.case(1))); % A conveniently defined case switch
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.Attributes.unit = 'g';
n_gust_cruise_plus = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.value;

Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise.value = calcngust(obj, rho_operative, ... % Standard atmosphere density
                                                                                            CLalfa, ...               % Normal force curve slope [1/rad]
                                                                                            KG, ...                   % Gust alleviation factor KG
                                                                                            Ude_cruise, ...           % Gust speed at cruise V = VC
                                                                                            WS, ...                   % Wing loading in SI units
                                                                                            VC, ...                   % Cruise speed from the V - n diagram
                                                                                            VD, ...                   % Dive speed from the V - n diagram 
                                                                                            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.Attributes.case(2)));  % A conveniently defined case switch
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise.Attributes.unit = 'g';
n_gust_cruise_neg = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise.value;

Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive.value = calcngust(obj, rho_operative, ...
                                                                                            CLalfa, ...
                                                                                            KG, ...
                                                                                            Ude_dive, ...
                                                                                            WS, ...
                                                                                            VC, ...
                                                                                            VD, ...
                                                                                            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.Attributes.case(1)));
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive.Attributes.unit = 'g';
n_gust_dive_plus = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive.value;

Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive.value = calcngust(obj, rho_operative, ...
                                                                                            CLalfa, ...
                                                                                            KG, ...
                                                                                            Ude_dive, ...
                                                                                            WS, ...
                                                                                            VC, ...
                                                                                            VD, ...
                                                                                            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.Attributes.case(2)));
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive.Attributes.unit = 'g';  
n_gust_dive_neg = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive.value;
% -------------------------------------------------------------------------

%% GUST ENVELOPE DIAGRAM
gust_envelope = figure; 
hold on; grid on; grid minor;
ylim([nmin-2.0 nmax+2.0])
xlim([0 VD+10])
plot(VSpos, npos, ':r', 'LineWidth', 0.25)
plot(VSneg, nneg, ':r', 'LineWidth', 0.25)
plot(V_from1toS, n_from1toS, '-r', 'LineWidth', 1)
plot(V_from1toS_inv, n_from1toS_inv, '-r', 'LineWidth', 1)
plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth', 1)
plot(V_fromAtoC, n_fromAtoC, '-r', 'LineWidth', 1)
plot(V_fromCtoD, n_fromCtoD, '-r', 'LineWidth', 1)
plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth', 1)
plot(V_fromDtoE, n_fromDtoE, '-r', 'LineWidth', 1)
plot(V_fromGtoE, n_fromGtoE, '-r', 'LineWidth', 1)
plot(V_fromStoG, n_fromStoG, '-r', 'LineWidth', 1)
plot(VA, nA, 'k.', 'MarkerSize', 12)
plot(VC, nC, 'k.', 'MarkerSize', 12)
plot(VD, nD, 'k.', 'MarkerSize', 12)
plot(VE, nE, 'k.', 'MarkerSize', 12)
plot(VF, nF, 'k.', 'MarkerSize', 12)
plot(VG, nG, 'k.', 'MarkerSize', 12)
plot(VS, nS, 'k.', 'MarkerSize', 12)
plot(VS_inv, nS_inv, 'k.', 'MarkerSize', 12)
% GUST VALUES
plot(V_gust_cruise, n_gust_cruise_plus, '--k', 'LineWidth', 0.25)
plot(V_gust_cruise, n_gust_cruise_neg, '--k', 'LineWidth', 0.25)
plot(V_gust_dive, n_gust_dive_plus, '--k', 'LineWidth', 0.25)
plot(V_gust_dive, n_gust_dive_neg, '--k', 'LineWidth', 0.25)
text(VA, nA, 'Point A', 'FontSize', 6)
text(VC, nC, 'Point C', 'FontSize', 6)
text(VD, nD, 'Point D', 'FontSize', 6)
text(VE, nE, 'Point E', 'FontSize', 6)
text(VF, nF, 'Point F', 'FontSize', 6)
text(VG, nG, 'Point G', 'FontSize', 6)
text(VS, nS, 'Point S', 'FontSize', 6)
text(VS_inv, nS_inv, 'Point S invert.', 'FontSize', 6)
xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
title("V~-~n diagram per ", Reg, "Interpreter", "latex")

% EXPORT FIGURE
exportgraphics(flight_envelope, 'Gustenvelope.pdf', 'ContentType', 'vector')
exportgraphics(flight_envelope, 'Gustenvelope.png', 'ContentType', 'vector')

% Saving figures inside correct folder
fprintf('Saving Gustenvelope.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile Gustenvelope.pdf Output
movefile Gustenvelope.png Output 
% -----------------------------------------------------------------

%% FINAL ENVELOPE
% Now we must develop a systematic approach to the final envelope, due to
% the fact that the gust envelope and flight envelope must be combined to
% obtain the final envelope. 

final_envelope = figure;
hold on; grid on; grid minor;

% POSITIVE SIDE OF THE FINAL ENVELOPE
syms a b c V
a        = rho0 * CLMAX_clean;
b        = rho_operative * CLalfa * KG * Ude_cruise;
c        = 2 * WS; 
eqn      = a * V^2 - b * V - c ;
Solution = vpasolve(eqn, V);

% WE MUST DECIDE THE NEW MANOEUVRING AIRSPEED
for i = 1:length(Solution) 
    if Solution(i) > 0 
        new_VA = cast(Solution(i), 'double');
        if VA > new_VA
            % Manoeuvring speed
            VA1 = VA;
            nA1 = nmax; 
            % Cruise speed
            VC  = VC;
            nC  = nGust(rho_operative, VC, CLalfa, KG, Ude_cruise, WS);
            
            % FROM 0 TO S
            n_from0toS = linspace(0.0, nS, numb)';
            V_from0toS = VS * ones(numb, 1); 
            
            % FROM S TO A1 
            n_fromStoA1 = linspace(nS, nA1, numb)';
            V_fromStoA1 = Vstall(WS, rho0, CLMAX_clean, n_fromStoA1);
            
            % POINT A1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.point_name.value    = 'Point A1';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.VA1.value            = VA1; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.VA1.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.nA1.value            = nA1; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.nA1.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.value          = CLmax_func(rho0, S, VA1, WS, abs(nA1));
            CL_A1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.Attribute.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.value            = 0.5*rho0*VA1^2;
            qA1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1.value            = CL_A1*qA1*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1.Attributes.unit  = "daN";
            
            if max(n_gust_cruise_plus) > nmax
                % FROM A1 TO C1 
                V_test       = linspace(VA1, VC, numb)';
                n_test       = nGust(rho_operative, V_test, CLalfa, KG, Ude_cruise, WS);
                tol          = 1e-3;
                for i = 1:length(V_test)
                    x = n_test(i);
                    y = x - nmax;
                    if abs(y) < tol
                        row = i; 
                        VC1 = V_test(row);
                    end
                end
                V_fromA1toC1 = linspace(VA1, VC1, numb)';
                n_fromA1toC1 = nmax*ones(length(V_fromA1toC1),1);
                nC1          = nmax;
            
                % FROM C1 TO C 
                V_fromC1toC = linspace(VC1, VC, numb)';
                n_fromC1toC = nGust(rho_operative, V_fromC1toC, CLalfa, KG, Ude_cruise, WS);
                VC          = V_fromC1toC(end);
                nC          = n_fromC1toC(end);

                % FROM C TO C2 
                p = polyfit([n_gust_cruise_plus(end) n_gust_dive_plus(end)], [V_gust_cruise(end) V_gust_dive(end)], 1);
                n_fromCtoC2 = linspace(nC, nmax, numb)';
                V_fromCtoC2 = polyval(p, n_fromCtoC2);
                VC2         = V_fromCtoC2(end);
                nC2         = n_fromCtoC2(end);

                % FROM C2 TO D 
                V_fromC2toD = linspace(VC2, VD, numb)';
                n_fromC2toD = nmax*ones(length(V_fromC2toD),1);
                VD          = V_fromC2toD(end);
                nD          = n_fromC2toD(end);

                % FROM D TO 0
                V_fromDto0 = VD*ones(numb, 1); 
                n_fromDto0 = linspace(nD, 0.0, numb)';
                
                % POINT C
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value    = 'Point C';
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value            = VC; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value            = nC; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value          = CLMAX_clean; % CLmax_func(rho0, S, VC, WS, nC);
                CL_C = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.Attribute.unit = "Non dimensional";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value            = 0.5*rho0*VC^2;
                qC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.Attributes.unit  = "Pa"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.value            = CL_C*qC*S*1e-1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.Attributes.unit  = "daN";
                
                % POINT C2
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.point_name.value    = 'Point C2';
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.VC2.value            = VC2; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.VC2.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.nC2.value            = nC2; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.nC2.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CL_C2.value          = CLmax_func(rho0, S, VC2, WS, abs(nC2));
                CL_C2 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CL_C2.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CL_C2.Attribute.unit = "Non dimensional";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.qC2.value            = 0.5*rho0*VC2^2;
                qC2 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.qC2.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.qC2.Attributes.unit  = "Pa"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.LC2.value            = CL_C2*qC2*S*1e-1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.LC2.Attributes.unit  = "daN";

                % POINT C1
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.point_name.value    = 'Point C1';
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.VC1.value            = VC1; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.VC1.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.nC1.value            = nC1; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.nC1.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CL_C1.value          = CLmax_func(rho0, S, VC1, WS, abs(nC1));
                CL_C1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CL_C1.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CL_C1.Attribute.unit = "Non dimensional";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.qC1.value            = 0.5*rho0*VC1^2;
                qC1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.qC1.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.qC1.Attributes.unit  = "Pa"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.LC1.value            = CL_C1*qC1*S*1e-1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.LC1.Attributes.unit  = "daN";
                
                
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value                 = V_fromStoA1; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.Attribute.unit        = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value           = n_fromStoA1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.Attribute.unit  = "g's";         
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC1.value                        = VC1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC1.Attributes.unit              = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC1.value                  = nmax;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC1.Attributes.unit        = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC2.value                        = VC2;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC2.Attributes.unit              = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC2.value                  = nmax;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC2.Attributes.unit        = "g's";
                
                figure(3);
                hold on;
                plot(V_fromStoA1, n_fromStoA1, 'r', 'LineWidth', 1.5)
                plot(V_fromA1toC1, n_fromA1toC1, 'r', 'LineWidth', 1.5)
                plot(V_fromC1toC, n_fromC1toC, 'r', 'LineWidth', 1.5)
                plot(V_fromCtoC2, n_fromCtoC2, 'r', 'LineWidth', 1.5)
                plot(V_fromC2toD, n_fromC2toD, 'r', 'LineWidth', 1.5)
                plot(V_fromDto0, n_fromDto0, 'r', 'LineWidth', 1.5)
                plot(V_from0toS, n_from0toS, 'r', 'LineWidth', 1.5)

                plot(VD, nD, 'k.', 'MarkerSize', 12);
                plot(VC2, nC2, 'k.', 'MarkerSize', 12)
                plot(VC1, nC1, 'k.', 'MarkerSize', 12)
                plot(VC, nC, 'k.', 'MarkerSize', 12)
                plot(VA1, nA1, 'k.', 'MarkerSize', 12)
                plot(VS, nS, 'k.', 'MarkerSize', 12)

                text(VD, nD, '  D', 'FontSize', 6)
                text(VC2, nC2, '  C2', 'FontSize', 6)
                text(VC, nC, '  C', 'FontSize', 6)
                text(VC1, nC1, '  C1', 'FontSize', 6)
                text(VA1, nA1, '  A', 'FontSize', 6)
                text(VS, nS, '  S', 'FontSize', 6)
                
            elseif abs(n_gust_cruise_plus) < nmax
                % FROM A1 TO C 
                V_fromA1toC = linspace(VA1, VC, numb)';
                n_fromA1toC = nmax*ones(numb, 1);
                VC          = V_fromA1toC(end);
                nC          = nmax;
                
                % FROM C TO D 
                V_fromCtoD = linspace(VC, VD, numb)';
                n_fromCtoD = nmax*ones(numb, 1);
                
                % FROM D TO 0
                V_fromDto0 = VD * ones(numb, 1);
                n_fromDto0 = linspace(nD, 0.0, numb)';
                
                figure(3);
                hold on; 
                plot(V_fromStoA1, n_fromStoA1, 'r', 'LineWidth', 1.5)
                plot(V_fromA1toC, n_fromA1toC, 'r', 'LineWidth', 1.5)
                plot(V_fromCtoD, n_fromCtoD, 'r', 'LineWidth', 1.5)
                plot(V_fromDto0, n_fromDto0, 'r', 'LineWidth', 1.5)
                plot(V_from0toS, n_from0toS, 'r', 'LineWidth', 1.5)

                plot(VD, nD, 'k.', 'MarkerSize', 12);
                plot(VC, nC, 'k.', 'MarkerSize', 12)
                plot(VA1, nA1, 'k.', 'MarkerSize', 12)
                plot(VA, nA, 'k.', 'MarkerSize', 12)
                plot(VS, nS, 'k.', 'MarkerSize', 12)

                text(VD, nD, '  D', 'FontSize', 6)
                text(VC, nC, '  C', 'FontSize', 6)
                text(VA1, nA1, '  A', 'FontSize', 6)
                text(VA, nA, '  A', 'FontSize', 6) 
                text(VS, nS, '  S', 'FontSize', 6)
                
                % POINT C
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value    = 'Point C';
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value            = VC; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value            = nC; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value          = CLMAX_clean; % CLmax_func(rho0, S, VC, WS, nC);
                CL_C = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.Attribute.unit = "Non dimensional";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value            = 0.5*rho0*VC^2;
                qC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.Attributes.unit  = "Pa"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.value            = CL_C*qC*S*1e-1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.Attributes.unit  = "daN";
                
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value                 = V_fromStoA1; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.Attribute.unit        = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value           = n_fromStoA1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.Attribute.unit  = "g's";
                
            end
 
        % -----------------------------------------------------------------------------------------------------------------------------

        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.value                         = VC;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.Attributes.unit               = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.value                   = nC; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.Attributes.unit         = "g's";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VD.value                         = VD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VD.Attributes.unit               = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value                   = nD;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.Attributes.unit         = "g's";
        % ----------------------------------------------------------------------------------------------------------------------------
            
        elseif VA < new_VA
            VA1 = new_VA;
            nA1 = nGust(rho_operative, VA1, CLalfa, KG, Ude_cruise, WS);
            disp(" ")
            % Input to the flight envelope
            Data1 = [  VA1, ...            % Resulting Design Manoeuvring airspeed VA
                       nA1];               % Load factor at V = VA
            disp(" ++++++++++ FINAL ENVELOPE - VA AND nA ++++++++++ ")
            format = ' %6.6f          %6.6f\n';
            label  = ' VA1                nA1\n';
            fprintf(label);
            fprintf(format, Data1.');
            disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

            % POSITIVE SIDE FIRST SEGMENT
            n_from0toS = linspace(0.0, nS, numb)';
            V_from0toS = VS * ones(numb, 1); 

            % POSITIVE STALL SPEED CALCULATION
            % n_fromStoA = linspace(nS, nA, numb)';
            % V_fromStoA = Vstall(WS, rho0, CLMAX_clean, n_fromStoA);

            % POSITIVE STALL SPEED CALCULATION
            n_fromStoA1 = linspace(nS, nA1, numb)';
            V_fromStoA1 = Vstall(WS, rho0, CLMAX_clean, n_fromStoA1);

            % FROM A1 TO C
            V_fromA1toC = linspace(VA1, VC, numb)';
            n_fromA1toC = nGust(rho_operative, V_fromA1toC, CLalfa, KG, Ude_cruise, WS);
            nC         = n_fromA1toC(end);

            % FROM C TO FINAL GUST 
            p = polyfit([n_gust_cruise_plus(end) n_gust_dive_plus(end)], [V_gust_cruise(end) V_gust_dive(end)], 1);
            n_fromCtoFG = linspace(n_gust_cruise_plus(end), nmax, numb)';
            V_fromCtoFG = polyval(p, n_fromCtoFG);

            % FROM FINAL GUST TO D
            n_fromFGtoD = nmax * ones(numb, 1); 
            V_fromFGtoD = linspace(V_fromCtoFG(end), VD, numb)';
            nFG         = nmax;
            VFG         = V_fromFGtoD(1);

            % FROM D TO ZERO 
            V_fromDto0 = VD * ones(numb, 1); 
            n_fromDto0 = linspace(nmax, 0.0, numb)';
            
            % -----------------------------------------------------------------------------------------------------------------------------
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value                 = [V_fromStoA1; V_fromA1toC]; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.Attribute.unit        = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value           = [n_fromStoA1; n_fromA1toC];
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.Attribute.unit  = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_speed.value                     = VFG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_speed.Attributes.unit           = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_load_factor.value               = nFG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_load_factor.Attributes.unit     = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.value                         = VC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.Attributes.unit               = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.value                   = nC; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.Attributes.unit         = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VD.value                         = VD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VD.Attributes.unit               = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value                   = nD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.Attributes.unit         = "g's";

            % POINT A1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.point_name.value    = 'Point A1';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.VA1.value            = VA1; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.VA1.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.nA1.value            = nA1; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.nA1.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.value          = CLmax_func(rho0, S, VA1, WS, abs(nA1));
            CL_A1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.Attribute.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.value            = 0.5*rho0*VA1^2;
            qA1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1.value            = CL_A1*qA1*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1.Attributes.unit  = "daN";

            % POINT C
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value    = 'Point C';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value            = VC; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value            = nC; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value          = CLMAX_clean; % CLmax_func(rho0, S, VC, WS, nC);
            CL_C = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.Attribute.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value            = 0.5*rho0*VC^2;
            qC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.value            = CL_C*qC*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.Attributes.unit  = "daN";

            % POINT FG
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointFG.point_name.value    = 'Point FG';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointFG.VFG.value            = VFG; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointFG.VFG.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointFG.nFG.value            = nFG; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointFG.nFG.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointFG.CL_FG.value          = CLmax_func(rho0, S, VFG, WS, nFG); % CLmax_func(rho0, S, VC, WS, nC);
            CL_FG = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointFG.CL_FG.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointFG.CL_FG.Attribute.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointFG.qFG.value            = 0.5*rho0*VFG^2;
            qFG = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointFG.qFG.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointFG.qFG.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointFG.LFG.value            = CL_FG*qFG*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointFG.LFG.Attributes.unit  = "daN";
            % -----------------------------------------------------------------------------------------------------------------------------
            
            figure(3);
            hold on; 
            plot(V_fromStoA1, n_fromStoA1, 'r', 'LineWidth', 1.5)
            plot(V_fromA1toC, n_fromA1toC, 'r', 'LineWidth', 1.5)
            plot(V_fromCtoFG, n_fromCtoFG, 'r', 'LineWidth', 1.5)
            plot(V_fromFGtoD, n_fromFGtoD, 'r', 'LineWidth', 1.5)
            plot(V_fromDto0, n_fromDto0, 'r', 'LineWidth', 1.5)
            plot(V_from0toS, n_from0toS, 'r', 'LineWidth', 1.5)

            plot(VD, nD, 'k.', 'MarkerSize', 12);
            plot(VFG, nFG, 'k.', 'MarkerSize', 12)
            plot(VC, nC, 'k.', 'MarkerSize', 12)
            plot(VA, nA, 'k.', 'MarkerSize', 12)
            plot(VA1, nA1, 'k.', 'MarkerSize', 12)
            plot(VS, nS, 'k.', 'MarkerSize', 12)
            
            text(VD, nD, '  D', 'FontSize', 6)
            text(VFG, nFG, '   C2', 'FontSize', 6)
            text(VC, nC, '  C', 'FontSize', 6)
            text(VA1, nA1, '  A1', 'FontSize', 6)
            text(VA, nA,   '  A', 'FontSize', 6)
            text(VS, nS, '  S', 'FontSize', 6)
        end
    end
end

% NEGATIVE SIDE OF THE FINAL ENVELOPE
syms a b c V
a        = rho0 * CLMAX_clean_inverted;
b        = rho_operative * CLalfa * KG * Ude_cruise;
c        = 2 * WS; 

if Mass < 150
    a = a;
elseif Mass > 150 
    a = abs(a);
end
eqn      = a * V^2 + b * V - c;  
Solution = vpasolve(eqn, V);
disp(" ")
% Input to the flight envelope
Data1 = [  cast(Solution(1), 'double'), ...            % Resulting Design Manoeuvring airspeed VA
           cast(Solution(2), 'double')];               % Load factor at V = VA
disp(" ++++++++++ FINAL ENVELOPE - VA AND nA ++++++++++ ")
format = ' %6.6f          %6.6f\n';
label  = ' Sol1                Sol2\n';
fprintf(label);
fprintf(format, Data1.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

% WE MUST DECIDE THE NEW MANOEUVRING AIRSPEED
for i = 1:length(Solution) 
    if (Solution(i) > 0) && (Solution(i) > VG)
        new_VG = cast(abs(Solution(i)), 'double');
        if VG > new_VG
            % Manoeuvring airspeed
            VG1 = VG;
            nG1 = nmin; 
            disp(" ")
            % Input to the flight envelope
            Data1 = [  VG1, ...            % Resulting Design Manoeuvring airspeed VA
                       nG1];               % Load factor at V = VA
            disp(" ++++++++++ FINAL ENVELOPE - VA AND nA ++++++++++ ")
            format = ' %6.6f          %6.6f\n';
            label  = ' VG1                nG1\n';
            fprintf(label);
            fprintf(format, Data1.');
            disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
            % Cruise airspeed
            VF  = VF; 
            nF  = nGust_inverted(rho_operative, VF, CLalfa, KG, Ude_cruise, WS);
            
            %FROM 0 TO S_inverted 
            nS_inv        = -1.0;
            n_from0toSinv = linspace(0.0, nS_inv, numb)';
            VS_inv        = Vstall(WS, rho0, abs(CLMAX_clean_inverted), abs(nS_inv));
            V_from0toSinv = VS_inv*ones(numb, 1);
            
            % FROM S_inverted TO G
            n_fromSinvtoG = linspace(nS_inv, nG, numb)';
            V_fromSinvtoG = Vstall(WS, rho0, abs(CLMAX_clean_inverted), abs(n_fromSinvtoG));
            
            if abs(min(n_gust_cruise_neg)) > abs(nmin)
                % FROM F TO F1
                V_test      = linspace(VG, VF, numb)';
                n_test      = nGust_inverted(rho_operative, V_test, CLalfa, KG, Ude_cruise, WS);
                tol         = 1e-3;
                for i = 1:length(V_test)
                    X           = abs(n_test(i) - nmin);
                    if X < tol 
                        row = i; 
                        VF1 = V_test(row);
                    end
                end
                V_fromGtoF1 = linspace(VG, VF1, numb)';
                n_fromGtoF1 = nmin*ones(length(V_fromGtoF1),1);
                nF1         = nmin;
                
                % FROM F1 TO F 
                V_fromF1toF = linspace(VF1, VF, numb)';
                n_fromF1toF = nGust_inverted(rho_operative, V_fromF1toF, CLalfa, KG, Ude_cruise, WS);
                nF          = n_fromF1toF(end);
                
                % FROM F TO F2 
                n_fromFtoF2 = linspace(nF, nmin, numb)';
                p           = polyfit([n_gust_cruise_neg(end) n_gust_dive_neg(end)], ...
                                 [V_gust_cruise(end) V_gust_dive(end)], 1);
                V_fromFtoF2 = polyval(p, n_fromFtoF2);   
                nF2         = n_fromFtoF2(end);
                VF2         = V_fromFtoF2(end);
                
                % FROM F2 TO E
                V_fromF2toE = linspace(VF2, VE, numb)';
                n_fromF2toE = nmin*ones(length(V_fromF2toE),1);
                
                % FROM E TO 0
                V_fromEto0 = VE * ones(numb, 1);
                n_fromEto0 = linspace(nE, 0.0, numb)';
                
                figure(3);
                hold on; 
                plot(V_from0toSinv, n_from0toSinv, 'r', 'LineWidth', 1.5)
                plot(V_fromSinvtoG, n_fromSinvtoG, 'r', 'LineWidth', 1.5)
                plot(V_fromGtoF1, n_fromGtoF1, 'r', 'LineWidth', 1.5)
                plot(V_fromF1toF, n_fromF1toF, 'r', 'LineWidth', 1.5)
                plot(V_fromFtoF2, n_fromFtoF2, 'r', 'LineWidth', 1.5)
                plot(V_fromF2toE, n_fromF2toE, 'r', 'LineWidth', 1.5)
                plot(V_fromEto0, n_fromEto0, 'r', 'LineWidth', 1.5)
                
                plot(VS_inv, nS_inv, 'k.', 'MarkerSize', 12);
                plot(VG, nG, 'k.', 'MarkerSize', 12)
                plot(VE, nE, 'k.', 'MarkerSize', 12)
                plot(VF, nF, 'k.', 'MarkerSize', 12)
                plot(VF1, nF1, 'k.', 'MarkerSize', 12)
                plot(VF2, nF2, 'k.', 'MarkerSize', 12)
                
                text(VF, nF, '  F', 'FontSize', 6)
                text(VF1, nF1, '  F1', 'FontSize', 6)
                text(VF2, nF2, '  F2', 'FontSize', 6)
                text(VG, nG, '  G', 'FontSize', 6)
                text(VE, nE, '  E', 'FontSize', 6)
                text(VS_inv, nS_inv, '  S inv.', 'FontSize', 6)
                
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VE.value                 = VE;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VE.Attributes.unit       = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value           = nE;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.Attributes.unit = "g's";
            
            % -----------------------------------------------------------------------------------------------------------------            
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value                 = V_fromSinvtoG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.Attributes.unit       = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.value           = n_fromSinvtoG;  
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.Attributes.unit = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF1.value                        = VF1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF1.Attributes.unit              = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF1.value                  = nF1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF1.Attributes.unit        = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF2.value                        = VF2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF2.Attributes.unit              = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF2.value                  = nF2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF2.Attributes.unit        = "g's";
            % ------------------------------------------------------------------------------------------------------------------
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.value                 = VF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.Attributes.unit       = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.value           = nF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.Attributes.unit = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VE.value                 = VE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VE.Attributes.unit       = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value           = nE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.Attributes.unit = "g's";

            % POINT F1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF1.point_name.value    = 'Point F1';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF1.VF1.value            = VF1; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF1.VF1.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF1.nF1.value            = nF1; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF1.nF1.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF1.CL_F1.value          = CLmax_func(rho0, S, VF1, WS, abs(nF1));
            CL_F1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF1.CL_F1.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF1.CL_F1.Attribute.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF1.qF1.value            = 0.5*rho0*VF1^2;
            qF1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF1.qF1.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF1.qF1.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF1.LF1.value            = CL_F1*qF1*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF1.LF1.Attributes.unit  = "daN";

            % POINT F2
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF2.point_name.value    = 'Point F2';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF2.VF2.value            = VF2; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF2.VF2.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF2.nF2.value            = nF2; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF2.nF2.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF2.CL_F2.value          = CLmax_func(rho0, S, VF2, WS, abs(nF2));
            CL_F2 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF2.CL_F2.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF2.CL_F2.Attribute.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF2.qF2.value            = 0.5*rho0*VF2^2;
            qF2 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF2.qF2.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF2.qF2.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF2.LF2.value            = CL_F2*qF2*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF2.LF2.Attributes.unit  = "daN";

            % POINT F
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value    = 'Point F';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.value            = VF; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.value            = nF; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value          = CLmax_func(rho0, S, VF, WS, abs(nF));
            CL_F = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.Attribute.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value            = 0.5*rho0*VF^2;
            qF = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.value            = CL_F*qF*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.Attributes.unit  = "daN";

% ------------------------------------------------------------------------------------------------------------------
            elseif abs(min(n_gust_cruise_neg)) < abs(nmin)
                
                % FROM G TO F
                V_fromGtoF = linspace(VG, VF, numb)';
                n_fromGtoF = nmin * ones(numb, 1); 
                VF         = V_fromGtoF(end);
                nF         = n_fromGtoF(end);
                
                % FROM F TO E 
                V_fromFtoE = linspace(VF, VD, numb)';
                n_fromFtoE = nmin * ones(numb, 1); 
                
                % FROM E TO 0 
                V_fromEto0 = VE * ones(numb, 1);
                n_fromEto0 = linspace(nmin, 0.0, numb)';
                
                % POINT F
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value    = 'Point F';
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.value            = VF; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.value            = nF; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value          = CLmax_func(rho0, S, VF, WS, abs(nF));
                CL_F = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.Attribute.unit = "Non dimensional";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value            = 0.5*rho0*VF^2;
                qF = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.Attributes.unit  = "Pa"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.value            = CL_F*qF*S*1e-1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.Attributes.unit  = "daN";
                
                figure(3);
                hold on; 
                plot(V_from0toSinv, n_from0toSinv, 'r', 'LineWidth', 1.5)
                plot(V_fromSinvtoG, n_fromSinvtoG, 'r', 'LineWidth', 1.5)
                plot(V_fromGtoF, n_fromGtoF, 'r', 'LineWidth', 1.5)
                plot(V_fromFtoE, n_fromFtoE, 'r', 'LineWidth', 1.5)
                plot(V_fromEto0, n_fromEto0, 'r', 'LineWidth', 1.5)

                plot(VS_inv, nS_inv, 'k.', 'MarkerSize', 12);
                plot(VG, nG, 'k.', 'MarkerSize', 12)
                plot(VE, nE, 'k.', 'MarkerSize', 12)
                plot(VF, nF, 'k.', 'MarkerSize', 12)

                text(VF, nF, '  F', 'FontSize', 6)
                text(VG, nG, '  G', 'FontSize', 6)
                text(VE, nE, '  E', 'FontSize', 6)
                text(VS_inv, nS_inv, '  S inv.', 'FontSize', 6)           
                
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VE.value                 = VE;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VE.Attributes.unit       = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value           = nE;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.Attributes.unit = "g's";
            end
            
        elseif VG < new_VG
            VG1 = new_VG;
            nG1 = nGust_inverted(rho_operative, VG1, CLalfa, KG, Ude_cruise, WS);
            disp(" ")
            % Input to the flight envelope
            Data1 = [  VG1, ...            % Resulting Design Manoeuvring airspeed VA
                       nG1];               % Load factor at V = VA
            disp(" ++++++++++ FINAL ENVELOPE - VG AND nG ++++++++++ ")
            format = ' %6.6f          %6.6f\n';
            label  = ' VG1                nG1\n';
            fprintf(label);
            fprintf(format, Data1.');
            disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

            % FROM 0 TO S_inverted
            nS_inv        = -1.0;
            n_from0toSinv = linspace(0.0, nS_inv, numb)';
            VS_inv        = Vstall(WS, rho0, abs(CLMAX_clean_inverted), abs(nS_inv));
            V_from0toSinv = VS_inv*ones(numb, 1);

            % FROM S_inverted TO G1
            n_fromSinvtoG1 = linspace(nS_inv, nG1, numb)';
            V_fromSinvtoG1 = Vstall(WS, rho0, abs(CLMAX_clean_inverted), abs(n_fromSinvtoG1));

            % FROM S_inverted TO G
            % n_fromGinvtoG1 = linspace(nG, nG1, numb)';
            % V_fromGinvtoG1 = Vstall(WS, rho0, abs(CLMAX_clean_inverted), abs(n_fromGinvtoG1));

            % FROM G TO F
            V_fromG1toF = linspace(VG1, VC, numb)';
            n_fromG1toF = nGust_inverted(rho_operative, V_fromG1toF, CLalfa, KG, Ude_cruise, WS);
            VF = VC;
            nF = n_fromG1toF(end);

            % FROM F TO E
            n_fromFtoE = linspace(n_fromG1toF(end), n_gust_dive_neg(end), numb)';
            nE         = n_fromFtoE(end);
            p          = polyfit([n_gust_cruise_neg(end) n_gust_dive_neg(end)], ...
                                 [V_gust_cruise(end) V_gust_dive(end)], 1);
            V_fromFtoE = polyval(p, double(n_fromFtoE));
            VE         = V_fromFtoE(end);

            % FROM E TO ZERO 
            n_fromEto0 = linspace(n_fromFtoE(end), 0.0, numb)';
            V_fromEto0 = VD*ones(numb, 1);
            
            % -----------------------------------------------------------------------------------------------------------------            
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value                 = [V_fromSinvtoG1; V_fromG1toF];
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.Attributes.unit       = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.value           = [n_fromSinvtoG1; n_fromG1toF]; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.Attributes.unit = "g's";
            % ------------------------------------------------------------------------------------------------------------------
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.value                 = VF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.Attributes.unit       = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.value           = nF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.Attributes.unit = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VE.value                 = VE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VE.Attributes.unit       = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value           = nE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.Attributes.unit = "g's";
            % ------------------------------------------------------------------------------------------------------------------

            % POINT F
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value    = 'Point F';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.value            = VF; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.value            = nF; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value          = CLmax_func(rho0, S, VF, WS, abs(nF));
            CL_F = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.Attribute.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value            = 0.5*rho0*VF^2;
            qF = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.value            = CL_F*qF*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.Attributes.unit  = "daN";

            % POINT G1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.point_name.value    = 'Point G1';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.VG1.value            = VG1; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.VG1.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.nG1.value            = nG1; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.nG1.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CL_G1.value          = CLmax_func(rho0, S, VG1, WS, abs(nG1));
            CL_G1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CL_G1.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CL_G1.Attribute.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.qG1.value            = 0.5*rho0*VG1^2;
            qG1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.qG1.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.qG1.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1.value            = CL_G1*qG1*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1.Attributes.unit  = "daN";
            
            figure(3);
            hold on;
            plot(V_from0toSinv, n_from0toSinv, 'r', 'LineWidth', 1.5)
            plot(V_fromSinvtoG1, n_fromSinvtoG1, 'r', 'LineWidth', 1.5)
            plot(V_fromG1toF, n_fromG1toF, 'r', 'LineWidth', 1.5)
            plot(V_fromFtoE, n_fromFtoE, 'r', 'LineWidth', 1.5)
            plot(V_fromEto0, n_fromEto0, 'r', 'LineWidth', 1.5)
            
            plot(VS_inv, nS_inv, 'k.', 'MarkerSize', 12);
            plot(VG, nG, 'k.', 'MarkerSize', 12)
            plot(VG1, nG1, 'k.', 'MarkerSize', 12)
            plot(VE, nE, 'k.', 'MarkerSize', 12)
            plot(VF, nF, 'k.', 'MarkerSize', 12)
            
            text(VF, nF, '  F', 'FontSize', 6)
            text(VG, nG, '  G', 'FontSize', 6)
            text(VG1, nG1, '  G1', 'FontSize', 6)
            text(VE, nE, '  E', 'FontSize', 6)
            text(VS_inv, nS_inv, '  S inv.', 'FontSize', 6)
            
        end
    end
end



% AIRCRAFT STRUCT VARIABLE FILLING 
% ------------------------------------------------------------------------------------------------------------------
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value                 = [V_fromStoA1; V_fromA1toC];
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.Attribute.unit        = "m/s";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value           = [n_fromStoA1; n_fromA1toC];
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.Attribute.unit  = "g's";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.value                         = VC;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.Attributes.unit               = "m/s";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_speed.value                     = V_fromCtoFG(end);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_speed.Attributes.unit           = "m/s";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_load_factor.value               = nmax;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_load_factor.Attributes.unit     = "g's"; 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value                 = [V_fromSinvtoG1; V_fromG1toF];
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.Attributes.unit       = "m/s";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.value           = [n_fromSinvtoG1; n_fromG1toF]; 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.Attributes.unit = "g's";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.value                   = nC; 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.Attributes.unit         = "g's";
% ------------------------------------------------------------------------------------------------------------------
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value                               = VD; 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.Attributes.unit                     = "m/s"; 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value           = nmax;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.Attributes.unit = "g's";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.value                 = VF;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.Attributes.unit       = "m/s";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.value           = nF;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.Attributes.unit = "g's";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VE.value                 = VE;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VE.Attributes.unit       = "m/s";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value           = nE;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.Attributes.unit = "g's";
% ------------------------------------------------------------------------------------------------------------------

%% FINAL ENVELOPE DIAGRAM 

% final_envelope = figure; 
% hold on; grid on; grid minor;
% ylim([nmin-2.0 nmax+2.0])
% xlim([0 VD+10])
% plot(V_from0toSinv, n_from0toSinv, 'r', 'LineWidth', 1.5)
% plot(V_fromStoA1, n_fromStoA1, 'r', 'LineWidth', 1.5)
% plot(V_fromA1toC, n_fromA1toC, 'r', 'LineWidth', 1.5)
% plot(V_fromCtoFG, n_fromCtoFG, 'r', 'LineWidth', 1.5)
% plot(V_fromFGtoD, n_fromFGtoD, 'r', 'LineWidth', 1.5)
% plot(V_fromDto0, n_fromDto0, 'r', 'LineWidth', 1.5)
% plot(V_from0toS, n_from0toS, 'r', 'LineWidth', 1.5)
% plot(V_fromSinvtoG1, n_fromSinvtoG1, 'r', 'LineWidth', 1.5)
% plot(V_fromG1toF, n_fromG1toF, 'r', 'LineWidth', 1.5)
% plot(V_fromFtoE, n_fromFtoE, 'r', 'LineWidth', 1.5)
% plot(V_fromEto0, n_fromEto0, 'r', 'LineWidth', 1.5)
% plot(VD, nD, 'k.', 'MarkerSize', 12);
% plot(V_fromCtoFG(end), n_fromCtoFG(end), 'k.', 'MarkerSize', 12)
% plot(VF, nF, 'k.', 'MarkerSize', 12);
% plot(VG, nG, 'k.', 'MarkerSize', 12);
% plot(VG1, nG1, 'k.', 'MarkerSize', 12);
% plot(VC, nC, 'k.', 'MarkerSize', 12)
% plot(VA, nA, 'k.', 'MarkerSize', 12)
% plot(VA1, nA1, 'k.', 'MarkerSize', 12)
% plot(VS, nS, 'k.', 'MarkerSize', 12)
% plot(VE, nE, 'k.', 'MarkerSize', 12)
% plot(VS_inv, nS_inv, 'k.', 'MarkerSize', 12)
% plot(V_from0toSinv(1), nS_inv, 'k.', 'MarkerSize', 12)
% text(VD, nD, '  D', 'FontSize', 6)
% text(VC, nC, '  C', 'FontSize', 6)
% text(VA, nA, '  A', 'FontSize', 6)
% text(VA1, nA1, '  A1', 'FontSize', 6)
% text(VS, nS, '  S', 'FontSize', 6)
% text(VF, nF, '  F', 'FontSize', 6)
% text(VG, nG, '  G', 'FontSize', 6)
% text(VG1, nG1, '  G1', 'FontSize', 6)
% text(VE, nE, '  E', 'FontSize', 6)
% text(VS_inv, nS_inv, '  S inv.', 'FontSize', 6)
% text(V_from0toSinv(1), nS_inv, '  S inv.', 'FontSize', 6)

ylim([nF-1.0 nC+1.0])
xlim([0 VD+10])
xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
title("V~-~n diagram per ", Reg, "Interpreter", "latex")

% EXPORT FIGURE
exportgraphics(final_envelope, 'Finalenvelope.pdf', 'ContentType', 'vector')
exportgraphics(final_envelope, 'Finalenvelope.png', 'ContentType', 'vector')

% Saving figures inside correct folder
fprintf('Saving Finalenvelope.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile Finalenvelope.pdf Output
movefile Finalenvelope.png Output 
% -----------------------------------------------------------------

% AIRCRAFT STRUCT VARIABLE FILLING

% POINT S
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value    = 'Point S';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.value            = VS; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.Attributes.unit  = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.nS.value            = 1.0; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.nS.Attributes.unit  = "g's"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value          = CLmax_func(rho0, S, VS, WS, 1.0);
CL_S = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.Attribute.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.value            = 0.5*rho0*VS^2;
qS = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.Attributes.unit  = "Pa"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS.value            = CL_S*qS*S*1e-1;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS.Attributes.unit  = "daN";

% POINT A
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.point_name.value    = 'Point A';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value            = VA; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.Attributes.unit  = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value            = nA; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.Attributes.unit  = "g's"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value          = CLmax_func(rho0, S, VA, WS, nA);
CL_A = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.Attribute.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value            = 0.5*rho0*VA^2;
qA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.Attributes.unit  = "Pa"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA.value            = CL_A*qA*S*1e-1;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA.Attributes.unit  = "daN";

% POINT D
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value = 'Point D';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value            = VD; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.Attributes.unit  = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value            = nmax; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.Attributes.unit  = "g's"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value          = CLmax_func(rho0, S, VD, WS, nmax);
CL_D = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.Attribute.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value            = 0.5*rho0*VD^2;
qD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.Attributes.unit  = "Pa"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD.value            = CL_D*qD*S*1e-1;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD.Attributes.unit  = "daN";

% POINT S_inv
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inverted.point_name.value             = 'Point S inv.';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inverted.VS_inverted.value            = VS_inv; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inverted.VS_inverted.Attributes.unit  = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inverted.nS_inverted.value            = nS_inv; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inverted.nS_inverted.Attributes.unit  = "g's"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inverted.CL_S_inverted.value          = abs(CLMAX_clean_inverted);
CL_S_inverted = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inverted.CL_S_inverted.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inverted.CL_S_inverted.Attribute.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inverted.qS_inverted.value            = 0.5*rho0*VS_inv^2;
qS_inverted = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inverted.qS_inverted.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inverted.qS_inverted.Attributes.unit  = "Pa"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inverted.LS_inverted.value            = CL_S_inverted*qS_inverted*S*1e-1;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inverted.LS_inverted.Attributes.unit  = "daN";

% POINT G
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.point_name.value    = 'Point G';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG.value            = VG; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG.Attributes.unit  = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.nG.value            = nG; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.nG.Attributes.unit  = "g's"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value          = CLmax_func(rho0, S, VG, WS, abs(nG));
CL_G = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.Attribute.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value            = 0.5*rho0*VG^2;
qG = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.Attributes.unit  = "Pa"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG.value            = CL_G*qG*S*1e-1;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG.Attributes.unit  = "daN";

% POINT E
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.point_name.value    = 'Point E';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.VE.value            = VE; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.VE.Attributes.unit  = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.nE.value            = nE; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.nE.Attributes.unit  = "g's"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.value          = CLmax_func(rho0, S, VE, WS, abs(nE));
CL_E = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.Attribute.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value            = 0.5*rho0*VD^2;
qE = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.Attributes.unit  = "Pa"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE.value            = CL_E*qE*S*1e-1;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE.Attributes.unit  = "daN";

