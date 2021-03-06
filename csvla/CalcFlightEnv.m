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

%% MEAN GEOMETRIC CHORD
b   = Aircraft.Geometry.Wing.b.value; 
S   = Aircraft.Geometry.Wing.S.value; 
MGC = S / b;

% STORE INSIDE THE STRUCT VARIABLE
Aircraft.Geometry.Wing.mgc.value           = MGC;
Aircraft.Geometry.Wing.mgc.Attributes.unit = "m";

%% MEAN AERODYNAMIC CHORD 
c_root      = Aircraft.Geometry.Wing.croot.value;
c_tip       = Aircraft.Geometry.Wing.ctip.value;
taper_ratio = c_tip / c_root;
MAC         = mean_aerodynamic_chord(c_root, taper_ratio);

% STORE INSIDE THE STRUCT VARIABLE
% Aircraft.Geometry.Wing.taper_ratio.value = taper_ratio;
% Aircraft.Geometry.Wing.taper_ratio.Attributes.unit = "Non dimensional";
Aircraft.Geometry.Wing.mac.value           = MAC;
Aircraft.Geometry.Wing.mac.Attributes.unit = "m";

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
Aircraft.Certification.ISA_Condition.Sea_Level.rho0.Attributes.unit = "kg/m^3";
Aircraft.Certification.ISA_Condition.Sea_Level.Altitude.value = h0;
Aircraft.Certification.ISA_Condition.Sea_Level.Altitude.Attributes.unit = "m";
Aircraft.Certification.ISA_Condition.Sea_Level.T0.value = T0; 
Aircraft.Certification.ISA_Condition.Sea_Level.T0.Attributes.unit = 'Kelvin'; 
Aircraft.Certification.ISA_Condition.Sea_Level.p0.value = p0; 
Aircraft.Certification.ISA_Condition.Sea_Level.p0.Attributes.unit = 'Pa';
Aircraft.Certification.ISA_Condition.Sea_Level.Speed_of_sound0.value = a0; 
Aircraft.Certification.ISA_Condition.Sea_Level.Speed_of_sound0.Attributes.unit = "m/s";
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
Aircraft.Certification.ISA_Condition.Operative_ceiling.rho0.Attributes.unit = "kg/m^3";
Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.value = h_operative;
Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.Attributes.unit = "m";
Aircraft.Certification.ISA_Condition.Operative_ceiling.T0.value = T_operative; 
Aircraft.Certification.ISA_Condition.Operative_ceiling.T0.Attributes.unit = 'K'; 
Aircraft.Certification.ISA_Condition.Operative_ceiling.p0.value = p_operative; 
Aircraft.Certification.ISA_Condition.Operative_ceiling.p0.Attributes.unit = 'Pa';
Aircraft.Certification.ISA_Condition.Operative_ceiling.Speed_of_sound0.value = a_operative; 
Aircraft.Certification.ISA_Condition.Operative_ceiling.Speed_of_sound0.Attributes.unit = "m/s";
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
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.rho0.Attributes.unit = "kg/m^3";
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.value = h_theoretical;
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.Attributes.unit = "m";
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.T0.value = T_theoretical; 
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.T0.Attributes.unit = 'K'; 
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.p0.value = p_theoretical; 
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.p0.Attributes.unit = 'Pa';
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Speed_of_sound0.value = a_theoretical; 
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Speed_of_sound0.Attributes.unit = "m/s";

%% NUMBER OF ELEMENTS
numb = 1e3;

% REGULATION APPLIED
Reg = Aircraft.Certification.Regulation.value;

% -------------------------------------------------------------------------
%% CS-VLA 333 Flight envelope
%  (a) GENERAL. Compliance with the strength requirements of this subpart
%      must be shown at any combination of airspeed and load factor on and
%      within the boundaries of a flight envelope (similar to the one in
%      sub-paragraph(d) of this paragraph) that represents the envelope of
%      the flight loading conditions specified by the manoeuvring and gust
%      criteria of sub-paragraphs(b) and (c) of this paragraph
%      respectively.
%
%  (b) MANOEUVRING ENVELOPE. Except where limited by maximum static lift
%      coefficients, the aeroplane is assumed to be subjected to
%      symmetrical manoeuvres resulting in the following limit load
%      factors: 
%      (1) the positive manoeuvring load factor specified in CS-VLA 337 at
%          speeds up to VD; 
%      (2) the negative manoeuvring load factor specified in CS-VLA 337 at
%          VC; 
%      (3) factors varying linearly with speed from the specified value at
%          VC to 0.0 at VD.
%
%  (c) GUST ENVELOPE.
%      (1) The aeroplane is assumed to be subjected to symmetrical vertical
%          gusts in level flight. The resulting limit load factors must
%          correspond to the conditions determined as follows: 
%          (i) positive (up) and negative (down) gusts of 15.24 m/s at VC
%              must be considered;
%         (ii) positive and negative gusts of 7.62 m/s at VD must be
%              considered.
%      (2) The following assumptions must be made: 
%          (i) the shape of the gust is 
%                  Ude   /          2*pi*S   \
%              U = --- * | 1 - cos( ------ ) |
%                   2    \          25*MGC   /
%              where 
%              S   = distance penetrated into gust (m);
%              MGC = mean geometric chord of wing (m);
%              Ude = derived gust velocity referred to in
%                    sub-paragraph(c)(1) (m/s)
%         (ii) gust load factors vary linearly with speed between VC and
%              VD.
%
%  (d) FLIGHT ENVELOPE. Figure inside the CS-VLA airworthiness
%      requirements. NOTE: point G need not be investigated when the
%      supplementary condition specified in CS-VLA 369 is investigated.

% Aircraft.Certification.Regulation.SubpartC.Flightloads.Attributes.cs = " 333 ";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.Attributes.cs = " 333(d) ";
% -------------------------------------------------------------------------
%% WING LOADING DEFINITION
% x = calcn(obj, nmax) - from csvla.m
% This function defines a vector with load factor values between two pre-
% scribed limits. Check the class file csvla.m to have a complete
% documentation.
g    = Aircraft.Constants.g.value;
S1   = 0.0;
S    = Aircraft.Geometry.Wing.S.value + S1;
x    = 0;
Mass = Aircraft.Weight.I_Level.W_maxTakeOff.value + x;
Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value = (Mass * g) / (S);
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
Aircraft.Certification.Regulation.SubpartC.Flightloads.Cruise_Speed_VC.Attributes.cs   = " 335(a)(1) ";
Aircraft.Certification.Regulation.SubpartC.Flightloads.nC.value = nmax; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.nC.Attributes.unit = "g's";
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
Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.Attributes.cs   = " 335(b)(1)/335(b)(2) ";
Aircraft.Certification.Regulation.SubpartC.Flightloads.nD.value = nmax; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.nD.Attributes.unit = "g's";
VD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value;
nD = nmax;

% INVERTED FLIGHT DIVE SPEED 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VE.value = VD; % Speed at points E and D are equal 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VE.Attributes.unit = "m/s";
Aircraft.Certification.Regulation.SubpartC.Flightloads.nE.value = nmin; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.nE.Attributes.unit = "g's";
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
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.Attributes.cs = " 335(c)(1)/335(c)(2) "; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.nA.value = nmax; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.nA.Attributes.unit = "g's";
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
Aircraft.Certification.Regulation.SubpartC.Flightloads.nG.Attributes.unit = "g's";
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

%% Point F definition 
VF = VC;
VE = VD;
V_fromGtoF = linspace(VG, VF, numb)';
n_fromGtoF = nmin*ones(numb, 1);

%% Point E definition 
V_fromFto0 = linspace(VF, VE, numb)';
n_fromFto0 = linspace(nmin, 0.0, numb)';

%% Flight envelope limit 
V_fromDto0 = VD*ones(numb, 1);
n_fromDto0 = linspace(nD, 0.0, numb)';

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

disp(" ")
% Input to the flight envelope
Data1 = [VC];                 % VG = VD on the negative side of V - n diagram
disp(" ++++++++++ DATA USED TO PLOT V - N DIAGRAM ++++++++++ ")
format = ' %6.6f          \n';
label  = ' VC\n';
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
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS.CL_S.Attributes.unit = "Non dimensional";
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
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointA.CL_A.Attributes.unit = "Non dimensional";
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
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointC.CL_C.Attributes.unit = "Non dimensional";
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
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointD.CL_D.Attributes.unit = "Non dimensional";
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
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointS_inverted.CL_S_inverted.Attributes.unit = "Non dimensional";
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
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointF.CL_F.Attributes.unit = "Non dimensional";
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
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointG.CL_G.Attributes.unit = "Non dimensional";
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
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointE.CL_E.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointE.qE.value = 0.5*rho0*VE^2;
qE = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointE.qE.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointE.qE.Attributes.unit = "Pa"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.PointE.LE.value = CL_E*qE*S;

flight_envelope = figure; 
hold on; grid on; grid minor;
ylim([nmin-1.0 nmax+1.0])
% xlim([0 VD+10])
% ylim 'padded';
xlim 'padded';
plot(VSpos, npos, ':r', 'LineWidth', 0.25)
plot(VSneg, nneg, ':r', 'LineWidth', 0.25)
plot(V_from1toS, n_from1toS, '-r', 'LineWidth', 1)
plot(V_from1toS_inv, n_from1toS_inv, '-r', 'LineWidth', 1)
plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth', 1)
plot(V_fromAtoC, n_fromAtoC, '-r', 'LineWidth', 1)
plot(V_fromCtoD, n_fromCtoD, '-r', 'LineWidth', 1)
plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth', 1)
plot(V_fromDto0, n_fromDto0, '-r', 'LineWidth', 1)
plot(V_fromGtoF, n_fromGtoF, '-r', 'LineWidth', 1)
plot(V_fromFto0, n_fromFto0, '-r', 'LineWidth', 1)
plot(V_fromStoG, n_fromStoG, '-r', 'LineWidth', 1)

plot(VA, nA, 'k.', 'MarkerSize', 12)
plot(VC, nC, 'k.', 'MarkerSize', 12)
plot(VD, nD, 'k.', 'MarkerSize', 12)
plot(VE, 0.0, 'k.', 'MarkerSize', 12)
plot(VF, nF, 'k.', 'MarkerSize', 12)
plot(VG, nG, 'k.', 'MarkerSize', 12)
plot(VS, nS, 'k.', 'MarkerSize', 12)
plot(VS_inv, nS_inv, 'k.', 'MarkerSize', 12)
text(VA, nA, 'Point A', 'FontSize', 6)
text(VC, nC, 'Point C', 'FontSize', 6)
text(VD, nD, 'Point D', 'FontSize', 6)
text(VE, 0.0, 'Point E', 'FontSize', 6)
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
% NOTE: It's important to remember that in this version of the code the air
%       density for all the wind gust calculations are referred to the
%       operational altitude of the selected aircraft, whilst the flight
%       envelope is referred to Sea Level, Standard Atmosphere air density.
% Vectors with airspeed values
% Two vectors with airspeed values within following ranges:
% 1st ---> [0, VC] where VC = MaxCruiseSpeed
% 2nd ---> [0, VD] where VD = MaxDiveSpeed
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value = linspace(0.0, VC, numb)'; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.Attributes.unit = "m/s";
V_gust_cruise = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value;

Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_dive.value = linspace(0.0, VD, numb)'; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_dive.Attributes.unit = "m/s";
V_gust_dive = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_dive.value;
% -------------------------------------------------------------------------

%% x = calcmug(obj, Wingloading, MAC, NormalForceCurveSlope, g)
% This function calculates the MASS RATIO for the selected airplane
% following the CS-VLA airworthiness prescriptions. To have a
% complete documentation check the class file csvla.m
% -----------------------------------------------------------------------------------------------------------------------
% SEA LEVEL
% -----------------------------------------------------------------------------------------------------------------------
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Mass_ratio_sea_level.value = calcmug(obj, ...
                                                                                WS, ...      % Wing loading in SI units
                                                                                MGC, ...     % Mean Geometric Chord in meters
                                                                                CLalfa, ...  % Normal force curve slope (practically equal to lift curve slope)
                                                                                rho0, ...     % Air density
                                                                                g);          % Gravity acceleration g
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Mass_ratio_sea_level.Attributes.unit = 'Non dimensional';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Mass_ratio_sea_level.Attributes.cs = " 341 ";
mu_g_sea_level = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Mass_ratio_sea_level.value;
% -----------------------------------------------------------------------------------------------------------------------
% OPERATIVE ALTITUDE
% -----------------------------------------------------------------------------------------------------------------------
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Mass_ratio.value = calcmug(obj, ...
                                                                                WS, ...      % Wing loading in SI units
                                                                                MGC, ...     % Mean Geometric Chord in meters
                                                                                CLalfa, ...  % Normal force curve slope (practically equal to lift curve slope)
                                                                                rho_operative, ...     % Air density
                                                                                g);          % Gravity acceleration g
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Mass_ratio.Attributes.unit = 'Non dimensional';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Mass_ratio.Attributes.cs = " 341 ";
mu_g = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Mass_ratio.value;

%% DENSITY: Sea Level - x = calckg(obj, MassRatio)
% This function calculates the GUST ALLEVIATION FACTOR for the
% selected airplane and flight conditions following the CS-VLA
% airworthiness prescriprions. To have a complete documentation
% check the class fil csvla.m
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor_sea_level.value = calckg(obj, mu_g_sea_level);  
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor_sea_level.Attributes.unit = 'Non dimensional'; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor_sea_level.Attributes.cs = " 341 ";
KG_sea_level = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor_sea_level.value;
% -------------------------------------------------------------------------

%% DENSITY: Operative altitude - x = calckg(obj, MassRatio)
% This function calculates the GUST ALLEVIATION FACTOR for the
% selected airplane and flight conditions following the CS-VLA
% airworthiness prescriprions. To have a complete documentation
% check the class fil csvla.m
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor.value = calckg(obj, mu_g);  
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor.Attributes.unit = 'Non dimensional'; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor.Attributes.cs = " 341 ";
KG = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor.value;
% -------------------------------------------------------------------------

% ---------------------------------------------------------------------------------------------------------------------------------------
% DENSITY: SEA LEVEL
% ---------------------------------------------------------------------------------------------------------------------------------------
%% x = calcngust(rho0, NormalForceCurveSlope, GustAlleviationFact, GustSpeedCruiseVect, WingLoading, CruiseSpeed, DiveSpeed, FlagToCalc)
% This function is able to calculates in any possible case a vector
% which contains gust load factors value, following CS-VLA
% airworthiness prescription. To have a complete documentation
% check the class file csvla.m

Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise_sea_level.value = calcngust(obj, rho0, ... % Standard atmosphere density
                                                                                            CLalfa, ...               % Normal force curve slope [1/rad]
                                                                                            KG_sea_level, ...         % Gust alleviation factor KG
                                                                                            Ude_cruise, ...           % Gust speed at cruise VC
                                                                                            WS, ...                   % Wing loading in SI units
                                                                                            VC, ...                   % Cruise speed from the V - n diagram 
                                                                                            VD, ...                   % Dive speed from the V - n diagram
                                                                                            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.Attributes.case(1))); % A conveniently defined case switch
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise_sea_level.Attributes.unit = "g's";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise_sea_level.Attributes.cs = " 341 ";
n_gust_cruise_plus_sea_level = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise_sea_level.value;

Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise_sea_level.value = calcngust(obj, rho0, ... % Standard atmosphere density
                                                                                            CLalfa, ...               % Normal force curve slope [1/rad]
                                                                                            KG_sea_level, ...         % Gust alleviation factor KG
                                                                                            Ude_cruise, ...           % Gust speed at cruise V = VC
                                                                                            WS, ...                   % Wing loading in SI units
                                                                                            VC, ...                   % Cruise speed from the V - n diagram
                                                                                            VD, ...                   % Dive speed from the V - n diagram 
                                                                                            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.Attributes.case(2)));  % A conveniently defined case switch
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise_sea_level.Attributes.unit = "g's";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise_sea_level.Attributes.cs = " 341 "; 
n_gust_cruise_neg_sea_level = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise_sea_level.value;

Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive_sea_level.value = calcngust(obj, rho0, ...
                                                                                            CLalfa, ...
                                                                                            KG_sea_level, ...
                                                                                            Ude_dive, ...
                                                                                            WS, ...
                                                                                            VC, ...
                                                                                            VD, ...
                                                                                            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.Attributes.case(1)));
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive_sea_level.Attributes.unit = "g's";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive_sea_level.Attributes.cs = " 341 "; 
n_gust_dive_plus_sea_level = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive_sea_level.value;

Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive_sea_level.value = calcngust(obj, rho0, ...
                                                                                            CLalfa, ...
                                                                                            KG_sea_level, ...
                                                                                            Ude_dive, ...
                                                                                            WS, ...
                                                                                            VC, ...
                                                                                            VD, ...
                                                                                            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.Attributes.case(2)));
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive_sea_level.Attributes.unit = "g's";  
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive_sea_level.Attributes.cs = " 341 "; 
n_gust_dive_neg_sea_level = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive_sea_level.value;
% -------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------

% ---------------------------------------------------------------------------------------------------------------------------------------
% DENSITY: OPERATIVE ALTITUDE
% ---------------------------------------------------------------------------------------------------------------------------------------
%% x = calcngust(rho0, NormalForceCurveSlope, GustAlleviationFact, GustSpeedCruiseVect, WingLoading, CruiseSpeed, DiveSpeed, FlagToCalc)
% This function is able to calculates in any possible case a vector
% which contains gust load factors value, following CS-VLA
% airworthiness prescription. To have a complete documentation
% check the class file csvla.m

Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.value = calcngust(obj, rho0, ... % Standard atmosphere density
                                                                                            CLalfa, ...               % Normal force curve slope [1/rad]
                                                                                            KG, ...                   % Gust alleviation factor KG
                                                                                            Ude_cruise, ...           % Gust speed at cruise VC
                                                                                            WS, ...                   % Wing loading in SI units
                                                                                            VC, ...                   % Cruise speed from the V - n diagram 
                                                                                            VD, ...                   % Dive speed from the V - n diagram
                                                                                            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.Attributes.case(1))); % A conveniently defined case switch
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.Attributes.unit = "g's";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.Attributes.cs = " 341 ";
n_gust_cruise_plus = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.value;

Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise.value = calcngust(obj, rho0, ... % Standard atmosphere density
                                                                                            CLalfa, ...               % Normal force curve slope [1/rad]
                                                                                            KG, ...                   % Gust alleviation factor KG
                                                                                            Ude_cruise, ...           % Gust speed at cruise V = VC
                                                                                            WS, ...                   % Wing loading in SI units
                                                                                            VC, ...                   % Cruise speed from the V - n diagram
                                                                                            VD, ...                   % Dive speed from the V - n diagram 
                                                                                            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.Attributes.case(2)));  % A conveniently defined case switch
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise.Attributes.unit = "g's";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise.Attributes.cs = " 341 "; 
n_gust_cruise_neg = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise.value;

Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive.value = calcngust(obj, rho0, ...
                                                                                            CLalfa, ...
                                                                                            KG, ...
                                                                                            Ude_dive, ...
                                                                                            WS, ...
                                                                                            VC, ...
                                                                                            VD, ...
                                                                                            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.Attributes.case(1)));
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive.Attributes.unit = "g's";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive.Attributes.cs = " 341 "; 
n_gust_dive_plus = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive.value;

Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive.value = calcngust(obj, rho0, ...
                                                                                            CLalfa, ...
                                                                                            KG, ...
                                                                                            Ude_dive, ...
                                                                                            WS, ...
                                                                                            VC, ...
                                                                                            VD, ...
                                                                                            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.Attributes.case(2)));
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive.Attributes.unit = "g's";  
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive.Attributes.cs = " 341 "; 
n_gust_dive_neg = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive.value;
% -------------------------------------------------------------------------

%% GUST ENVELOPE DIAGRAM
gust_envelope = figure; 
hold on; grid on; grid minor;
% ylim([nmin-2.0 nmax+2.0])
% xlim([0 VD+10])
ylim 'padded';
xlim 'padded';
plot(VSpos, npos, ':r', 'LineWidth', 0.25)
plot(VSneg, nneg, ':r', 'LineWidth', 0.25)
plot(V_from1toS, n_from1toS, '-r', 'LineWidth', 1)
plot(V_from1toS_inv, n_from1toS_inv, '-r', 'LineWidth', 1)
plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth', 1)
plot(V_fromAtoC, n_fromAtoC, '-r', 'LineWidth', 1)
plot(V_fromCtoD, n_fromCtoD, '-r', 'LineWidth', 1)
plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth', 1)
plot(V_fromDto0, n_fromDto0, '-r', 'LineWidth', 1)
plot(V_fromGtoF, n_fromGtoF, '-r', 'LineWidth', 1)
plot(V_fromFto0, n_fromFto0, '-r', 'LineWidth', 1)
plot(V_fromStoG, n_fromStoG, '-r', 'LineWidth', 1)
plot(VA, nA, 'k.', 'MarkerSize', 12)
plot(VC, nC, 'k.', 'MarkerSize', 12)
plot(VD, nD, 'k.', 'MarkerSize', 12)
plot(VE, 0.0, 'k.', 'MarkerSize', 12)
plot(VF, nF, 'k.', 'MarkerSize', 12)
plot(VG, nG, 'k.', 'MarkerSize', 12)
plot(VS, nS, 'k.', 'MarkerSize', 12)
plot(VS_inv, nS_inv, 'k.', 'MarkerSize', 12)
% GUST VALUES
plot(V_gust_cruise, n_gust_cruise_plus, '--k', 'LineWidth', 0.25)
plot(V_gust_cruise, n_gust_cruise_neg, '--k', 'LineWidth', 0.25)
plot(V_gust_dive, n_gust_dive_plus, '--k', 'LineWidth', 0.25)
plot(V_gust_dive, n_gust_dive_neg, '--k', 'LineWidth', 0.25)
plot( [V_gust_cruise(end) V_gust_dive(end)], ...
      [n_gust_cruise_plus(end) n_gust_dive_plus(end)], '-.k', 'LineWidth', 0.2)
plot( [V_gust_cruise(end) V_gust_dive(end)], ...
      [n_gust_cruise_neg(end) n_gust_dive_neg(end)], '-.k', 'LineWidth', 0.2)
text(VA, nA, 'Point A', 'FontSize', 6)
text(VC, nC, 'Point C', 'FontSize', 6)
text(VD, nD, 'Point D', 'FontSize', 6)
text(VE, 0.0, 'Point E', 'FontSize', 6)
text(VF, nF, 'Point F', 'FontSize', 6)
text(VG, nG, 'Point G', 'FontSize', 6)
text(VS, nS, 'Point S', 'FontSize', 6)
text(VS_inv, nS_inv, 'Point S invert.', 'FontSize', 6)
xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
title("V~-~n diagram per ", Reg, "Interpreter", "latex")

% EXPORT FIGURE
exportgraphics(gust_envelope, 'Gustenvelope.pdf', 'ContentType', 'vector')
exportgraphics(gust_envelope, 'Gustenvelope.png', 'ContentType', 'vector')

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
b        = rho0 * CLalfa * KG * Ude_cruise;
c        = 2 * WS; 
eqn      = a * V^2 - b * V - c ;
Solution = vpasolve(eqn, V);

% WE MUST DECIDE THE NEW MANOEUVRING AIRSPEED
for i = 1:length(Solution) 
    if Solution(i) > 0 
        new_VA = cast(Solution(i), 'double');
        if VA > new_VA
            
            % CASE 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Straight_flight.value = 'Case 1';
            
            % Manoeuvring speed
            VA1 = VA;
            nA1 = nmax; 
            % Cruise speed
            VC  = VC;
            nC  = nGust(rho0, VC, CLalfa, KG, Ude_cruise, WS);
            
            % FROM 0 TO S
            n_from0toS = linspace(0.0, nS, numb)';
            V_from0toS = VS * ones(numb, 1); 
            
            % FROM S TO A1 
            n_fromStoA1 = linspace(nS, nA1, numb)';
            V_fromStoA1 = Vstall(WS, rho0, CLMAX_clean, n_fromStoA1);
            
            % POINT S
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value    = 'Point S';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.value            = VS; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.nS.value            = 1.0; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.nS.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value          = CLmax_func(rho0, S, VS, WS, 1.0);
            CL_S = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.Attributes.unit = "Non dimensional";
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
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value          = CLmax_func(rho0, S, VA, WS, abs(nA));
            CL_A = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value            = 0.5*rho0*VA^2;
            qA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA.value            = CL_A*qA*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA.Attributes.unit  = "daN";            
            
            % POINT A1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.point_name.value    = 'Point A1';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.VA1.value            = VA1; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.VA1.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.nA1.value            = nA1; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.nA1.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.value          = CLmax_func(rho0, S, VA1, WS, abs(nA1));
            CL_A1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.value            = 0.5*rho0*VA1^2;
            qA1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1.value            = CL_A1*qA1*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LA1.Attributes.unit  = "daN";
            
            if max(n_gust_cruise_plus) > nmax
                % FROM A1 TO C1 
                V_test       = linspace(VA1, VC, numb)';
                n_test       = nGust(rho0, V_test, CLalfa, KG, Ude_cruise, WS);
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
                n_fromC1toC = nGust(rho0, V_fromC1toC, CLalfa, KG, Ude_cruise, WS);
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
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value          = CLmax_func(rho0, S, VC, WS, nC); %CLMAX_clean; % CLmax_func(rho0, S, VC, WS, nC);
                CL_C = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.Attributes.unit = "Non dimensional";
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
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC2.CL_C2.Attributes.unit = "Non dimensional";
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
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.CL_C1.Attributes.unit = "Non dimensional";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.qC1.value            = 0.5*rho0*VC1^2;
                qC1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.qC1.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.qC1.Attributes.unit  = "Pa"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.LC1.value            = CL_C1*qC1*S*1e-1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC1.LC1.Attributes.unit  = "daN";
                
                % POINT D
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value    = 'Point D';
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value            = VD; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value            = nD; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value          = CLmax_func(rho0, S, VD, WS, abs(nD));
                CL_D = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.Attributes.unit = "Non dimensional";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value            = 0.5*rho0*VD^2;
                qD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.Attributes.unit  = "Pa"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD.value            = CL_D*qD*S*1e-1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD.Attributes.unit  = "daN";
                
                % VALUES TO STORE
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_from0toS.value                           = V_from0toS; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_from0toS.Attributes.unit                 = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_from0toS.value                           = n_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_from0toS.Attributes.unit                 = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value                 = V_fromStoA1; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.Attributes.unit       = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value           = n_fromStoA1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.Attributes.unit = "g's";    
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromA1toC1.value                         = V_fromA1toC1; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromA1toC1.Attributes.unit               = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromA1toC1.value                         = n_fromA1toC1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromA1toC1.Attributes.unit               = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromC1toC.value                          = V_fromC1toC; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromC1toC.Attributes.unit                = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromC1toC.value                          = n_fromC1toC;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromC1toC.Attributes.unit                = "g's";        
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromCtoC2.value                          = V_fromCtoC2; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromCtoC2.Attributes.unit                = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromCtoC2.value                          = n_fromCtoC2;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromCtoC2.Attributes.unit                = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromC2toD.value                          = V_fromC2toD; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromC2toD.Attributes.unit                = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromC2toD.value                          = n_fromC2toD;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromC2toD.Attributes.unit                = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromDto0.value                           = V_fromDto0; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromDto0.Attributes.unit                 = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromDto0.value                           = n_fromDto0;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromDto0.Attributes.unit                 = "g's";
                % ----------------------------------------------------------------------------------------------------------------------------------
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromDto0.Attributes.unit                 = "g's";   
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC1.value                        = VC1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC1.Attributes.unit              = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC1.value                  = nmax;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC1.Attributes.unit        = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC2.value                        = VC2;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC2.Attributes.unit              = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC2.value                  = nmax;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC2.Attributes.unit        = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC1.Attributes.unit        = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.value                         = VC;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.Attributes.unit               = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.value                   = nC;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.Attributes.unit         = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VA1.value                        = VA1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VA1.Attributes.unit              = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nA1.value                  = nA1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nA1.Attributes.unit        = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VD.value                         = VD;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VD.Attributes.unit               = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value                   = nD;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.Attributes.unit         = "g's";
                
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
                
            elseif max(n_gust_cruise_plus) < nmax
                
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
                plot(VA, nA, 'k.', 'MarkerSize', 12)
                plot(VS, nS, 'k.', 'MarkerSize', 12)

                text(VD, nD, '  D', 'FontSize', 6)
                text(VC, nC, '  C', 'FontSize', 6)
                text(VA, nA, '  A', 'FontSize', 6) 
                text(VS, nS, '  S', 'FontSize', 6)
                
                % POINT S
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value    = 'Point S';
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.value            = VS; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.nS.value            = nS; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.nS.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value          = CLmax_func(rho0, S, VS, WS, nS);
                CL_S = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.Attributes.unit = "Non dimensional";
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
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.Attributes.unit = "Non dimensional";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value            = 0.5*rho0*VA^2;
                qA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.Attributes.unit  = "Pa"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA.value            = CL_A*qA*S*1e-1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA.Attributes.unit  = "daN";                
                
                % POINT C
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value    = 'Point C';
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value            = VC; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value            = nC; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value          = CLmax_func(rho0, S, VC, WS, nC); % CLMAX_clean; % CLmax_func(rho0, S, VC, WS, nC);
                CL_C = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.Attributes.unit = "Non dimensional";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value            = 0.5*rho0*VC^2;
                qC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.Attributes.unit  = "Pa"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.value            = CL_C*qC*S*1e-1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.Attributes.unit  = "daN";
                
                % POINT D
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value    = 'Point D';
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value            = VD; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value            = nD; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value          = CLmax_func(rho0, S, VD, WS, nD);
                CL_D = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.Attributes.unit = "Non dimensional";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value            = 0.5*rho0*VD^2;
                qD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.Attributes.unit  = "Pa"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD.value            = CL_D*qD*S*1e-1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD.Attributes.unit  = "daN";                
                
                % VALUES TO STORE
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_from0toS.value                           = V_from0toS; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_from0toS.Attributes.unit                 = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_from0toS.value                           = n_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_from0toS.Attributes.unit                 = "g's";                 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value                 = V_fromStoA1; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.Attributes.unit       = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value           = n_fromStoA1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.Attributes.unit = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromA1toC.value                          = V_fromA1toC; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromA1toC.Attributes.unit                = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromA1toC.value                          = n_fromA1toC;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromA1toC.Attributes.unit                = "g's";                 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromCtoD.value                           = V_fromCtoD; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromCtoD.Attributes.unit                 = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromCtoD.value                           = n_fromCtoD;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromCtoD.Attributes.unit                 = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromDto0.value                           = V_fromDto0; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromDto0.Attributes.unit                 = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromDto0.value                           = n_fromDto0;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromDto0.Attributes.unit                 = "g's";                                
                % ----------------------------------------------------------------------------------------------------------------------------
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VA1.value                        = VA1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VA1.Attributes.unit              = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nA1.value                  = nA1;
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nA1.Attributes.unit        = "g's";               
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
                
            end
            
        elseif VA < new_VA
            
            % CASE 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Straight_flight.value = 'Case 2';
            
            VA1 = new_VA;
            nA1 = nGust(rho0, VA1, CLalfa, KG, Ude_cruise, WS);
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
            n_fromStoA1 = linspace(nS, nA1, numb)';
            V_fromStoA1 = Vstall(WS, rho0, CLMAX_clean, n_fromStoA1);

            % FROM A1 TO C
            V_fromA1toC = linspace(VA1, VC, numb)';
            n_fromA1toC = nGust(rho0, V_fromA1toC, CLalfa, KG, Ude_cruise, WS);
            nC          = n_fromA1toC(end);

            % FROM C TO A2 
            p = polyfit([n_gust_cruise_plus(end) n_gust_dive_plus(end)], [V_gust_cruise(end) V_gust_dive(end)], 1);
            n_fromCtoA2 = linspace(n_gust_cruise_plus(end), nmax, numb)';
            V_fromCtoA2 = polyval(p, n_fromCtoA2);

            % FROM FINAL GUST TO D
            n_fromA2toD = nmax * ones(numb, 1); 
            V_fromA2toD = linspace(V_fromCtoA2(end), VD, numb)';
            nA2         = nmax;
            VA2         = V_fromA2toD(1);

            % FROM D TO ZERO 
            V_fromDto0 = VD * ones(numb, 1); 
            n_fromDto0 = linspace(nmax, 0.0, numb)';
            
            % VALUES TO STORE
            % -----------------------------------------------------------------------------------------------------------------------------
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_from0toS.value                           = V_from0toS; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_from0toS.Attributes.unit                 = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_from0toS.value                           = n_from0toS;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_from0toS.Attributes.unit                 = "g's";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value                 = [V_fromStoA1; V_fromA1toC]; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.Attributes.unit       = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value           = [n_fromStoA1; n_fromA1toC];
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.Attributes.unit = "g's";           
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromCtoA2.value                          = V_fromCtoA2; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromCtoA2.Attributes.unit                = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromCtoA2.value                          = n_fromCtoA2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromCtoA2.Attributes.unit                = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromA2toD.value                          = V_fromA2toD; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromA2toD.Attributes.unit                = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromA2toD.value                          = n_fromA2toD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromA2toD.Attributes.unit                = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromDto0.value                           = V_fromDto0; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromDto0.Attributes.unit                 = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromDto0.value                           = n_fromDto0;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromDto0.Attributes.unit                 = "g's";                                          
            % ----------------------------------------------------------------------------------------------------------------------------            
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VA2.value                        = VA2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VA2.Attributes.unit              = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nA2.value                  = nA2;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nA2.Attributes.unit        = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.value                         = VC;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.Attributes.unit               = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.value                   = nC; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.Attributes.unit         = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VD.value                         = VD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VD.Attributes.unit               = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value                   = nD;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.Attributes.unit         = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VA1.value                        = VA1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VA1.Attributes.unit              = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nA1.value                  = nA1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nA1.Attributes.unit        = "g's";
            
            % POINT S
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value    = 'Point S';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.value            = VS; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.nS.value            = nS; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.nS.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value          = CLmax_func(rho0, S, VS, WS, nS);
            CL_S = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.Attributes.unit = "Non dimensional";
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
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value          = CLmax_func(rho0, S, VA, WS, abs(nA));
            CL_A = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value            = 0.5*rho0*VA^2;
            qA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA.value            = CL_A*qA*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA.Attributes.unit  = "daN";              

            % POINT A1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.point_name.value    = 'Point A1';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.VA1.value            = VA1; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.VA1.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.nA1.value            = nA1; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.nA1.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.value          = CLmax_func(rho0, S, VA1, WS, abs(nA1));
            CL_A1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.CL_A1.Attributes.unit = "Non dimensional";
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
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value          = CLmax_func(rho0, S, VC, WS, nC); % CLMAX_clean; % CLmax_func(rho0, S, VC, WS, nC);
            CL_C = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value            = 0.5*rho0*VC^2;
            qC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.value            = CL_C*qC*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.Attributes.unit  = "daN";

            % POINT A2
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.point_name.value    = 'Point A2';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.VA2.value            = VA2; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.VA2.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.nA2.value            = nA2; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.nA2.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CL_A2.value          = CLmax_func(rho0, S, VA2, WS, nA2);
            CL_A2 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CL_A2.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.CL_A2.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.qA2.value            = 0.5*rho0*VA2^2;
            qA2 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.qA2.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.qA2.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.LA2.value            = CL_A2*qA2*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA2.LA2.Attributes.unit  = "daN";
            
            % POINT D
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value    = 'Point D';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value            = VD; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value            = nD; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value          = CLmax_func(rho0, S, VD, WS, nD);
            CL_D = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value            = 0.5*rho0*VD^2;
            qD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD.value            = CL_D*qD*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD.Attributes.unit  = "daN";               
            % -----------------------------------------------------------------------------------------------------------------------------
            
            figure(3);
            hold on; 
            plot(V_fromStoA1, n_fromStoA1, 'r', 'LineWidth', 1.5)
            plot(V_fromA1toC, n_fromA1toC, 'r', 'LineWidth', 1.5)
            plot(V_fromCtoA2, n_fromCtoA2, 'r', 'LineWidth', 1.5)
            plot(V_fromA2toD, n_fromA2toD, 'r', 'LineWidth', 1.5)
            plot(V_fromDto0, n_fromDto0, 'r', 'LineWidth', 1.5)
            plot(V_from0toS, n_from0toS, 'r', 'LineWidth', 1.5)

            plot(VD, nD, 'k.', 'MarkerSize', 12);
            plot(VA2, nA2, 'k.', 'MarkerSize', 12)
            plot(VC, nC, 'k.', 'MarkerSize', 12)
            plot(VA, nA, 'k.', 'MarkerSize', 12)
            plot(VA1, nA1, 'k.', 'MarkerSize', 12)
            plot(VS, nS, 'k.', 'MarkerSize', 12)
            
            text(VD, nD, '  D', 'FontSize', 6)
            text(VA2, nA2, '   A2', 'FontSize', 6)
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
b        = rho0 * CLalfa * KG * Ude_cruise;
c        = 2 * WS; 
eqn      = a * V^2 + b * V - c;  
Solution = vpasolve(eqn, V);
check_s  = isreal(Solution);
disp(" ")
disp(" CHECK: ")
if check_s == 0
    disp(check_s)
    disp(" WARNING: There are no solutions in the Real numbers") 
    
    % CASE 
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Inverted_flight.value = 'Case 1';
    
    if abs(min(n_gust_cruise_neg)) > abs(nmin)
        VG = VG; 
        nG = nmin; 
        % Cruise airspeed
        VF  = VF; 
        nF  = nGust_inverted(rho0, VF, CLalfa, KG, Ude_cruise, WS);

        % FROM 0 TO S_INVERTED 
        nS_inv        = -1.0;
        n_from0toSinv = linspace(0.0, nS_inv, numb)';
        VS_inv        = Vstall(WS, rho0, abs(CLMAX_clean_inverted), abs(nS_inv));
        V_from0toSinv = VS_inv*ones(numb, 1);

        % FROM S_INVERTED TO G
        n_fromSinvtoG = linspace(nS_inv, nG, numb)';
        V_fromSinvtoG = Vstall(WS, rho0, abs(CLMAX_clean_inverted), abs(n_fromSinvtoG));

        % FROM G TO G1 
        V_test = linspace(VG, VF, numb)';
        n_test = nGust_inverted(rho0, V_test, CLalfa, KG, Ude_cruise, WS);
        tol    = 1e-3;
        for i = 1:length(V_test) 
            if abs(abs(n_test(i)) - abs(nmin)) < tol
                row = i;
                VG1 = V_test(row);
            end
        end
        n_fromGtoG1 = nmin * ones(numb, 1);
        V_fromGtoG1 = linspace(VG, VG1, numb)';

        % FROM G1 TO F 
        V_fromG1toF = linspace(VG1, VF, numb)';
        n_fromG1toF = nGust_inverted(rho0, V_fromG1toF, CLalfa, KG, Ude_cruise, WS);
        nG1         = n_fromG1toF(1);

        % FROM F TO G2
        n_fromFtoG2 = linspace(n_fromG1toF(end), nmin, numb)';
        nG2         = n_fromFtoG2(end);
        p           = polyfit([n_gust_cruise_neg(end) n_gust_dive_neg(end)], ...
                             [V_gust_cruise(end) V_gust_dive(end)], 1);
        V_fromFtoG2 = polyval(p, n_fromFtoG2);
        VG2         = V_fromFtoG2(end);   

        % FROM G2 TO E 
        V_fromG2toE = linspace(VG2, VE, numb)';
        n_fromG2toE = nmin * ones(numb, 1);

        % FROM E TO 0
        n_fromEto0 = linspace(nE, 0.0, numb)';
        V_fromEto0 = VE * ones(numb, 1);

        % -----------------------------------------------------------------------------------------------------------------   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_from0toSinv.value                        = V_from0toSinv; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_from0toSinv.Attributes.unit              = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_from0toSinv.value                        = n_from0toSinv; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_from0toSinv.Attributes.unit              = "g's";        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value                 = V_fromSinvtoG; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.Attributes.unit       = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.value           = n_fromSinvtoG; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.Attributes.unit = "g's";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromGtoG1.value                          = V_fromGtoG1; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromGtoG1.Attributes.unit                = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromGtoG1.value                          = n_fromGtoG1; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromGtoG1.Attributes.unit                = "g's"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromG1toF.value                          = V_fromG1toF; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromG1toF.Attributes.unit                = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromG1toF.value                          = n_fromG1toF; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromG1toF.Attributes.unit                = "g's"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromFtoG2.value                          = V_fromFtoG2; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromFtoG2.Attributes.unit                = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromFtoG2.value                          = n_fromFtoG2; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromFtoG2.Attributes.unit                = "g's";  
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromG2toE.value                          = V_fromG2toE; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromG2toE.Attributes.unit                = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromG2toE.value                          = n_fromG2toE; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromG2toE.Attributes.unit                = "g's"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromEto0.value                           = V_fromEto0; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromEto0.Attributes.unit                 = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromEto0.value                           = n_fromEto0; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromEto0.Attributes.unit                 = "g's";          
        % ------------------------------------------------------------------------------------------------------------------
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.value                  = VF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.Attributes.unit        = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.value            = nF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.Attributes.unit  = "g's";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VE.value                  = VE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VE.Attributes.unit        = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value            = nE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.Attributes.unit  = "g's";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VG1.value                 = VG1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VG1.Attributes.unit       = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nG1.value           = nG1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nG1.Attributes.unit = "g's";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VG2.value                 = VG2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VG2.Attributes.unit       = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nG2.value           = nG2;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nG2.Attributes.unit = "g's";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VG.value                  = VG;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VG.Attributes.unit        = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nG.value            = nG;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nG.Attributes.unit  = "g's";        
        % ------------------------------------------------------------------------------------------------------------------

        % POINT S_INV
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.point_name.value        = 'Point S inv';
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.VS_inv.value            = VS_inv; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.VS_inv.Attributes.unit  = "m/s"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.nS_inv.value            = nS_inv; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.nS_inv.Attributes.unit  = "g's"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv.value          = CLmax_func(rho0, S, VS_inv, WS, nS_inv);
        CL_S_inv = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.value            = 0.5*rho0*VS_inv^2;
        qS_inv = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.Attributes.unit  = "Pa"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv.value            = CL_S_inv*qS_inv*S*1e-1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv.Attributes.unit  = "daN";          
        
        % POINT G1
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.point_name.value    = 'Point G1';
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.VG1.value            = VG1; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.VG1.Attributes.unit  = "m/s"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.nG1.value            = nG1; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.nG1.Attributes.unit  = "g's"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CL_G1.value          = CLmax_func(rho0, S, VG1, WS, abs(nG1));
        CL_G1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CL_G1.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CL_G1.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.qG1.value            = 0.5*rho0*VG1^2;
        qG1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.qG1.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.qG1.Attributes.unit  = "Pa"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1.value            = CL_G1*qG1*S*1e-1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1.Attributes.unit  = "daN";
        
        % POINT G2
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.point_name.value    = 'Point G2';
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.VG2.value            = VG2; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.VG2.Attributes.unit  = "m/s"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.nG2.value            = nG2; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.nG2.Attributes.unit  = "g's"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CL_G2.value          = CLmax_func(rho0, S, VG2, WS, abs(nG2));
        CL_G2 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CL_G2.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.CL_G2.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.qG2.value            = 0.5*rho0*VG2^2;
        qG2 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.qG2.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.qG2.Attributes.unit  = "Pa"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.LG2.value            = CL_G2*qG2*S*1e-1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG2.LG2.Attributes.unit  = "daN";      
        
        % POINT F
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value    = 'Point F';
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.value            = VF; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.Attributes.unit  = "m/s"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.value            = nF; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.Attributes.unit  = "g's"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value          = CLmax_func(rho0, S, VF, WS, abs(nF));
        CL_F = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value            = 0.5*rho0*VF^2;
        qF = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.Attributes.unit  = "Pa"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.value            = CL_F*qF*S*1e-1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.Attributes.unit  = "daN";      
        
        % POINT G
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.point_name.value    = 'Point G';
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG.value            = VG; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG.Attributes.unit  = "m/s"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.nG.value            = nG; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.nG.Attributes.unit  = "g's"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value          = CLmax_func(rho0, S, VG, WS, abs(nG));
        CL_G = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.Attributes.unit = "Non dimensional";
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
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value            = 0.5*rho0*VE^2;
        qE = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.Attributes.unit  = "Pa"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE.value            = CL_E*qE*S*1e-1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE.Attributes.unit  = "daN";
        
        figure(3);
        hold on; 
        plot(V_from0toSinv, n_from0toSinv, 'r', 'LineWidth', 1.5)
        plot(V_fromSinvtoG, n_fromSinvtoG, 'r', 'LineWidth', 1.5)
        plot(V_fromGtoG1, n_fromGtoG1, 'r', 'LineWidth', 1.5)
        plot(V_fromG1toF, n_fromG1toF, 'r', 'LineWidth', 1.5)
        plot(V_fromFtoG2, n_fromFtoG2, 'r', 'LineWidth', 1.5)
        plot(V_fromG2toE, n_fromG2toE, 'r', 'LineWidth', 1.5)
        plot(V_fromEto0, n_fromEto0, 'r', 'LineWidth', 1.5)

        plot(VS_inv, nS_inv, 'k.', 'MarkerSize', 12);
        plot(VG1, nG1, 'k.', 'MarkerSize', 12)
        plot(VG2, nG2, 'k.', 'MarkerSize', 12)
        plot(VG, nG, 'k.', 'MarkerSize', 12)
        plot(VE, nE, 'k.', 'MarkerSize', 12)
        plot(VF, nF, 'k.', 'MarkerSize', 12)

        text(VF, nF, '  F', 'FontSize', 6)
        text(VG, nG, '  G', 'FontSize', 6)
        text(VG1, nG1, '  G1', 'FontSize', 6)
        text(VG2, nG2, '  G2', 'FontSize', 6)
        text(VE, nE, '  E', 'FontSize', 6)
        text(VS_inv, nS_inv, '  S inv.', 'FontSize', 6)
        
    elseif abs(min(n_gust_cruise_neg)) < abs(nmin)
        VG = VG; 
        nG = nmin; 
        % Cruise airspeed
        VF  = VF; 
        nF  = nmin;

        % FROM 0 TO S_INVERTED 
        nS_inv        = -1.0;
        n_from0toSinv = linspace(0.0, nS_inv, numb)';
        VS_inv        = Vstall(WS, rho0, abs(CLMAX_clean_inverted), abs(nS_inv));
        V_from0toSinv = VS_inv*ones(numb, 1);

        % FROM S_INVERTED TO G
        n_fromSinvtoG = linspace(nS_inv, nG, numb)';
        V_fromSinvtoG = Vstall(WS, rho0, abs(CLMAX_clean_inverted), abs(n_fromSinvtoG));
        
        % FROM G TO F
        n_fromGtoF = nmin * ones(numb, 1);
        V_fromGtoF = linspace(VG, VF, numb)';   
        
        % FROM F TO E
        n_fromFtoE = nmin * ones(numb, 1);
        V_fromFtoE = linspace(VF, VE, numb)';  
        
        % FROM E TO 0
        n_fromEto0 = linspace(nE, 0.0, numb)';
        V_fromEto0 = VE * ones(numb, 1);   
        
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
        
        % -----------------------------------------------------------------------------------------------------------------    
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_from0toSinv.value                        = V_from0toSinv; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_from0toSinv.Attributes.unit              = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_from0toSinv.value                        = n_from0toSinv; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_from0toSinv.Attributes.unit              = "g's"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value                 = V_fromSinvtoG; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.Attributes.unit       = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.value           = n_fromSinvtoG; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.Attributes.unit = "g's";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromGtoF.value                           = V_fromGtoF; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromGtoF.Attributes.unit                 = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromGtoF.value                           = n_fromGtoF; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromGtoF.Attributes.unit                 = "g's";   
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromFtoE.value                           = V_fromFtoE; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromFtoE.Attributes.unit                 = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromFtoE.value                           = n_fromFtoE; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromFtoE.Attributes.unit                 = "g's";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromEto0.value                           = V_fromEto0; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromEto0.Attributes.unit                 = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromEto0.value                           = n_fromEto0; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromEto0.Attributes.unit                 = "g's";        
        % ------------------------------------------------------------------------------------------------------------------
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.value                  = VF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.Attributes.unit        = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.value            = nF;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.Attributes.unit  = "g's";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VE.value                  = VE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VE.Attributes.unit        = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value            = nE;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.Attributes.unit  = "g's";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VG.value                  = VG;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VG.Attributes.unit        = "m/s";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nG.value            = nG;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nG.Attributes.unit  = "g's";        
        % ------------------------------------------------------------------------------------------------------------------

        % POINT S_INV
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.point_name.value        = 'Point S inv';
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.VS_inv.value            = VS_inv; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.VS_inv.Attributes.unit  = "m/s"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.nS_inv.value            = nS_inv; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.nS_inv.Attributes.unit  = "g's"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv.value          = CLmax_func(rho0, S, VS_inv, WS, nS_inv);
        CL_S_inv = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.value            = 0.5*rho0*VS_inv^2;
        qS_inv = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.Attributes.unit  = "Pa"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv.value            = CL_S_inv*qS_inv*S*1e-1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv.Attributes.unit  = "daN";         
        
        % POINT G
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.point_name.value    = 'Point G';
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG.value            = VG; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG.Attributes.unit  = "m/s"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.nG.value            = nG; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.nG.Attributes.unit  = "g's"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value          = CLmax_func(rho0, S, VG, WS, abs(nG));
        CL_G = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.Attributes.unit = "Non dimensional";
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
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value            = 0.5*rho0*VE^2;
        qE = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.Attributes.unit  = "Pa"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE.value            = CL_E*qE*S*1e-1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE.Attributes.unit  = "daN";        
        
        % POINT F
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value    = 'Point F';
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.value            = VF; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.Attributes.unit  = "m/s"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.value            = nF; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.Attributes.unit  = "g's"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value          = CLmax_func(rho0, S, VF, WS, abs(nF));
        CL_F = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value            = 0.5*rho0*VF^2;
        qF = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.Attributes.unit  = "Pa"; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.value            = CL_F*qF*S*1e-1;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.Attributes.unit  = "daN";        
        
    end
    
elseif check_s == 1
    
    % CASE 
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Inverted_flight.value = 'Case 2';
    
    disp(" Solutions are real numbers")
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
    for i = 1:length(Solution)
        if (Solution(i) > 0) && (Solution(i) > VG)
            new_VG = cast(abs(Solution(i)), 'double');
            VG1    = new_VG;
            nG1    = nGust_inverted(rho0, VG1, CLalfa, KG, Ude_cruise, WS);
            % Cruise airspeed
            VF  = VF; 
            nF  = nGust_inverted(rho0, VF, CLalfa, KG, Ude_cruise, WS);

            % FROM 0 TO S_INVERTED 
            nS_inv        = -1.0;
            n_from0toSinv = linspace(0.0, nS_inv, numb)';
            VS_inv        = Vstall(WS, rho0, abs(CLMAX_clean_inverted), abs(nS_inv));
            V_from0toSinv = VS_inv*ones(numb, 1);

            % FROM S_INVERTED TO G1
            n_fromSinvtoG1 = linspace(nS_inv, nG1, numb)';
            V_fromSinvtoG1 = Vstall(WS, rho0, abs(CLMAX_clean_inverted), abs(n_fromSinvtoG1));

            % FROM G1 TO F 
            V_fromG1toF = linspace(VG1, VF, numb)';
            n_fromG1toF = nGust_inverted(rho0, V_fromG1toF, CLalfa, KG, Ude_cruise, WS);

            % FROM F TO E
            n_fromFtoE = linspace(n_fromG1toF(end), n_gust_dive_neg(end), numb)';
            nE         = n_fromFtoE(end);
            p          = polyfit([n_gust_cruise_neg(end) n_gust_dive_neg(end)], ...
                                 [V_gust_cruise(end) V_gust_dive(end)], 1);
            V_fromFtoE = polyval(p, n_fromFtoE);
            VE         = V_fromFtoE(end);

            % FROM E TO 0
            V_fromEto0 = VE * ones(numb, 1);
            n_fromEto0 = linspace(nE, 0.0, numb)';

            figure(3);
            hold on; 
            plot(V_from0toSinv, n_from0toSinv, 'r', 'LineWidth', 1.5)
            plot(V_fromSinvtoG1, n_fromSinvtoG1, 'r', 'LineWidth', 1.5)
            plot(V_fromG1toF, n_fromG1toF, 'r', 'LineWidth', 1.5)
            plot(V_fromFtoE, n_fromFtoE, 'r', 'LineWidth', 1.5)
            plot(V_fromEto0, n_fromEto0, 'r', 'LineWidth', 1.5)

            plot(VS_inv, nS_inv, 'k.', 'MarkerSize', 12);
            plot(VG1, nG1, 'k.', 'MarkerSize', 12)
            plot(VG, nG, 'k.', 'MarkerSize', 12)
            plot(VE, nE, 'k.', 'MarkerSize', 12)
            plot(VF, nF, 'k.', 'MarkerSize', 12)

            text(VF, nF, '  F', 'FontSize', 6)
            text(VG, nG, '  G', 'FontSize', 6)
            text(VG1, nG1, '  G1', 'FontSize', 6)
            text(VE, nE, '  E', 'FontSize', 6)
            text(VS_inv, nS_inv, '  S inv.', 'FontSize', 6)    

            % -----------------------------------------------------------------------------------------------------------------      
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_from0toSinv.value                        = V_from0toSinv; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_from0toSinv.Attributes.unit              = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_from0toSinv.value                        = n_from0toSinv; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_from0toSinv.Attributes.unit              = "g's";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value                 = [V_fromSinvtoG1; V_fromG1toF]; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.Attributes.unit       = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.value           = [n_fromSinvtoG1; n_fromG1toF]; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.Attributes.unit = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromFtoE.value                           = V_fromFtoE; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromFtoE.Attributes.unit                 = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromFtoE.value                           = n_fromFtoE; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromFtoE.Attributes.unit                 = "g's";   
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromEto0.value                           = V_fromEto0; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromEto0.Attributes.unit                 = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromEto0.value                           = n_fromEto0; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromEto0.Attributes.unit                 = "g's";               
            % ------------------------------------------------------------------------------------------------------------------
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.value                  = VF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.Attributes.unit        = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.value            = nF;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.Attributes.unit  = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VE.value                  = VE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VE.Attributes.unit        = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value            = nE;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.Attributes.unit  = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VG.value                  = VG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VG.Attributes.unit        = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nG.value            = nG;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nG.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VG1.value                 = VG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VG1.Attributes.unit       = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nG1.value           = nG1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nG1.Attributes.unit = "g's";
            % ------------------------------------------------------------------------------------------------------------------            

            % POINT S_INV
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.point_name.value        = 'Point S inv';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.VS_inv.value            = VS_inv; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.VS_inv.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.nS_inv.value            = nS_inv; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.nS_inv.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv.value          = CLmax_func(rho0, S, VS_inv, WS, nS_inv);
            CL_S_inv = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.CL_S_inv.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.value            = 0.5*rho0*VS_inv^2;
            qS_inv = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.qS_inv.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv.value            = CL_S_inv*qS_inv*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.LS_inv.Attributes.unit  = "daN";             
            
            % POINT G
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.point_name.value    = 'Point G';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG.value            = VG; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.nG.value            = nG; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.nG.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value          = CLmax_func(rho0, S, VG, WS, abs(nG));
            CL_G = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.Attributes.unit = "Non dimensional";
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
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value            = 0.5*rho0*VE^2;
            qE = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE.value            = CL_E*qE*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE.Attributes.unit  = "daN";        

            % POINT F
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value    = 'Point F';
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.value            = VF; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.Attributes.unit  = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.value            = nF; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.Attributes.unit  = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value          = CLmax_func(rho0, S, VF, WS, abs(nF));
            CL_F = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.Attributes.unit = "Non dimensional";
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
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.CL_G1.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.qG1.value            = 0.5*rho0*VG1^2;
            qG1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.qG1.value;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.qG1.Attributes.unit  = "Pa"; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1.value            = CL_G1*qG1*S*1e-1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG1.LG1.Attributes.unit  = "daN";            
            
        end
    end
end

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