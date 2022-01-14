%% CONTROL HINGE AIRLOADS 
%  In this script we want to evaluate the airloads acting on the hinge of
%  the movable surfaces of the aircraft. Hinge moments are extremely hard
%  to predict with theoretical methods, but they are essential in the
%  sizing process of all the movable surfaces hinges and support
%  structures.
%
%% A11 CONTROL SURFACE LOADS 
%  
%  (a) GENERAL. Each control surface load must be determined using the
%      criteria of sub-paragraph (b) of this paragraph and must lie within 
%      the simplified loading of sub-paragraph (c) of this paragraph.
%  
%  (b) LIMIT PILOT FORCES. In each control surface loading condition
%      described in sub-paragraphs (c) to (e) of this paragraph, the
%      airloads on the movable surfaces and the corresponding deflections
%      need not exceed those which could be obtained in flight by employing
%      the maximum limit pilot forces specified in the table in CS - VLA
%      397 (b). If the surface loads are limited by these maximum travel in
%      the direction which would assist the pilot or the deflection must
%      correspond to the maximum degree of 'out of trim' expected at the
%      speed for the condition under consideration. The tab load, however,
%      need not exceed the value specified in Table 2 of this Appendix. 
%
%  (c) SURFACE LOADING CONDITIONS. Each surface loading condition must be
%      investigated as follows:
%
%      (1) Simplified limit surface loadings and distributions for the
%          horizontal tail, vertical tail, aileron, wing flaps and trim 
%          tabs are specified in Table 2 and figures A4 and A5 of this
%          Appendix. If more than one distribution is given, each
%          distribution must be investigated. Figure A4 is limited to use
%          with vertical tails with aspect ratios less than 2.5 and
%          horizontal tails with aspect ratios less than 5 and tail volumes
%          greater than 0.4. 
%  (d) OUTBOARD FINS. Outboard fins must meet the requirements of CS - VLA 
%      445.
%
%  (e) T AND V TAILS. T and V tails must meet the requirements of CS - VLA 
%      427. 
% 
%  (f) SPECIAL DEVICES. Special devices must meet the requirements of CS -
%      VLA 459. 
%
%% A13 CONTROL SYSTEM LOADS 
% 
%  (a) PRIMARY FLIGHT CONTROLS AND SYSTEMS. Each primary flight control and
%      system must be designed as follows:
%
%      (1) The flight control system and its supporting structure must be
%          designed for loads corresponding to 125 % of the computed hing
%          moments of the movable control surface in the conditions
%          prescribed in paragraph A11 of this Appendix. In addition, 
%      
%          (i) The system limit loads need not exceed those that could be
%              produced by the pilot and automatic devices operating the
%              controls; and 
%  
%         (ii) The design must provide a rugged system for service use,
%              including jamming, ground gusts, taxying downwind, control 
%              inertia and friction. 
%   
%      (2) Acceptable maximum and minimum limit pilot forces for elevator,
%          aileron and rudder controls are shown in the table in CS - VLA
%          387 (b). These pilots loads must be assumed to act at the
%          appropriate control grips or pads as they would under flight
%          conditions and to be reacted at the attachments of the control
%          surface horn. 
%
%  (b) DUAL CONTROLS. If there are dual controls, the systems must be 
%      designed for pilots operating in opposition, using individual pilot
%      loads equal to 75 % of those obtained in accordance with
%      sub-paragraph (a) of this paragraph, except that individual pilot
%      loads may not be less than the minimum limit pilot forces shown in
%      the table in CS - VLA 397 (b). 
% 
%  (c) GROUND GUST CONDITIONS. Ground gust conditions must meet the
%      requirements of CS - VLA 415.
% 
%  (d) SECONDARY CONTROLS AND SYSTEMS. Secondary controls and systems must
%      meet the requirements of CS - VLA 405.
%
%% INITIALIZATION OF THE CALCULATION 

% STORING INSIDE LOCAL VARIABLES
S_vertical            = Aircraft.Geometry.Vertical.S.value;
S_rudder              = Aircraft.Geometry.Rudder.S.value;
chord_rudder          = Aircraft.Geometry.Rudder.chord.value;
chord_ratio_cf_c      = Aircraft.Geometry.Rudder.chord_ratio_cf_c.value;
overhang_rudder       = Aircraft.Geometry.Rudder.overhang.value;
span_ratio_rudder     = Aircraft.Geometry.Rudder.span_ratio.value;
max_deflection_rudder = Aircraft.Geometry.Rudder.max_deflection.value;

% GRAVITY ACCELERATION
g = Aircraft.Constants.g.value;

% HINGE MOMENT COEFFICIENTS DERIVATIS IN 1 / RAD 
C_h_delta_rad = Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_delta_rad.value;
C_h_alfa_rad  = Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_alfa_rad.value;  

% CONVERSION FACTOR
conversion_factor = 180.0 / pi;

% HINGE MOMENT COEFFICIENT IN 1/DEG - RUDDER
C_h_delta_deg = C_h_delta_rad / conversion_factor;
C_h_alfa_deg  = C_h_alfa_rad / conversion_factor;
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_delta_deg.value = C_h_delta_deg; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_delta_deg.Attributes.unit = "1/deg";
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_alfa_deg.value = C_h_alfa_deg; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_alfa_deg.Attributes.unit = "1/deg";

% RUDDER GLOBAL ANGLE OF ATTACK
beta_deg = 0;

% TOTAL HINGE MOMENT COEFFICIENT - RUDDER
C_h_total_deg = C_h_delta_deg * max_deflection_rudder + C_h_alfa_deg * beta_deg;
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_total_deg.value = C_h_total_deg; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_total_deg.Attributes.unit = "1/deg";
C_h_total_rad = C_h_total_deg * conversion_factor;
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_total_rad.value = C_h_total_rad; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_total_rad.Attributes.unit = "1/deg";

% DYNAMIC PRESSURE AT VA
qA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;

% HINGE MOMENT IN NEWTON 
HA_newton = qA * C_h_total_deg * (2 * S_rudder) * chord_rudder;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder.value = HA_newton; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder.Attributes.unit = "N * m";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder_kg.value = HA_newton/g; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder_kg.Attributes.unit = "kg * m";

% 1.25 * HINGE MOMENT 
HA_newton_125 = HA_newton * 1.25; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder_125.value = HA_newton_125; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder_125.Attributes.unit = "N * m";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder_kg_125.value = (HA_newton_125)/(g); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder_kg_125.Attributes.unit = "kg * m";

% TAKING INTO ACCOUNT THE TWO FINS
% Total_hinge_moment_125_newton = 2 * HA_newton_125;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Total_hinge_moment_125_newton.value = Total_hinge_moment_125_newton; 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Total_hinge_moment_125_newton.Attributes.unit = "N * m";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Total_hinge_moment_125_kg.value = (Total_hinge_moment_125_newton) / g; 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Total_hinge_moment_125_kg.Attributes.unit = "kg * m";

% MOMENT ARM 
moment_arm_rudder = chord_rudder * (0.25 - overhang_rudder);

% TOTAL LOADS ON THE RUDDER 
total_rudder_loads =  HA_newton_125 / moment_arm_rudder; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.total_rudder_loads.value = total_rudder_loads; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.total_rudder_loads.Attributes.unit = "N * m";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.total_rudder_loads_kg.value = total_rudder_loads / g; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.total_rudder_loads_kg.Attributes.unit = "N * m";
