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
%          designed for loads corresponding to 125 % of the computed hinge
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
%% CALCULATION METHODS APPLICABLE TO HINGE STRENGTH CALCULATION 
%  1. Roskam;
%  2. NACA / McCormick

%% INITIALIZATION OF THE CALCULATION 

% WING TIP AND ROOT CHORD 
c_root = Aircraft.Geometry.Wing.croot.value;
c_tip  = Aircraft.Geometry.Wing.ctip.value;

% AILERON'S LIMIT 
y_inner = Aircraft.Geometry.Aileron.y_inner.value;
y_outer = Aircraft.Geometry.Aileron.y_outer.value; 

% CHORD RATIOS, INNER AND OUTER
ca_c_inner  = Aircraft.Geometry.Aileron.ca_c_root.value; 
ca_c_outer  = Aircraft.Geometry.Aileron.ca_c_tip.value;

% ROOT CHORD AILERON
c_aileron_root = c_root * ca_c_inner;
c_aileron_tip  = c_tip * ca_c_outer;

% AILERON SURFACE
S_aileron = c_aileron_root * ( y_outer - y_inner );

% STORE VALUES INSIDE THE STRUCT VARIABLE
Aircraft.Geometry.Aileron.croot.value = c_aileron_root;
Aircraft.Geometry.Aileron.croot.Attributes.unit = "m"; 
Aircraft.Geometry.Aileron.ctip.value = c_aileron_tip;
Aircraft.Geometry.Aileron.ctip.Attributes.unit = "m";
Aircraft.Geometry.Aileron.S.value = S_aileron;
Aircraft.Geometry.Aileron.S.Attributes.unit = "m^2";

if c_aileron_root == c_aileron_tip
    ca = c_aileron_root;
elseif c_aileron_root ~= c_aileron_tip
    ca = 0.5 * ( c_aileron_root + c_aileron_tip );
end

% AILERON CHORD 
Aircraft.Geometry.Aileron.ca.value = ca;
Aircraft.Geometry.Aileron.ca.Attributes.unit = "m";

% CALCULATIONS OF THE MOMENT ARM
cb = Aircraft.Geometry.Aileron.cb.value; % OVERHANG
cf = ca - cb; 
Aircraft.Geometry.Aileron.cf.value = cf; 
Aircraft.Geometry.Aileron.cf.Attributes.unit = "m";

% MOMENT ARM CALCULATION 
c_ratio    = cb / cf;
moment_arm = ca * (0.25 - c_ratio); 
Aircraft.Geometry.Aileron.moment_arm.value = moment_arm;
Aircraft.Geometry.Aileron.moment_arm.Attributes.unit = "m";

% HINGE MOMENT COEFFICIENTS IN 1/RAD
C_h_delta_rad = Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_delta_rad.value; 
C_h_alfa_rad  = Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_alfa_rad.value; 

% CONVERSION FACTOR
conversion_factor = 180.0 / pi;

% HINGE MOMENT COEFFICIENT IN 1/DEG
C_h_delta_deg = C_h_delta_rad / conversion_factor;
C_h_alfa_deg  = C_h_alfa_rad / conversion_factor;
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_delta_deg.value = C_h_delta_deg; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_delta_deg.Attributes.unit = "1/deg";
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_alfa_deg.value = C_h_alfa_deg; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_alfa_deg.Attributes.unit = "1/deg";

% MAXIMUM DEFLECTION OF THE AILERON (TIMES TWO FOR DIFF. DEFLECTION)
delta_max_deg = 2 * Aircraft.Geometry.Aileron.Max_deflection.value;
 
% MANOUEVRE SPEED
VA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;

% ANGLE OF ATTACK AT MANOEUVRE SPEED 
% alfa_A_deg = real(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.value);

switch (Straight_flight_Case)
    % CASE 1: VA greater than the intercept
    case 'Case 1'
        % MANOUEVRE SPEED
        VA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.VA1.value;
        
        % ANGLE OF ATTACK AT MANOEUVRE SPEED 
        alfa_A_deg = real(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.alfaA1.value);
    case 'Case 2'
        % MANOUEVRE SPEED
        VA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
        
        % ANGLE OF ATTACK AT MANOEUVRE SPEED 
        alfa_A_deg = real(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.value);
end

% TOTAL HINGE MOMENT COEFFICIENT CH = CH_DELTA * DELTA + CH_ALFA * ALFA
C_h_total_deg = C_h_delta_deg * delta_max_deg + C_h_alfa_deg * alfa_A_deg; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_total_deg.value = C_h_total_deg; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_total_deg.Attributes.unit = "1/deg";
C_h_total_rad = C_h_total_deg / conversion_factor;
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_total_rad.value = C_h_total_rad; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_total_rad.Attributes.unit = "1/rad";

% DYNAMIC PRESSURE AT VA
qA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;

% HINGE MOMENT 
HA = C_h_total_deg * qA * S_aileron * cf;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA.value = HA; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA.Attributes.unit = "N * m";

% GRAVITY ACCELERATION
g = Aircraft.Constants.g.value; 

% CONVERSION OF HINGE MOMENT IN KG * M
HA_converted = HA / g;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_converted.value = HA_converted;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_converted.Attributes.unit = "kg * m";

% 125 % FACTOR TO MAGNIFY THE CALCULATED HINGE MOMENTS
HA_125           = HA * 1.25; 
HA_converted_125 = HA_converted * 1.25;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_125.value = HA_125;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_125.Attributes.unit = "N * m";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_converted_125.value = HA_converted_125;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_converted_125.Attributes.unit = "kg * m";

% HA_125 TIMES TWO TO HAVE THE FULL AILERON LOAD 
HA_total_SI        = HA_125 * 2; 
HA_total_converted = HA_converted_125 * 2;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_total.value = HA_total_SI;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_total.Attributes.unit = "N * m";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_total_converted.value = HA_total_converted;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_total_converted.Attributes.unit = "kg * m";

% AILERON LOAD 
aileron_load_SI        = HA_total_SI / moment_arm;
aileron_load_converted = HA_total_converted / moment_arm;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Aileron_load_SI.value = aileron_load_SI;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Aileron_load_SI.Attributes.unit = " N ";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Aileron_load_SI.Attributes.cs = " 395 ";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Aileron_load_converted.value = aileron_load_converted;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Aileron_load_converted.Attributes.unit = " kg ";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Aileron_load_converted.Attributes.cs = " 395 ";

% Total horizontal tail increment
disp(" ++++ A11 CONTROL SURFACE LOADS - HINGE MOMENTS ++++ ")
disp(" ------------------------------------------------------ ")
Total = [ HA, ...
         HA_converted, ...
         HA_125, ...
         HA_total_SI, ...
         HA_total_converted];
disp(" +++++++++++++++++ Hinge Moments - AILERON +++++++++++++++++ ")
format = ' %6.6f       %6.6f                  %6.6f           %6.6f                %6.6f\n';
label  = ' HA [N * m]  HA_converted [kg * m]       HA_125 [N * m]      HA_total_SI [N * m]      HA_total_converted [kg * m]\n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")