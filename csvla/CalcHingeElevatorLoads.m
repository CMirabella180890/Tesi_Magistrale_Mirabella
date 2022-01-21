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

% INIZIALIZATION OF CHORDS VALUES 
b_horizontal         = Aircraft.Geometry.Horizontal.b.value;
c_root_horizontal    = Aircraft.Geometry.Horizontal.croot.value;
c_tip_horizontal     = Aircraft.Geometry.Horizontal.ctip.value;
ce_c_root_horizontal = Aircraft.Geometry.Horizontal.ce_c_root.value;
ce_c_tip_horizontal  = Aircraft.Geometry.Horizontal.ce_c_tip.value;
eta_inner_elevator   = Aircraft.Geometry.Elevator.eta_inner.value;
eta_outer_elevator   = Aircraft.Geometry.Elevator.eta_outer.value;
y_inner_elevator     = eta_inner_elevator * ( b_horizontal / 2 );
y_outer_elevator     = eta_outer_elevator * ( b_horizontal / 2 );
c_root_elevator      = ce_c_root_horizontal * c_root_horizontal;
c_tip_elevator       = ce_c_tip_horizontal * c_tip_horizontal;

% STORE INSIDE THE STRUCT
Aircraft.Geometry.Elevator.y_inner.value = y_inner_elevator;
Aircraft.Geometry.Elevator.y_inner.Attributes.unit = "m";
Aircraft.Geometry.Elevator.y_outer.value = y_outer_elevator;
Aircraft.Geometry.Elevator.y_outer.Attributes.unit = "m";

% CHORD 
if c_root_elevator == c_tip_elevator
    ce = c_root_elevator;
elseif c_root_elevator ~= c_tip_elevator
    ce = 0.5 * ( c_root_elevator + c_tip_elevator );
end

% CHORD RATIO
if ce_c_root_horizontal == ce_c_tip_horizontal
    ce_c_elevator = ce_c_root_horizontal;
elseif ce_c_root_horizontal ~= ce_c_tip_horizontal
    ce_c_elevator = 0.5 * ( ce_c_root_horizontal + ce_c_tip_horizontal );
end

% INSIDE STRUCT
Aircraft.Geometry.Elevator.chord.value = ce;
Aircraft.Geometry.Elevator.chord.Attributes.unit = "m";
Aircraft.Geometry.Elevator.chord_ratio_ce_c.value = ce_c_elevator;
Aircraft.Geometry.Elevator.chord_ratio_ce_c.Attributes.unit = "Non dimensional";

% ELEVATOR S_elevator
S_elevator = ( (c_root_elevator + c_tip_elevator) * (y_outer_elevator - y_inner_elevator) / 2) * 2;
Aircraft.Geometry.Elevator.S.value = S_elevator;
Aircraft.Geometry.Elevator.S.Attributes.unit = "m^2";

% ELEVATOR C_f 
chord_elevator       = Aircraft.Geometry.Elevator.chord.value;
chord_ratio_elevator = Aircraft.Geometry.Elevator.chord_ratio_ce_c.value; 
overhang             = Aircraft.Geometry.Elevator.overhang.value;
cf_elevator          = chord_elevator - chord_elevator * overhang;

% STORING THE DATA
Aircraft.Geometry.Elevator.cf.value           = cf_elevator; 
Aircraft.Geometry.Elevator.cf.Attributes.unit = "m";

% ELEVATOR MAX DEFLECTION
elevator_max_deflection = Aircraft.Geometry.Elevator.max_deflection.value;

% MANOUEVRE SPEED
VA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;

% ANGLE OF ATTACK AT MANOEUVRE SPEED 
alfa_A_deg = real(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.value);

% D EPSILON / D ALFA 
depsilon_dalfa = 0.3;

% HINGE MOMENT COEFFICIENTS IN 1/RAD
C_h_delta_rad = Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_delta_rad.value; 
C_h_alfa_rad  = Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_alfa_rad.value; 

% CONVERSION FACTOR
conversion_factor = 180.0 / pi;

% HINGE MOMENT COEFFICIENT IN 1/DEG
C_h_delta_deg = C_h_delta_rad / conversion_factor;
C_h_alfa_deg  = C_h_alfa_rad / conversion_factor;
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_delta_deg.value = C_h_delta_deg; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_delta_deg.Attributes.unit = "1/deg";
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_alfa_deg.value = C_h_alfa_deg; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_alfa_deg.Attributes.unit = "1/deg";

% MAXIMUM DEFLECTION OF THE ELEVATOR (TIMES TWO FOR DIFF. DEFLECTION)
delta_max_deg = Aircraft.Geometry.Elevator.max_deflection.value;

% TOTAL HINGE MOMENT COEFFICIENT CH = CH_DELTA * DELTA + CH_ALFA * ALFA * ( 1 - d EPSILON / d ALFA )
C_h_total_deg = C_h_delta_deg * delta_max_deg + C_h_alfa_deg * alfa_A_deg * (1 - depsilon_dalfa); 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_total_deg.value = C_h_total_deg; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_total_deg.Attributes.unit = "1/deg";
C_h_total_rad = C_h_total_deg * conversion_factor;
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_total_rad.value = C_h_total_rad; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_total_rad.Attributes.unit = "1/deg";

% DYNAMIC PRESSURE AT VA
qA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;

% GRAVITY ACCELERATION
g = Aircraft.Constants.g.value; 

% HINGE MOMENT 
HA           = C_h_total_deg * qA * S_elevator * cf_elevator;
HA_converted = HA/g;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_elevator.value = HA; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_elevator.Attributes.unit = "N * m";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_elevator_kg.value = HA/g; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_elevator_kg.Attributes.unit = "kg * m";

% HA_125 TIMES TWO TO HAVE THE FULL ELEVATOR LOAD 
HA_125           = HA * 1.25;
HA_125_converted = HA_125/g;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_elevator_125.value = HA_125; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_elevator_125.Attributes.unit = "N * m";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_elevator_125_kg.value = HA_125/g; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_elevator_125_kg.Attributes.unit = "kg * m";

% TOTAL HINGE MOMENT 
HA_Total = HA_125 * 2;
HA_Total_converted = HA_Total / g;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_Total.value = HA_Total; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_Total.Attributes.unit = "N * m";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_Total_converted.value = HA_Total_converted; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_Total_converted.Attributes.unit = "kg * m";

% HINGE MOMENT ARM 
elevator_moment_arm = chord_elevator * (0.25 - overhang);
Aircraft.Geometry.Elevator.moment_arm.value = elevator_moment_arm; 
Aircraft.Geometry.Elevator.moment_arm.Attributes.unit = "m";

% TOTAL HINGE MOMENT 
Total_elevator_loads_elevator_kg = (( HA_125 * 2 )/g) / (elevator_moment_arm);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Total_elevator_loads_kg.value = Total_elevator_loads_elevator_kg; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Total_elevator_loads_kg.Attributes.unit = "kg";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Total_elevator_loads_kg.value = Total_elevator_loads_elevator_kg * g; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Total_elevator_loads_kg.Attributes.unit = "N";

% Total horizontal tail increment
disp(" ++++ A11 CONTROL SURFACE LOADS - HINGE MOMENTS ++++ ")
disp(" ------------------------------------------------------ ")
Total = [ HA, ...
         HA_converted, ...
         HA_125, ...
         HA_125_converted, ...
         HA_Total, ...
         HA_Total_converted ];
disp(" +++++++++++++++++ Hinge Moments - ELEVATOR +++++++++++++++++ ")
format = ' %6.6f       %6.6f                  %6.6f           %6.6f           %6.6f           %6.6f\n';
label  = ' HE [N * m]  HE_converted [kg * m]       HE_125 [N * m]   HE_125_convert [kg * m]  HE_Total [N * m]    HE_Total_convert [kg * m]\n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")