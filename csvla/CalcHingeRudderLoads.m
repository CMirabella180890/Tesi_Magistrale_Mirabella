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

% VERTICAL - CHORDS
croot_vertical = Aircraft.Geometry.Vertical.croot.value;
ctip_vertical  = Aircraft.Geometry.Vertical.ctip.value;
cr_c_root      = Aircraft.Geometry.Rudder.cr_c_root.value;
cr_c_tip       = Aircraft.Geometry.Rudder.cr_c_tip.value;
croot_rudder   = cr_c_root * croot_vertical; 
ctip_rudder    = cr_c_tip * ctip_vertical; 

% STORE INSIDE STRUCT VARIABLE
Aircraft.Geometry.Rudder.croot.value = croot_rudder; 
Aircraft.Geometry.Rudder.croot.Attributes.unit = "m"; 
Aircraft.Geometry.Rudder.ctip.value = ctip_rudder;
Aircraft.Geometry.Rudder.ctip.Attributes.unit = "m";

% VERTICAL - SPAN 
eta_inner_rudder = Aircraft.Geometry.Rudder.eta_inner.value;
eta_outer_rudder = Aircraft.Geometry.Rudder.eta_outer.value;
b_vertical       = Aircraft.Geometry.Vertical.b.value;
b_half           = b_vertical / 2;
y_inner_rudder   = b_half * eta_inner_rudder;
y_outer_rudder    = b_half * eta_outer_rudder;
Aircraft.Geometry.Rudder.y_inner.value = y_inner_rudder; 
Aircraft.Geometry.Rudder.y_inner.Attributes.unit = "m";
Aircraft.Geometry.Rudder.y_outer.value = y_outer_rudder;
Aircraft.Geometry.Rudder.y_outer.Attributes.unit = "m"; 

% S_RUDDER CALCULATIONS 
S_rudder = ( ((croot_rudder + ctip_rudder) * (y_outer_rudder - y_inner_rudder)) / 2 );
Aircraft.Geometry.Rudder.S.value = S_rudder;
Aircraft.Geometry.Rudder.S.Attributes.unit = "m^2";

% CHORD 
if croot_rudder == ctip_rudder
    cr = croot_rudder;
elseif croot_rudder ~= ctip_rudder
    cr = 0.5 * ( croot_rudder + ctip_rudder );
end
Aircraft.Geometry.Rudder.chord.value = cr;
Aircraft.Geometry.Rudder.chord.Attributes.unit = "m";

% STORING INSIDE LOCAL VARIABLES
S_vertical            = Aircraft.Geometry.Vertical.S.value;
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
switch (Aircraft.Geometry.Vertical.empennage_flag.value)
    case 'Multiple fin'
        if Aircraft.Geometry.Vertical.empennage_flag.Attributes.number_of_fin ~ NaN 
            n = Aircraft.Geometry.Vertical.empennage_flag.Attributes.number_of_fin;
        end
        HA_newton = qA * C_h_total_deg * (n * S_rudder) * chord_rudder;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder.value = HA_newton; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder.Attributes.unit = "N * m";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder_kg.value = HA_newton/g; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder_kg.Attributes.unit = "kg * m";
    case 'Double fin'
        HA_newton = qA * C_h_total_deg * (2 * S_rudder) * chord_rudder;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder.value = HA_newton; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder.Attributes.unit = "N * m";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder_kg.value = HA_newton/g; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder_kg.Attributes.unit = "kg * m";
    case 'Single fin' 
        HA_newton = qA * C_h_total_deg * S_rudder * chord_rudder;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder.value = HA_newton; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder.Attributes.unit = "N * m";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder_kg.value = HA_newton/g; 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder_kg.Attributes.unit = "kg * m";
end

% 1.25 * HINGE MOMENT 
HA_newton_125 = HA_newton * 1.25; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder_125.value = HA_newton_125; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder_125.Attributes.unit = "N * m";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder_kg_125.value = (HA_newton_125)/(g); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA_rudder_kg_125.Attributes.unit = "kg * m";

% MOMENT ARM 
moment_arm_rudder = chord_rudder * (0.25 - overhang_rudder);
Aircraft.Geometry.Rudder.moment_arm.value = moment_arm_rudder; 
Aircraft.Geometry.Rudder.moment_arm.Attributes.unit = "m";

% TOTAL LOADS ON THE RUDDER 
total_rudder_loads =  HA_newton_125 / moment_arm_rudder; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.total_rudder_loads.value = total_rudder_loads; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.total_rudder_loads.Attributes.unit = "N";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.total_rudder_loads_kg.value = total_rudder_loads / g; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.total_rudder_loads_kg.Attributes.unit = "kg";

% Total horizontal tail increment
disp(" ++++ A11 CONTROL SURFACE LOADS - HINGE MOMENTS ++++ ")
disp(" ------------------------------------------------------ ")
Total = [ HA_newton, ...
         HA_newton/g, ...
         HA_newton_125, ...
         (HA_newton_125)/(g) ];
disp(" +++++++++++++++++ Hinge Moments - RUDDER +++++++++++++++++ ")
format = ' %6.6f       %6.6f                  %6.6f           %6.6f\n';
label  = ' HR [N * m]  HR_converted [kg * m]       HR_125 [N * m]   HR_125_convert [kg * m]\n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")