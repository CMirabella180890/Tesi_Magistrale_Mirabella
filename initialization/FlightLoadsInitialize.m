function Aircraft = FlightLoadsInitialize()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initilization of Aircraft structure %%
% This function initialize the matlab structure funtion to be used in next
% steps. All the variables to be added in main structure, has to be here
% declared before using.
%
% The declared main object is:  Aircraft
% Aircraft contains:
% 1) TLAR       -> Aircraft.TLAR.value = which are the Top Level Aircraft Requirements
% 2) Geometry   -> Aircraft.Geometry.value = which are the main aircraft
% geometrical data
%
% Variables values could be here added if you want neglect statistical
% calculations and use your proper input values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INIZIALIZATION OF AircraftSTRUCT VARIABLE (Just the minimal values)
Aircraft.Certification.Aircraft_Name.value = NaN;
Aircraft.Certification.Aircraft_Name.Attributes.type = NaN;
Aircraft.Certification.Aircraft_Name.Attributes.designer = NaN;
Aircraft.Certification.Regulation.value = NaN;                                      % Selected regulation
Aircraft.Certification.Regulation.Attributes.Date = NaN;                            % Optional field, to keep note of the actual regulation implemented
Aircraft.Certification.Regulation.Attributes.Amendment = NaN;                       % -- 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Attributes.cs = " 321 ";
Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value = NaN;            % FIRST OUTPUT FROM REGULATION: Maximum load factor
Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.Attributes.unit = "g's";  % Load factor are non dimensional number: L = n*W;
Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.Attributes.cs = " 337(a) ";
Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value = NaN;            % FIRST OUTPUT FROM REGULATION: Minimum load factor
Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.Attributes.unit = "g's";  % Load factor are non dimensional number: L = n*W;
Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.Attributes.cs = " 337(b) ";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Max_Continuous_Power_Speed_VH.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Max_Continuous_Power_Speed_VH.Attributes.unit = "m/s";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Min_Design_Cruise_Speed.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Min_Design_Cruise_Speed.Attributes.unit = "m/s";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.Attributes.unit = "Positive g";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.Attributes.unit = "Negative g";
% ---------------------------------------------------------------------------------------------------------
% FLAG FOR SPANWISE AIRLOADS CALCULATIONS
% OPTIONS: 
% -------- 
% 1. OPEN VSP: The code will run open vsp if available 
% 2. SCHRENK: The code will perform airloads calculation with the Schrenk's
%             method, but pitching moment will not be available and a model
%             for the drag coefficiente must be developed (probably a
%             parabolic approximation will be the choice). 
% ---------------------------------------------------------------------------------------------------------
Aircraft.Certification.Regulation.SubpartC.Flightloads.Airload_case.Attributes.case = NaN;
% ---------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_VS.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_VS.Attributes.unit = "m/s";                                                                              
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_VS.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_VS.Attributes.unit = "m/s"; 
% -------------------------------------------------------------------------
Aircraft.Certification.Regulation.SubpartC.Flightloads.Cruise_Speed_VC.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Cruise_Speed_VC.Attributes.unit = "m/s";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.Attributes.unit = "m/s";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VG.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VG.Attributes.unit = "m/s";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_Design_manoeuvring_speed_VE.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_Design_manoeuvring_speed_VE.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.V_n_diagram.value = NaN;      % Store here the V - N diagram
Aircraft.Certification.Regulation.SubpartC.Flightloads.V_n_diagram.Attributes = NaN; % Store here other metadata or additional informations
% -------------------------------------------------------------------------
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.value = NaN; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.Attributes.unit = 'm/s';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.Attributes.case = {'positive_cruise', 'negative_cruise'};
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.value = NaN; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.Attributes.unit = 'm/s';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.Attributes.case = {'positive_cruise', 'negative_cruise'};
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.Attributes.unit = 'g';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise.Attributes.unit = 'g';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive.Attributes.unit = 'g';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive.Attributes.unit = 'g';
% -------------------------------------------------------------------------
% AIRFOIL NAME
% -------------------------------------------------------------------------
Aircraft.Certification.Aerodynamic_data.airfoil_first_panel.value = NaN; 
Aircraft.Certification.Aerodynamic_data.airfoil_first_panel.Attributes.unit = NaN;
Aircraft.Certification.Aerodynamic_data.airfoil_second_panel.value = NaN; 
Aircraft.Certification.Aerodynamic_data.airfoil_second_panel.Attributes.unit = NaN;
Aircraft.Certification.Aerodynamic_data.airfoil_third_panel.value = NaN; 
Aircraft.Certification.Aerodynamic_data.airfoil_third_panel.Attributes.unit = NaN;
Aircraft.Certification.Aerodynamic_data.airfoil_fourth_panel.value = NaN; 
Aircraft.Certification.Aerodynamic_data.airfoil_fourth_panel.Attributes.unit = NaN;
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value = NaN; 
Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.Attributes.unit = "Non dimensional"; 
Aircraft.Certification.Aerodynamic_data.Min_Lift_Coefficient.value = NaN; 
Aircraft.Certification.Aerodynamic_data.Min_Lift_Coefficient.Attributes.unit = "Non dimensional"; 
Aircraft.Certification.Aerodynamic_data.Max_Inverted_Lift_Coefficient.value = 1.0;
Aircraft.Certification.Aerodynamic_data.Max_Inverted_Lift_Coefficient.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value = NaN; 
Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.Attributes.unit = '1/rad';
Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value = NaN; 
Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.Attributes.unit = '1/deg';
Aircraft.Certification.Aerodynamic_data.alpha.value = NaN; % A vector which contains AoA values 
Aircraft.Certification.Aerodynamic_data.alpha.Attributes.unit = "deg";
Aircraft.Certification.Aerodynamic_data.CL.value = NaN;    % A vector which contains CL values 
Aircraft.Certification.Aerodynamic_data.CL.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.CD.value = NaN;    % A vector which contains CD values 
Aircraft.Certification.Aerodynamic_data.CD.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.CM.value = NaN;    % A vector which contains CM values 
Aircraft.Certification.Aerodynamic_data.CM.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.CD_landing_gear.value = NaN;    % A vector which contains CD values 
Aircraft.Certification.Aerodynamic_data.CD_landing_gear.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.CD0.value = NaN; % Zero lift drag coefficient 
Aircraft.Certification.Aerodynamic_data.CD0.Attributes.unit = "Non dimensional"; 
Aircraft.Certification.Aerodynamic_data.e.value = NaN; % Oswaldt efficiency factor
Aircraft.Certification.Aerodynamic_data.e.Attributes.unit = "Non dimensional"; 
Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value = NaN; % Coefficient inside an expression for the CD in polynomial form
Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value = NaN; % Coefficient inside an expressione for the CD in polynomial form
Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.Attributes.unit = "Non dimensional"; 
Aircraft.Certification.Aerodynamic_data.CL0.value = NaN; % Zero lift coefficient 
Aircraft.Certification.Aerodynamic_data.CL0.Attributes.unit = "Non dimensional"; 
Aircraft.Certification.Aerodynamic_data.CM0.value = NaN; % Zero lift pitching mom. coefficient 
Aircraft.Certification.Aerodynamic_data.CM0.Attributes.unit = "Non dimensional"; 
Aircraft.Certification.Aerodynamic_data.CMCL.value = NaN; % Pitching mom. curve slope 
Aircraft.Certification.Aerodynamic_data.CMCL.Attributes.unit = "Non dimensional"; 
Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value = NaN; 
Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.Attributes.unit = '1/deg';
Aircraft.Certification.Aerodynamic_data.CL_star.value = NaN; % Starting from this lift coefficient, CL vs AoA is a non linear curve 
Aircraft.Certification.Aerodynamic_data.CL_star.Attributes.unit = "Non dimensional"; 
Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a.value = NaN; 
Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b.value = NaN; 
Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c.value = NaN;
Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c.Attributes.unit = "Non dimensional"; 
Aircraft.Certification.Aerodynamic_data.Ideal_cl.value = NaN;
Aircraft.Certification.Aerodynamic_data.Ideal_cl.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value = rad2deg(0.083)*1.0;
Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.Attributes.unit = "1/rad";
Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value = 0.3;
Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.Horizontal.damping_factor.value = NaN; 
Aircraft.Certification.Aerodynamic_data.Horizontal.damping_factor.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.Horizontal.tau.value = NaN;
Aircraft.Certification.Aerodynamic_data.Horizontal.tau.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.Horizontal.tau.Attributes.flag = "Conventional";
Aircraft.Certification.Aerodynamic_data.Horizontal.CL_delta_elevator.value = 2.378;
Aircraft.Certification.Aerodynamic_data.Horizontal.CL_delta_elevator.Attributes.unit = "1/rad";
Aircraft.Certification.Aerodynamic_data.Horizontal.CM_q.value = -17.40; 
Aircraft.Certification.Aerodynamic_data.Horizontal.CM_q.Attributes.unit = "1/rad"; 
Aircraft.Certification.Aerodynamic_data.Horizontal.CM_alpha_dot.value = -5.23;
Aircraft.Certification.Aerodynamic_data.Horizontal.CM_alpha_dot.Attributes.unit = "1/rad"; 
Aircraft.Certification.Aerodynamic_data.Horizontal.eta_horizontal.value = 1.0; 
Aircraft.Certification.Aerodynamic_data.Horizontal.eta_horizontal.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.Vertical.a_vt.value = 3.6528578;
Aircraft.Certification.Aerodynamic_data.Vertical.a_vt.Attributes.unit = "1/rad";
Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_takeoff.value = NaN; % 1.9
Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_takeoff.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_landing.value = NaN; % 2.1
Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_landing.Attributes.unit = "Non dimensional";
% -------------------------------------------------------------------------
% AERODYNAMIC HINGE MOMENTS
% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ELEVATOR
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_delta_rad.value = NaN;
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_delta_rad.Attributes.unit = "1/rad";
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_alfa_rad.value = NaN;
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_alfa_rad.Attributes.unit = "1/rad";
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_delta_deg.value = NaN;
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_delta_deg.Attributes.unit = "1/deg";
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_alfa_deg.value = NaN;
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_alfa_deg.Attributes.unit = "1/deg";
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% RUDDER
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_delta_rad.value = NaN; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_delta_rad.Attributes.unit = "1/rad"; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_alfa_rad.value = NaN; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_alfa_rad.Attributes.unit = "1/rad";
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_delta_deg.value = NaN;
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_delta_deg.Attributes.unit = "1/deg";
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_alfa_deg.value = NaN;
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_alfa_deg.Attributes.unit = "1/deg";
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% AILERON
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_delta_rad.value = NaN; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_delta_rad.Attributes.unit = "1/rad"; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_alfa_rad.value = NaN; 
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_alfa_rad.Attributes.unit = "1/rad";
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_delta_deg.value = NaN;
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_delta_deg.Attributes.unit = "1/deg";
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_alfa_deg.value = NaN;
Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_alfa_deg.Attributes.unit = "1/deg";
% -------------------------------------------------------------------------
% GENERAL
% -------------------------------------------------------------------------
Aircraft.Geometry.General.X_cg.value = 0.0;
Aircraft.Geometry.General.X_cg.Attributes.unit = "m";
Aircraft.Geometry.General.xac.value = NaN;
Aircraft.Geometry.General.xac.Attributes.unit = "m"; % Measured from the aircraft nose
Aircraft.Geometry.General.yac.value = NaN;
Aircraft.Geometry.General.yac.Attributes.unit = "m"; % Measured from the aircraft nose
Aircraft.Geometry.General.zac.value = NaN;
Aircraft.Geometry.General.zac.Attributes.unit = "m"; % Measured from the aircraft nose
Aircraft.Geometry.General.xcg.value = NaN;
Aircraft.Geometry.General.xcg.Attributes.unit = "m"; % Measured from the aircraft nose
Aircraft.Geometry.General.ycg.value = NaN;
Aircraft.Geometry.General.ycg.Attributes.unit = "m"; % Measured from the aircraft nose
Aircraft.Geometry.General.zcg.value = NaN;
Aircraft.Geometry.General.zcg.Attributes.unit = "m"; % Measured from the aircraft nose
Aircraft.Geometry.General.XAC_nondim.value = NaN;
Aircraft.Geometry.General.XAC_nondim.Attributes.unit = "Non dimensional"; % xac/M.A.C.
Aircraft.Geometry.General.XCG_nondim.value = NaN;
Aircraft.Geometry.General.XCG_nondim.Attributes.unit = "Non dimensional";
Aircraft.Geometry.General.bcg.value = NaN;
Aircraft.Geometry.General.bcg.Attributes.unit = "m"; % c.g. distance from the Aerodynamic center from the Z - axis
% -------------------------------------------------------------------------
% Aileron
% -------------------------------------------------------------------------
Aircraft.Geometry.Aileron.S.value = NaN; 
Aircraft.Geometry.Aileron.S.Attributes.unit = "m";
Aircraft.Geometry.Aileron.b.value = NaN; 
Aircraft.Geometry.Aileron.b.Attributes.unit = "m";
Aircraft.Geometry.Aileron.ca.value = NaN;
Aircraft.Geometry.Aileron.ca.Attributes.unit = "m";
Aircraft.Geometry.Aileron.cb.value = NaN;
Aircraft.Geometry.Aileron.cb.Attributes.unit = "m";
Aircraft.Geometry.Aileron.y_inner.value = NaN;
Aircraft.Geometry.Aileron.y_inner.Attributes.unit = "m";
Aircraft.Geometry.Aileron.y_outer.value = NaN;
Aircraft.Geometry.Aileron.y_outer.Attributes.unit = "m";
Aircraft.Geometry.Aileron.eta_inner.value = NaN;
Aircraft.Geometry.Aileron.eta_inner.Attributes.unit = "Non dimensional";
Aircraft.Geometry.Aileron.eta_outer.value = NaN; 
Aircraft.Geometry.Aileron.eta_outer.Attributes.unit = "Non dimensional";
% -------------------------------------------------------------------------
% Elevator
% -------------------------------------------------------------------------
Aircraft.Geometry.Elevator.S.value = NaN; 
Aircraft.Geometry.Elevator.S.Attributes.unit = "m^2";
Aircraft.Geometry.Elevator.chord.value = NaN; 
Aircraft.Geometry.Elevator.chord.Attributes.unit = "m^2";
Aircraft.Geometry.Elevator.chord_ratio_ce_c.value = NaN;
Aircraft.Geometry.Elevator.chord_ratio_ce_c.Attributes.unit = "Non dimensional";
Aircraft.Geometry.Elevator.overhang.value = NaN;
Aircraft.Geometry.Elevator.overhang.Attributes.unit = "Non dimensional";
Aircraft.Geometry.Elevator.span_ratio.value = NaN;
Aircraft.Geometry.Elevator.span_ratio.Attributes.unit = "Non dimensional";
Aircraft.Geometry.Elevator.S_hinge.value = NaN;
Aircraft.Geometry.Elevator.S_hinge.Attributes.unit = "m^2";
% -------------------------------------------------------------------------
% Wing
% -------------------------------------------------------------------------
Aircraft.Geometry.Wing.b.value = NaN;       % Wing span m
Aircraft.Geometry.Wing.b.Attributes.unit = 'm';
Aircraft.Geometry.Wing.S.value = NaN;        % Wing span m2
Aircraft.Geometry.Wing.S.Attributes.unit = "m^2";
Aircraft.Geometry.Wing.AR.value = NaN;
Aircraft.Geometry.Wing.taper.value = NaN;    % taper ratio
% -------------------------------------------------------------------------
Aircraft.Geometry.Wing.sweep_first.value = NaN;     % sweep angle 1/4 c deg.
Aircraft.Geometry.Wing.sweep_first.Attributes.unit = 'deg';
Aircraft.Geometry.Wing.sweep_second.value = NaN;     % sweep angle 1/4 c deg.
Aircraft.Geometry.Wing.sweep_second.Attributes.unit = 'deg';
Aircraft.Geometry.Wing.sweep_third.value = NaN;     % sweep angle 1/4 c deg.
Aircraft.Geometry.Wing.sweep_third.Attributes.unit = 'deg';
% -------------------------------------------------------------------------
Aircraft.Geometry.Wing.sweep_location.value = NaN;     
Aircraft.Geometry.Wing.sweep_location.Attributes.unit = 'percentage';
Aircraft.Geometry.Wing.secondary_sweep_location.value = NaN;     
Aircraft.Geometry.Wing.secondary_sweep_location.Attributes.unit = 'percentage';
Aircraft.Geometry.Wing.croot.value = NaN;     % root chord m
Aircraft.Geometry.Wing.croot.Attributes.unit = 'm';
Aircraft.Geometry.Wing.ctip.value = NaN;% tip chord m
Aircraft.Geometry.Wing.ctip.Attributes.unit = 'm';
Aircraft.Geometry.Wing.xle.value = NaN;      % wing leading edge as fraction of fuselage lenght in the simmetry plane
Aircraft.Geometry.Wing.xle.Attributes.unit = '% fuselage length';
Aircraft.Geometry.Wing.yle.value = NaN;      % wing leading edge as fraction of fuselage lenght in the simmetry plane
Aircraft.Geometry.Wing.yle.Attributes.unit = '% fuselage length';
Aircraft.Geometry.Wing.zle.value = NaN;      % wing leading edge as fraction of fuselage lenght in the simmetry plane
Aircraft.Geometry.Wing.zle.Attributes.unit = '% fuselage length';
Aircraft.Geometry.Wing.xtip_le.value = NaN; % leading edge of tip chord in % of fuselage lenght
Aircraft.Geometry.Wing.xtip_le.Attributes.unit = '% fuselage length';
% -------------------------------------------------------------------------
Aircraft.Geometry.Wing.dihedral_first.value = NaN; % geometric dihedral angle at c/4 in deg.
Aircraft.Geometry.Wing.dihedral_first.Attributes.unit = 'deg';
Aircraft.Geometry.Wing.dihedral_second.value = NaN; % geometric dihedral angle at c/4 in deg.
Aircraft.Geometry.Wing.dihedral_second.Attributes.unit = 'deg';
Aircraft.Geometry.Wing.dihedral_third.value = NaN; % geometric dihedral angle at c/4 in deg.
Aircraft.Geometry.Wing.dihedral_third.Attributes.unit = 'deg';
% -------------------------------------------------------------------------
Aircraft.Geometry.Wing.mac.value = NaN;            % mean aerodynamic chord in meters
Aircraft.Geometry.Wing.mac.Attributes.unit = 'm';
Aircraft.Geometry.Wing.xmac.value = NaN;           % x mac coordinate as function of fuselage length
Aircraft.Geometry.Wing.xmac.Attributes.unit = '% fuselage length';
Aircraft.Geometry.Wing.ymac.value = NaN;           % y mac coordinate as function of semispan
Aircraft.Geometry.Wing.ymac.Attributes.unit = '% semispan';
Aircraft.Geometry.Wing.ypos.value = NaN;           % y position of root chord as wing semispan fraction
Aircraft.Geometry.Wing.ypos.Attributes.unit = '% semispan';
Aircraft.Geometry.Wing.zpos.value = NaN;           % z position of root chord as fuselage diameter fraction
Aircraft.Geometry.Wing.zpos.Attributes.unit = '% fuselage diameter'; % +1 high wing; -1 low wing;
Aircraft.Geometry.Wing.camberloc.value = NaN;
Aircraft.Geometry.Wing.camberloc.Attributes.unit = "percentage";
Aircraft.Geometry.Wing.thickchord.value = NaN; 
Aircraft.Geometry.Wing.thickchord.Attributes.unit = "percentage";
Aircraft.Geometry.Wing.type.value = NaN; % <--- Flag could be: Rectangular, With_kinks
Aircraft.Geometry.Wing.type.Attributes.unit = "flag"; 
% -------------------------------------------------------------------------
Aircraft.Geometry.Wing.twist_angle_first.value = NaN;
Aircraft.Geometry.Wing.twist_angle_first.Attributes.unit = 'deg';
Aircraft.Geometry.Wing.twist_angle_second.value = NaN;
Aircraft.Geometry.Wing.twist_angle_second.Attributes.unit = 'deg';
Aircraft.Geometry.Wing.twist_angle_third.value = NaN;
Aircraft.Geometry.Wing.twist_angle_third.Attributes.unit = 'deg';
Aircraft.Geometry.Wing.twist_angle_fourth.value = NaN;
Aircraft.Geometry.Wing.twist_angle_fourth.Attributes.unit = 'deg';
% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Vertical
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Aircraft.Geometry.Vertical.xle.value = NaN; %of fuselage lenght
Aircraft.Geometry.Vertical.xle.Attributes.unit = "% of fuselage length"; 
Aircraft.Geometry.Vertical.croot.value = NaN; %m
Aircraft.Geometry.Vertical.croot.Attributes.unit = "m";
Aircraft.Geometry.Vertical.ctip.value = NaN; %m
Aircraft.Geometry.Vertical.ctip.Attributes.unit = "m";
Aircraft.Geometry.Vertical.xtip_le.value = NaN; %of fuselage lenght
Aircraft.Geometry.Vertical.xtip_le.Attributes.unit = "% of fuselage length"; 
xtip_le_v = Aircraft.Geometry.Vertical.xtip_le.value;
Aircraft.Geometry.Vertical.b.value = NaN; %m
Aircraft.Geometry.Vertical.b.Attributes.unit = "m";
Aircraft.Geometry.Vertical.zpos.value = NaN;
Aircraft.Geometry.Vertical.zpos.Attributes.unit = "% of fuselage diameter";
Aircraft.Geometry.Vertical.S.value = NaN; 
Aircraft.Geometry.Vertical.S.Attributes.unit = "m^2";
Aircraft.Geometry.Vertical.chord.value = NaN; 
Aircraft.Geometry.Vertical.chord.Attributes.unit = "m";
Aircraft.Geometry.Vertical.sweep.value       = NaN; % 20;
Aircraft.Geometry.Vertical.sweep.Attributes.unit = "deg";
Aircraft.Geometry.Vertical.sweeploc.value    = NaN; % 0;
Aircraft.Geometry.Vertical.sweeploc.Attributes.unit = "percentage";
Aircraft.Geometry.Vertical.secsweeploc.value = NaN;
Aircraft.Geometry.Vertical.secsweeploc.Attributes.unit = "percentage";
Aircraft.Geometry.Vertical.dihedral.value       = NaN; % 0;
Aircraft.Geometry.Vertical.dihedral.Attributes.unit = "deg";
Aircraft.Geometry.Vertical.twist.value       = NaN; % 0;
Aircraft.Geometry.Vertical.twist.Attributes.unit = "deg";
Aircraft.Geometry.Vertical.twistloc.value       = NaN; % 0;
Aircraft.Geometry.Vertical.twistloc.Attributes.unit = "deg";
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Rudder
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Aircraft.Geometry.Rudder.S.value = NaN;
Aircraft.Geometry.Rudder.S.Attributes.unit = "m^2";
Aircraft.Geometry.Rudder.chord.value = NaN;
Aircraft.Geometry.Rudder.chord.Attributes.unit = "m";
Aircraft.Geometry.Rudder.chord_ratio_cf_c.value = NaN;
Aircraft.Geometry.Rudder.chord_ratio_cf_c.Attributes.unit = "Non dimensional";
Aircraft.Geometry.Rudder.overhang.value = NaN;
Aircraft.Geometry.Rudder.overhang.Attributes.unit = "Non dimensional";
Aircraft.Geometry.Rudder.span_ratio.value = NaN;
Aircraft.Geometry.Rudder.span_ratio.Attributes.unit = "Non dimensional";
Aircraft.Geometry.Rudder.max_deflection.value = NaN;
Aircraft.Geometry.Rudder.max_deflection.Attributes.unit = "deg";
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% FLAPS
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Aircraft.Geometry.Flaps.S.value = NaN; 
Aircraft.Geometry.Flaps.S.Attributes.unit = "m";
Aircraft.Geometry.Flaps.b.value = NaN; 
Aircraft.Geometry.Flaps.b.Attributes.unit = "m";
Aircraft.Geometry.Flaps.ca.value = NaN;
Aircraft.Geometry.Flaps.ca.Attributes.unit = "m";
Aircraft.Geometry.Flaps.cb.value = NaN;
Aircraft.Geometry.Flaps.cb.Attributes.unit = "m";
Aircraft.Geometry.Flaps.y_inner.value = NaN;
Aircraft.Geometry.Flaps.y_inner.Attributes.unit = "m";
Aircraft.Geometry.Flaps.y_outer.value = NaN;
Aircraft.Geometry.Flaps.y_outer.Attributes.unit = "m";
Aircraft.Geometry.Flaps.eta_inner.value = NaN; % 0.0 percentage
Aircraft.Geometry.Flaps.eta_inner.Attributes.unit = "Non dimensional";
Aircraft.Geometry.Flaps.eta_outer.value = NaN; % 0.6269 percentage
Aircraft.Geometry.Flaps.eta_outer.Attributes.unit = "Non dimensional";
Aircraft.Geometry.Flaps.cf_c_root.value = NaN;
Aircraft.Geometry.Flaps.cf_c_root.Attributes.unit = "Non dimensional";
Aircraft.Geometry.Flaps.cf_c_tip.value = NaN;
Aircraft.Geometry.Flaps.cf_c_tip.Attributes.unit = "Non dimensional";
Aircraft.Geometry.Flaps.croot.value = NaN;
Aircraft.Geometry.Flaps.croot.Attributes.unit = "m";
Aircraft.Geometry.Flaps.ctip.value = NaN;
Aircraft.Geometry.Flaps.ctip.Attributes.unit = "m"; 
Aircraft.Geometry.Flaps.cf.value = NaN;
Aircraft.Geometry.Flaps.cf.Attributes.unit = "m";

% Aircraft.Geometry.Flaps.croot.value = NaN; %m
% Aircraft.Geometry.Flaps.croot.Attributes.unit = "m";
% Aircraft.Geometry.Flaps.ctip.value = NaN; %m
% Aircraft.Geometry.Flaps.ctip.Attributes.unit = "m";
% Aircraft.Geometry.Flaps.xtip_le.value = NaN; %of fuselage lenght
% Aircraft.Geometry.Flaps.xtip_le.Attributes.unit = "% of fuselage length"; 
% Aircraft.Geometry.Flaps.b.value = NaN; %m
% Aircraft.Geometry.Flaps.b.Attributes.unit = "m";
% Aircraft.Geometry.Flaps.zpos.value = NaN; % % of df
% Aircraft.Geometry.Flaps.zpos.Attributes.unit = "% of fuselage diameter";
% Aircraft.Geometry.Flaps.S.value = NaN; 
% Aircraft.Geometry.Flaps.S.Attributes.unit = "m^2";
% Aircraft.Geometry.Flaps.chord.value = NaN; 
% Aircraft.Geometry.Flaps.chord.Attributes.unit = "m";
% Aircraft.Geometry.Flaps.sweep.value       = NaN; % 20;
% Aircraft.Geometry.Flaps.sweep.Attributes.unit = "deg";
% Aircraft.Geometry.Flaps.sweeploc.value    = NaN; % 0;
% Aircraft.Geometry.Flaps.sweeploc.Attributes.unit = "percentage";
% Aircraft.Geometry.Flaps.secsweeploc.value = NaN;
% Aircraft.Geometry.Flaps.secsweeploc.Attributes.unit = "percentage";
% Aircraft.Geometry.Flaps.dihedral.value       = NaN; % 0;
% Aircraft.Geometry.Flaps.dihedral.Attributes.unit = "deg";
% Aircraft.Geometry.Flaps.twist.value       = NaN; % 0;
% Aircraft.Geometry.Flaps.twist.Attributes.unit = "deg";
% Aircraft.Geometry.Flaps.twistloc.value       = NaN; % 0;
% Aircraft.Geometry.Flaps.twistloc.Attributes.unit = "deg";
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------
% Fuselage
% -------------------------------------------------------------------------
Aircraft.Geometry.Fuselage.id   = "Fuselage";
Aircraft.Geometry.Fuselage.type = 'TransportFuse';     % OpenVSP component type
% Here we can select a flag to properly model the tail empennage: 
% - Double fin 
% - Single fin
% - T tail 
% - Others
Aircraft.Geometry.Fuselage.length.value = NaN;
Aircraft.Geometry.Fuselage.length.Attributes.unit = "m"; % length
Aircraft.Geometry.Fuselage.diameter.value = NaN;              % diameter
Aircraft.Geometry.Fuselage.diameter.Attributes.unit = "m";
Aircraft.Geometry.Fuselage.Non_dim_radius_of_gyration.value = 0.34;
Aircraft.Geometry.Fuselage.Non_dim_radius_of_gyration.Attributes.unit = "Non dimensional";
Aircraft.Geometry.Fuselage.Radius_of_gyration.value = Aircraft.Geometry.Fuselage.length.value*Aircraft.Geometry.Fuselage.Non_dim_radius_of_gyration.value*0.5;
Aircraft.Geometry.Fuselage.Radius_of_gyration.Attributes.unit = "m";
% -------------------------------------------------------------------------
% Horizontal
% -------------------------------------------------------------------------
Aircraft.Geometry.Horizontal.S.value = NaN;     % Horizontal span m2
Aircraft.Geometry.Horizontal.S.Attributes.unit = "m^2";
Aircraft.Geometry.Horizontal.l.value = NaN;     % tail arm in meters
Aircraft.Geometry.Horizontal.l.Attributes.unit = 'm';
Aircraft.Geometry.Horizontal.camber.value      = NaN; % [0 0];
Aircraft.Geometry.Horizontal.camber.Attributes.unit = "percentage";
Aircraft.Geometry.Horizontal.camberloc.value  = NaN; % [0.2 0.2];
Aircraft.Geometry.Horizontal.camberloc.Attributes.unit = "percentage";
Aircraft.Geometry.Horizontal.thickchord.value = NaN; % [0.12 0.12];
Aircraft.Geometry.Horizontal.thickchord.Attributes.unit = "percentage";
Aircraft.Geometry.Horizontal.twist.value       = NaN; % [0 0];
Aircraft.Geometry.Horizontal.twist.Attributes.unit = "deg";
Aircraft.Geometry.Horizontal.twistloc.value    = NaN; % [0.25 0.25];
Aircraft.Geometry.Horizontal.twistloc.Attributes.unit = "percentage";
Aircraft.Geometry.Horizontal.xloc0.value       = NaN; % 1.49;
Aircraft.Geometry.Horizontal.xloc0.Attributes.unit = "m";
Aircraft.Geometry.Horizontal.xloc.value        = NaN; % Aircraft.Geometry.Horizontal.Horizontal.xloc0.value + Aircraft.Geometry.Wing.xle.value; % 1.49+1.638;
Aircraft.Geometry.Horizontal.xloc.Attributes.unit = "m";
Aircraft.Geometry.Horizontal.yloc.value        = 0.0;
Aircraft.Geometry.Horizontal.yloc.Attributes.unit = "m";
Aircraft.Geometry.Horizontal.zloc.value        = NaN; % 0.15;
Aircraft.Geometry.Horizontal.zloc.Attributes.unit = "m";
Aircraft.Geometry.Horizontal.xrot.value        = 0.0;
Aircraft.Geometry.Horizontal.xrot.Attributes.unit = "m";
Aircraft.Geometry.Horizontal.yrot.value        = 0.0;
Aircraft.Geometry.Horizontal.yrot.Attributes.unit = "m";
Aircraft.Geometry.Horizontal.zrot.value        = 0.0;
Aircraft.Geometry.Horizontal.zrot.Attributes.unit = "m";
Aircraft.Geometry.Horizontal.b.value           = NaN; % 1.496;
Aircraft.Geometry.Horizontal.b.Attributes.unit = "m";
Aircraft.Geometry.Horizontal.ctip.value        = NaN; % 0.3136;
Aircraft.Geometry.Horizontal.ctip.Attributes.unit = "m";
Aircraft.Geometry.Horizontal.croot.value       = NaN; % 0.3929;
Aircraft.Geometry.Horizontal.sweep.value       = NaN; % 15;
Aircraft.Geometry.Horizontal.sweep.Attributes.unit = "deg";
Aircraft.Geometry.Horizontal.sweeploc.value    = NaN; % 0;
Aircraft.Geometry.Horizontal.sweeploc.Attributes.unit = "percentage";
Aircraft.Geometry.Horizontal.secsweeploc.value = 1.0;
Aircraft.Geometry.Horizontal.secsweeploc.Attributes.unit = "percentage";
Aircraft.Geometry.Horizontal.dihedral.value    = NaN; % 0;
Aircraft.Geometry.Horizontal.dihedral.Attributes.unit = "deg";
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Aircraft.Geometry.Horizontal.Movable.eta_inner.value = NaN;
% Aircraft.Geometry.Horizontal.Movable.eta_inner.Attributes.unit = "percentage";
% Aircraft.Geometry.Horizontal.Movable.eta_outer.value = NaN;
% Aircraft.Geometry.Horizontal.Movable.eta_outer.Attributes.unit = "percentage";
% Aircraft.Geometry.Horizontal.Movable.cf_c_inner.value = NaN;
% Aircraft.Geometry.Horizontal.Movable.cf_c_inner.Attributes.unit = "percentage";
% Aircraft.Geometry.Horizontal.Movable.cf_c_outer.value = NaN;
% Aircraft.Geometry.Horizontal.Movable.cf_c_outer.Attributes.unit = "percentage";
% Aircraft.Geometry.Horizontal.Movable.max_deflection.value = 25.0;
% Aircraft.Geometry.Horizontal.Movable.max_deflection.Attributes.unit = "deg";
% Aircraft.Geometry.Horizontal.Movable.total_deflection_time.value = NaN;
% Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes.unit = "seconds";
% Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes.flag1 = "Commuter";
% Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes.flag2 = "Wheel";
% Una possibile soluzione alternativa a quella trovata qui di seguito
% potrebbe essere quella di creare dei campi come: 
%  Aircraft.Geometry.Elevator. ...
%  Aircraft.Geometry.Rudder. ...
% In questo modo si semplifica la struttura e si tengono ben separate le
% quantit?? di interesse. 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ELEVATOR
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Aircraft.Geometry.Elevator.eta_inner.value = NaN;
Aircraft.Geometry.Elevator.eta_inner.Attributes.unit = "percentage";
Aircraft.Geometry.Elevator.eta_outer.value = NaN;
Aircraft.Geometry.Elevator.eta_outer.Attributes.unit = "percentage";
Aircraft.Geometry.Elevator.cf_c_inner.value = NaN;
Aircraft.Geometry.Elevator.cf_c_inner.Attributes.unit = "percentage";
Aircraft.Geometry.Elevator.cf_c_outer.value = NaN;
Aircraft.Geometry.Elevator.cf_c_outer.Attributes.unit = "percentage";
Aircraft.Geometry.Elevator.max_deflection.value = 25.0;
Aircraft.Geometry.Elevator.max_deflection.Attributes.unit = "deg";
Aircraft.Geometry.Elevator.total_deflection_time.value = NaN;
Aircraft.Geometry.Elevator.total_deflection_time.Attributes.unit = "seconds";
Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 = "Normal";
Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag2 = "Wheel";
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Aircraft.Geometry.Movable.Horizontal.eta_inner.value = NaN;
% Aircraft.Geometry.Movable.Horizontal.eta_inner.Attributes.unit = "percentage";
% Aircraft.Geometry.Movable.Horizontal.eta_outer.value = NaN;
% Aircraft.Geometry.Movable.Horizontal.eta_outer.Attributes.unit = "percentage";
% Aircraft.Geometry.Movable.Horizontal.cf_c_inner.value = NaN;
% Aircraft.Geometry.Movable.Horizontal.cf_c_inner.Attributes.unit = "percentage";
% Aircraft.Geometry.Movable.Horizontal.cf_c_outer.value = NaN;
% Aircraft.Geometry.Movable.Horizontal.cf_c_outer.Attributes.unit = "percentage";
% Aircraft.Geometry.Movable.Horizontal.max_deflection.value = 25.0;
% Aircraft.Geometry.Movable.Horizontal.max_deflection.Attributes.unit = "deg";
% Aircraft.Geometry.Movable.Horizontal.total_deflection_time.value = NaN;
% Aircraft.Geometry.Movable.Horizontal.total_deflection_time.Attributes.unit = "seconds";
% Aircraft.Geometry.Movable.Horizontal.total_deflection_time.Attributes.flag1 = "Commuter";
% Aircraft.Geometry.Movable.Horizontal.total_deflection_time.Attributes.flag2 = "Wheel";
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% AILERON
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% MAX AILERON FLAP DEFLECTION
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Aircraft.Geometry.Aileron.Max_deflection.value = 15.0; 
Aircraft.Geometry.Aileron.Max_deflection.Attributes.unit = "deg";
% ----------------------------------------------------------------------------------
% Aircraft.Geometry.Aileron.Hinge_coefficients.C_h_delta_rad.value = NaN; 
% Aircraft.Geometry.Aileron.Hinge_coefficients.C_h_delta_rad.Attributes.unit = "1/rad";
% Aircraft.Geometry.Aileron.Hinge_coefficients.C_h_alfa_rad.value = NaN; 
% Aircraft.Geometry.Aileron.Hinge_coefficients.C_h_alfa_rad.Attributes.unit = "1/rad";
% Aircraft.Geometry.Aileron.Hinge_coefficients.C_h_delta_deg.value = NaN; 
% Aircraft.Geometry.Aileron.Hinge_coefficients.C_h_delta_deg.Attributes.unit = "1/deg";
% Aircraft.Geometry.Aileron.Hinge_coefficients.C_h_alfa_deg.value = NaN; 
% Aircraft.Geometry.Aileron.Hinge_coefficients.C_h_alfa_deg.Attributes.unit = "1/deg";
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Aircraft.Geometry.Vertical.MAC.value = 0.2335363;
Aircraft.Geometry.Vertical.MAC.Attributes.unit = "m";
Aircraft.Geometry.Vertical.l_vt.value = 1.650;
Aircraft.Geometry.Vertical.l_vt.Attributes.unit = "m";
Aircraft.Geometry.Vertical.empennage_flag.value = "Double fin"; 
Aircraft.Geometry.Vertical.empennage_flag.Attributes.unit = NaN;
Aircraft.Geometry.Vertical.empennage_flag.Attributes.number_of_fin = 2;
% -------------------------------------------------------------------------
% Engine
% -------------------------------------------------------------------------
Aircraft.Geometry.Engine.Primary.Thrust_axes.value = NaN; % Engine Thrust line of action
Aircraft.Geometry.Engine.Primary.Thrust_axes.Attributes.unit = 'm';
% -------------------------------------------------------------------------
% Weight
% -------------------------------------------------------------------------
% Aircraft.Weight.I_Level.W_maxTakeOff.value = NaN;
% Aircraft.Weight.I_Level.W_maxTakeOff.Attributes.unit = 'kg';
% Aircraft.Weight.I_Level.W_OperativeEmpty.value = NaN;
% Aircraft.Weight.I_Level.W_OperativeEmpty.Attributes.unit = 'kg';
% Aircraft.Weight.I_Level.W_Payload.value = NaN;
% Aircraft.Weight.I_Level.W_Payload.Attributes.unit = 'kg';
% Aircraft.Weight.I_Level.W_Fuel.value = NaN;
% Aircraft.Weight.I_Level.W_Fuel.Attributes.unit = 'kg';
% Aircraft.Weight.I_Level.W_Crew.value = NaN;
% Aircraft.Weight.I_Level.W_Crew.Attributes.unit = 'kg';
% Aircraft.Weight.I_Level.IY.value = 100.0;
% Aircraft.Weight.I_Level.IY.Attributes.unit = "kg * m^2";
% -------------------------------------------------------------------------
Aircraft.Certification.ISA_Condition.Sea_level.Altitude.value = NaN;
Aircraft.Certification.ISA_Condition.Sea_level.Altitude.Attributes.unit = "m";
% -------------------------------------------------------------------------
Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.value = NaN;
Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.Attributes.unit = "m";
% -------------------------------------------------------------------------
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.value = NaN;
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.Attributes.unit = "m";
% -------------------------------------------------------------------------
Aircraft.Certification.Performance.I_Level.Wing_loading_Eng.value = NaN;
Aircraft.Certification.Performance.I_Level.Wing_loading_Eng.Attributes.unit = "psf";
Aircraft.Certification.Performance.I_Level.Power_loading_Eng.value = NaN;
Aircraft.Certification.Performance.I_Level.Power_loading_Eng.Attributes.unit = "lb/hp";
Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value = NaN;
Aircraft.Certification.Performance.I_Level.Wing_loading_SI.Attributes.unit = "Pa";
Aircraft.Certification.Performance.I_Level.Power_loading_SI.value = NaN;
Aircraft.Certification.Performance.I_Level.Power_loading_SI.Attributes.unit = "kg/kW";
% -------------------------------------------------------------------------
Aircraft.Constants.g.value = 9.80665;
Aircraft.Constants.g.Attributes.unit = 'm/s^2'; 
% -------------------------------------------------------------------------
Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value = 2.0;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.Attributes.unit = "g's";
%% INTERPOLATION 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Interpolated_Cm.value = NaN; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Interpolated_Cm.Attributes.unit = 'Non dimensional';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Interpolated_Cd.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Interpolated_Cd.Attributes.unit = 'Non dimensional';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.Attributes.unit = "Non dimensional";

%% ENGINE 
Aircraft.Engine.Takeoff.Power.value = NaN; 
Aircraft.Engine.Takeoff.Power.Attributes.unit = "kW"; 
Aircraft.Engine.Takeoff.RPM.value = NaN; 
Aircraft.Engine.Takeoff.RPM.Attributes.unit = "RPM"; 
Aircraft.Engine.Max_Continous.Power.value = NaN; 
Aircraft.Engine.Max_Continous.Power.Attributes.unit = "kW"; 
Aircraft.Engine.Max_Continous.RPM.value = NaN; 
Aircraft.Engine.Max_Continous.RPM.Attributes.unit = "RPM"; 
Aircraft.Engine.Correction_factor.value = NaN; 
Aircraft.Engine.Correction_factor.Attributes.unit = "Non dimensional";
Aircraft.Engine.Limit_side_load.value = NaN;
Aircraft.Engine.Limit_side_load.Attributes.unit = "Non dimensional";
Aircraft.Engine.Engine_mount_mass.value = NaN;
Aircraft.Engine.Engine_mount_mass.Attributes.unit = "kg"; 
Aircraft.Engine.Engine_accessories_mass.value = NaN;
Aircraft.Engine.Engine_accessories_mass.Attributes.unit = "kg"; 
Aircraft.Engine.Propeller_spinner_mass.value = NaN;
Aircraft.Engine.Propeller_spinner_mass.Attributes.unit = "kg"; 
Aircraft.Engine.Propeller_polar_moment.value = NaN; 
Aircraft.Engine.Propeller_polar_moment.Attributes.unit = "kg * m^2";
Aircraft.Engine.Pitch_speed.value = NaN;
Aircraft.Engine.Pitch_speed.Attributes.unit = "rad/s"; 
Aircraft.Engine.Yaw_speed.value = NaN;
Aircraft.Engine.Yaw_speed.Attributes.unit = "rad/s"; 
Aircraft.Engine.Propeller_blade_number.value = NaN;
Aircraft.Engine.Propeller_blade_number.Attributes.unit = "Pure number";
Aircraft.Engine.Engine_normal_load_factor.value = NaN;
Aircraft.Engine.Engine_normal_load_factor.Attributes.unit = "Non dimensional"; 

% ENGINE FLAGS 
% FLAG 1 
% 1st value: FOUR STROKE
% 2nd value: TWO STROKE
Aircraft.Engine.Correction_factor.Attributes.flag1 = NaN;
% FLAG 2: NUMBER OF CYLINDERS 
Aircraft.Engine.Correction_factor.Attributes.flag2 = NaN;
% REDUCTION RATIO: It's the ratio between the engine rpm and the propeller
%                  rpm
Aircraft.Engine.Reduction_ratio.value = NaN;
Aircraft.Engine.Reduction_ratio.Attributes.unit = "Non dimensional";

end

