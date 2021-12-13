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
Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value = NaN;            % FIRST OUTPUT FROM REGULATION: Maximum load factor
Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.Attributes.unit = "g";  % Load factor are non dimensional number: L = n*W;
Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value = NaN;            % FIRST OUTPUT FROM REGULATION: Minimum load factor
Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.Attributes.unit = "g";  % Load factor are non dimensional number: L = n*W;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Max_Continuous_Power_Speed_VH.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Max_Continuous_Power_Speed_VH.Attributes.unit = "m/s";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Min_Design_Cruise_Speed.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Min_Design_Cruise_Speed.Attributes.unit = "m/s";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.Attributes.unit = "Positive g";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.Attributes.unit = "Negative g";
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
Aircraft.Certification.Aerodynamic_data.alpha.Attributes.unit = "degrees";
Aircraft.Certification.Aerodynamic_data.CL.value = NaN;    % A vector which contains CL values 
Aircraft.Certification.Aerodynamic_data.CL.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.CD.value = NaN;    % A vector which contains CD values 
Aircraft.Certification.Aerodynamic_data.CD.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.CM.value = NaN;    % A vector which contains CM values 
Aircraft.Certification.Aerodynamic_data.CM.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.xac.value = NaN;
Aircraft.Certification.Aerodynamic_data.xac.Attributes.unit = "meters"; % Measured from the aircraft nose
Aircraft.Certification.Aerodynamic_data.yac.value = NaN;
Aircraft.Certification.Aerodynamic_data.yac.Attributes.unit = "meters"; % Measured from the aircraft nose
Aircraft.Certification.Aerodynamic_data.zac.value = NaN;
Aircraft.Certification.Aerodynamic_data.zac.Attributes.unit = "meters"; % Measured from the aircraft nose
Aircraft.Certification.Aerodynamic_data.xcg.value = NaN;
Aircraft.Certification.Aerodynamic_data.xcg.Attributes.unit = "meters"; % Measured from the aircraft nose
Aircraft.Certification.Aerodynamic_data.ycg.value = NaN;
Aircraft.Certification.Aerodynamic_data.ycg.Attributes.unit = "meters"; % Measured from the aircraft nose
Aircraft.Certification.Aerodynamic_data.zcg.value = NaN;
Aircraft.Certification.Aerodynamic_data.zcg.Attributes.unit = "meters"; % Measured from the aircraft nose
Aircraft.Certification.Aerodynamic_data.XAC_nondim.value = NaN;
Aircraft.Certification.Aerodynamic_data.XAC_nondim.Attributes.unit = "Non dimensional"; % xac/M.A.C.
Aircraft.Certification.Aerodynamic_data.bcg.value = NaN;
Aircraft.Certification.Aerodynamic_data.bcg.Attributes.unit = "meters"; % c.g. distance from the Aerodynamic center from the Z - axis
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
% Aircraft.Certification.Aerodynamic_data.Horizontal.tau.value = 0.50; 
Aircraft.Certification.Aerodynamic_data.Horizontal.tau.value = NaN;
Aircraft.Certification.Aerodynamic_data.Horizontal.tau.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.Horizontal.tau.Attributes.flag = "Conventional";
Aircraft.Certification.Aerodynamic_data.Horizontal.CL_delta_elevator.value = 2.378;
Aircraft.Certification.Aerodynamic_data.Horizontal.CL_delta_elevator.Attributes.unit = "1/radians";
Aircraft.Certification.Aerodynamic_data.Horizontal.CM_q.value = -17.40; 
Aircraft.Certification.Aerodynamic_data.Horizontal.CM_q.Attributes.unit = "1/radians"; 
Aircraft.Certification.Aerodynamic_data.Horizontal.CM_alpha_dot.value = -5.23;
Aircraft.Certification.Aerodynamic_data.Horizontal.CM_alpha_dot.Attributes.unit = "1/radians"; 
Aircraft.Certification.Aerodynamic_data.Horizontal.eta_horizontal.value = 1.0; 
Aircraft.Certification.Aerodynamic_data.Horizontal.eta_horizontal.Attributes.unit = "Non dimensional";
Aircraft.Certification.Aerodynamic_data.Vertical.a_vt.value = 3.6528578;
Aircraft.Certification.Aerodynamic_data.Vertical.a_vt.Attributes.unit = "1/rad";
Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_flaps.value = 1.9; 
Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_flaps.Attributes.unit = "Non dimensional";
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
Aircraft.Geometry.Aileron.y_iniziale.value = NaN;
Aircraft.Geometry.Aileron.y_iniziale.Attributes.unit = "m";
Aircraft.Geometry.Aileron.y_finale.value = NaN;
Aircraft.Geometry.Aileron.y_finale.Attributes.unit = "m";
Aircraft.Geometry.Aileron.eta_iniziale.value = NaN;
Aircraft.Geometry.Aileron.eta_iniziale.Attributes.unit = "Non dimensional";
Aircraft.Geometry.Aileron.eta_finale.value = NaN; 
Aircraft.Geometry.Aileron.eta_finale.Attributes.unit = "Non dimensional";
% -------------------------------------------------------------------------
% Wing
% -------------------------------------------------------------------------
Aircraft.Geometry.Wing.b.value = NaN;       % Wing span m
Aircraft.Geometry.Wing.b.Attributes.unit = 'meters';
Aircraft.Geometry.Wing.S.value = NaN;        % Wing span m2
Aircraft.Geometry.Wing.S.Attributes.unit = 'square meters';
Aircraft.Geometry.Wing.AR.value = NaN;
Aircraft.Geometry.Wing.taper.value = NaN;    % taper ratio
Aircraft.Geometry.Wing.sweep.value = NaN;     % sweep angle 1/4 c deg.
Aircraft.Geometry.Wing.sweep.Attributes.unit = 'degrees';
Aircraft.Geometry.Wing.sweep_location.value = NaN;     
Aircraft.Geometry.Wing.sweep_location.Attributes.unit = 'percentage';
Aircraft.Geometry.Wing.secondary_sweep_location.value = NaN;     
Aircraft.Geometry.Wing.secondary_sweep_location.Attributes.unit = 'percentage';
Aircraft.Geometry.Wing.croot.value = NaN;     % root chord m
Aircraft.Geometry.Wing.croot.Attributes.unit = 'meters';
Aircraft.Geometry.Wing.ctip.value = NaN;% tip chord m
Aircraft.Geometry.Wing.ctip.Attributes.unit = 'meters';
Aircraft.Geometry.Wing.xle.value = NaN;      % wing leading edge as fraction of fuselage lenght in the simmetry plane
Aircraft.Geometry.Wing.xle.Attributes.unit = '% fuselage length';
Aircraft.Geometry.Wing.yle.value = NaN;      % wing leading edge as fraction of fuselage lenght in the simmetry plane
Aircraft.Geometry.Wing.yle.Attributes.unit = '% fuselage length';
Aircraft.Geometry.Wing.zle.value = NaN;      % wing leading edge as fraction of fuselage lenght in the simmetry plane
Aircraft.Geometry.Wing.zle.Attributes.unit = '% fuselage length';
Aircraft.Geometry.Wing.xtip_le.value = NaN; % leading edge of tip chord in % of fuselage lenght
Aircraft.Geometry.Wing.xtip_le.Attributes.unit = '% fuselage length';
Aircraft.Geometry.Wing.dihedral.value = NaN; % geometric dihedral angle at c/4 in deg.
Aircraft.Geometry.Wing.dihedral.Attributes.unit = 'degrees';
Aircraft.Geometry.Wing.mac.value = NaN;            % mean aerodynamic chord in meters
Aircraft.Geometry.Wing.mac.Attributes.unit = 'meters';
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

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Vertical
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Aircraft.Geometry.Vertical.xle.value = 0.95; %of fuselage lenght
Aircraft.Geometry.Vertical.xle.Attributes.unit = "% of fuselage length"; 
Aircraft.Geometry.Vertical.croot.value = 0.3136; %m
Aircraft.Geometry.Vertical.croot.Attributes.unit = "m";
Aircraft.Geometry.Vertical.ctip.value = 0.1534725; %m
Aircraft.Geometry.Vertical.ctip.Attributes.unit = "m";
Aircraft.Geometry.Vertical.xtip_le.value = 1.0; %of fuselage lenght
Aircraft.Geometry.Vertical.xtip_le.Attributes.unit = "% of fuselage length"; 
xtip_le_v = Aircraft.Geometry.Vertical.xtip_le.value;
Aircraft.Geometry.Vertical.b.value = 0.437502; %m
Aircraft.Geometry.Vertical.b.Attributes.unit = "m";
Aircraft.Geometry.Vertical.zpos.value = 1.0; % % of df
Aircraft.Geometry.Vertical.zpos.Attributes.unit = "% of df";
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------
% Fuselage
% -------------------------------------------------------------------------
Aircraft.Geometry.Fuselage.id   = "Fuselage";
Aircraft.Geometry.Fuselage.type = 'TransportFuse';     % OpenVSP component type
% Here we can select a flag to properly model the tail empennage: 
% - Double fin 
% - Conventional
% - T tail 
% - Others
Aircraft.Geometry.Fuselage.empennage = "Double fin"; 
Aircraft.Geometry.Fuselage.length.value = NaN;
Aircraft.Geometry.Fuselage.length.Attributes.unit = "meters"; % length
Aircraft.Geometry.Fuselage.diameter.value = NaN;              % diameter
Aircraft.Geometry.Fuselage.diameter.Attributes.unit = "meters";
Aircraft.Geometry.Fuselage.Non_dim_radius_of_gyration.value = 0.34;
Aircraft.Geometry.Fuselage.Non_dim_radius_of_gyration.Attributes.unit = "Non dimensional";
Aircraft.Geometry.Fuselage.Radius_of_gyration.value = Aircraft.Geometry.Fuselage.length.value*Aircraft.Geometry.Fuselage.Non_dim_radius_of_gyration.value*0.5;
Aircraft.Geometry.Fuselage.Radius_of_gyration.Attributes.unit = "meters";
% -------------------------------------------------------------------------
% Horizontal
% -------------------------------------------------------------------------
Aircraft.Geometry.Horizontal.S.value = NaN;     % Horizontal span m2
Aircraft.Geometry.Horizontal.S.Attributes.unit = 'square meters';
Aircraft.Geometry.Horizontal.l.value = NaN;     % tail arm in meters
Aircraft.Geometry.Horizontal.l.Attributes.unit = 'meters';
Aircraft.Geometry.Horizontal.camber.value      = NaN; % [0 0];
Aircraft.Geometry.Horizontal.camber.Attributes.unit = "percentage";
Aircraft.Geometry.Horizontal.camberloc.value  = NaN; % [0.2 0.2];
Aircraft.Geometry.Horizontal.camberloc.Attributes.unit = "percentage";
Aircraft.Geometry.Horizontal.thickchord.value = NaN; % [0.12 0.12];
Aircraft.Geometry.Horizontal.thickchord.Attributes.unit = "percentage";
Aircraft.Geometry.Horizontal.twist.value       = NaN; % [0 0];
Aircraft.Geometry.Horizontal.twist.Attributes.unit = "degrees";
Aircraft.Geometry.Horizontal.twistloc.value    = NaN; % [0.25 0.25];
Aircraft.Geometry.Horizontal.twistloc.Attributes.unit = "percentage";
Aircraft.Geometry.Horizontal.xloc0.value       = NaN; % 1.49;
Aircraft.Geometry.Horizontal.xloc0.Attributes.unit = "meters";
Aircraft.Geometry.Horizontal.xloc.value        = NaN; % Aircraft.Geometry.Horizontal.Horizontal.xloc0.value + Aircraft.Geometry.Wing.xle.value; % 1.49+1.638;
Aircraft.Geometry.Horizontal.xloc.Attributes.unit = "meters";
Aircraft.Geometry.Horizontal.yloc.value        = 0.0;
Aircraft.Geometry.Horizontal.yloc.Attributes.unit = "meters";
Aircraft.Geometry.Horizontal.zloc.value        = NaN; % 0.15;
Aircraft.Geometry.Horizontal.zloc.Attributes.unit = "meters";
Aircraft.Geometry.Horizontal.xrot.value        = 0.0;
Aircraft.Geometry.Horizontal.xrot.Attributes.unit = "meters";
Aircraft.Geometry.Horizontal.yrot.value        = 0.0;
Aircraft.Geometry.Horizontal.yrot.Attributes.unit = "meters";
Aircraft.Geometry.Horizontal.zrot.value        = 0.0;
Aircraft.Geometry.Horizontal.zrot.Attributes.unit = "meters";
Aircraft.Geometry.Horizontal.b.value           = NaN; % 1.496;
Aircraft.Geometry.Horizontal.b.Attributes.unit = "meters";
Aircraft.Geometry.Horizontal.ctip.value        = NaN; % 0.3136;
Aircraft.Geometry.Horizontal.ctip.Attributes.unit = "meters";
Aircraft.Geometry.Horizontal.croot.value       = NaN; % 0.3929;
Aircraft.Geometry.Horizontal.sweep.value       = NaN; % 15;
Aircraft.Geometry.Horizontal.sweep.Attributes.unit = "degrees";
Aircraft.Geometry.Horizontal.sweeploc.value    = NaN; % 0;
Aircraft.Geometry.Horizontal.sweeploc.Attributes.unit = "percentage";
Aircraft.Geometry.Horizontal.secsweeploc.value = 1.0;
Aircraft.Geometry.Horizontal.secsweeploc.Attributes.unit = "percentage";
Aircraft.Geometry.Horizontal.dihedral.value    = NaN; % 0;
Aircraft.Geometry.Horizontal.dihedral.Attributes.unit = "degrees";
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
% Aircraft.Geometry.Horizontal.Movable.max_deflection.Attributes.unit = "degrees";
% Aircraft.Geometry.Horizontal.Movable.total_deflection_time.value = NaN;
% Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes.unit = "seconds";
% Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes.flag1 = "Commuter";
% Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes.flag2 = "Wheel";
% Una possibile soluzione alternativa a quella trovata qui di seguito
% potrebbe essere quella di creare dei campi come: 
%  Aircraft.Geometry.Elevator. ...
%  Aircraft.Geometry.Rudder. ...
% In questo modo si semplifica la struttura e si tengono ben separate le
% quantità di interesse. 
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
Aircraft.Geometry.Elevator.max_deflection.Attributes.unit = "degrees";
Aircraft.Geometry.Elevator.total_deflection_time.value = NaN;
Aircraft.Geometry.Elevator.total_deflection_time.Attributes.unit = "seconds";
Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 = "Commuter";
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
% Aircraft.Geometry.Movable.Horizontal.max_deflection.Attributes.unit = "degrees";
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
Aircraft.Geometry.Aileron.Max_deflection.Attributes.unit = "degrees";
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Aircraft.Geometry.Vertical.MAC.value = 0.2335363;
Aircraft.Geometry.Vertical.MAC.Attributes.unit = "m";
Aircraft.Geometry.Vertical.l_vt.value = 1.650;
Aircraft.Geometry.Vertical.l_vt.Attributes.unit = "m";
% -------------------------------------------------------------------------
% Engine
% -------------------------------------------------------------------------
Aircraft.Geometry.Engine.Primary.Thrust_axes.value = NaN; % Engine Thrust line of action
Aircraft.Geometry.Engine.Primary.Thrust_axes.Attributes.unit = 'meters';
% -------------------------------------------------------------------------
% Weight
% -------------------------------------------------------------------------
Aircraft.Weight.I_Level.W_maxTakeOff.value = NaN;
Aircraft.Weight.I_Level.W_maxTakeOff.Attributes.unit = 'kg';
Aircraft.Weight.I_Level.W_OperativeEmpty.value = NaN;
Aircraft.Weight.I_Level.W_OperativeEmpty.Attributes.unit = 'kg';
Aircraft.Weight.I_Level.W_Payload.value = NaN;
Aircraft.Weight.I_Level.W_Payload.Attributes.unit = 'kg';
Aircraft.Weight.I_Level.W_Fuel.value = NaN;
Aircraft.Weight.I_Level.W_Fuel.Attributes.unit = 'kg';
Aircraft.Weight.I_Level.W_Crew.value = NaN;
Aircraft.Weight.I_Level.W_Crew.Attributes.unit = 'kg';
Aircraft.Weight.I_Level.X_cg.value = 0.0;
Aircraft.Weight.I_Level.X_cg.Attributes.unit = 'meters';
Aircraft.Weight.I_Level.IY.value = 100.0;
Aircraft.Weight.I_Level.IY.Attributes.unit = "kg * m^2";
% -------------------------------------------------------------------------
Aircraft.Certification.ISA_Condition.Sea_level.Altitude.value = NaN;
Aircraft.Certification.ISA_Condition.Sea_level.Altitude.Attribute.unit = "m";
% -------------------------------------------------------------------------
Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.value = NaN;
Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.Attribute.unit = "m";
% -------------------------------------------------------------------------
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.value = NaN;
Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.Attribute.unit = "m";
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

end

