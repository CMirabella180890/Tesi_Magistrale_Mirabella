function Aircraft = FromFileCertification_fun(Aircraft,varargin)
%% READ INPUT FILE
%% READING EXCEL
% cd ../_Utilities
fprintf('\n');
fprintf('Input File: ');
if nargin==1
    filename='Aircraft.xlsx';
else
    filename=[varargin{1},'.xlsx'];
end
fprintf(filename);

%reading excel spreadsheet with readtable function with options
opts = spreadsheetImportOptions("NumVariables", 3);
% Specify sheet and range
%opts.Sheet = "Sheet1";
%opts.DataRange = "A2:C25";

% Specify column names and types
opts.VariableNames = ["field", "value", "unit"];
opts.VariableTypes = ["char", "char", "char"];

% Specify variable properties
%opts = setvaropts(opts, ["field", "value", "unit"], "WhitespaceRule", "preserve");
%opts = setvaropts(opts, ["field", "value", "unit"], "EmptyFieldRule", "auto");
Table = readtable(filename, opts, "UseExcel", false);

value=Table(:,2);
unit=Table(:,3);
label=Table{:,1};
fprintf('\n');
% cd ../_PreDesign

VarText={'Aircraft.Certification.Aircraft_Name.value', ...
         'Aircraft.Certification.Aircraft_Name.Attributes.type', ...
         'Aircraft.Certification.Aircraft_Name.Attributes.designer', ...
         'Aircraft.Certification.Regulation.value', ...
         'Aircraft.Certification.Regulation.Attributes.Date', ...
         'Aircraft.Certification.Regulation.Attributes.Amendment', ...
         'Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax', ...
         'Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin', ...
         'Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient', ...
         'Aircraft.Certification.Aerodynamic_data.Min_Lift_Coefficient', ...
         'Aircraft.Certification.Regulation.SubpartC.Flightloads.Max_Continous_Power_Speed_VH', ...
         'Aircraft.Certification.Regulation.SubpartC.Flightloads.Min_Design_Cruise_Speed', ...
         'Aircraft.Certification.Regulation.SubpartC.Gustloads.Gust_speed_cruise', ...
         'Aircraft.Certification.Regulation.SubpartC.Gustloads.Gust_speed_dive', ...
         'Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope', ...
         'Aircraft.Certification.Aerodynamic_data.alpha.value', ...
         'Aircraft.Certification.Aerodynamic_data.CL.value', ...
         'Aircraft.Certification.Aerodynamic_data.CD.value', ...
         'Aircraft.Certification.Aerodynamic_data.CM.value', ...
         'Aircraft.Geometry.Wing.S.value', ...
         'Aircraft.Geometry.Wing.b.value', ...
         'Aircraft.Geometry.Wing.AR.value', ...
         'Aircraft.Weight.I_Level.W_maxTakeOff.value', ...
         'Aircraft.Certification.Aerodynamic_data.XAC_nondim', ...
         'Aircraft.Certification.Aerodynamic_data.XCG_nondim', ...
         'Aircraft.Geometry.Horizontal.S.value', ...
         'Aircraft.Geometry.Horizontal.l.value', ...
         'Aircraft.Geometry.Engine.Primary.Thrust_axes.value', ...
         'Aircraft.Certification.Aerodynamic_data.CD_landing_gear.value', ...
         'Aircraft.Certification.Aerodynamic_data.xcg.value', ...
         'Aircraft.Certification.Aerodynamic_data.ycg.value', ...
         'Aircraft.Certification.Aerodynamic_data.zcg.value', ...
         'Aircraft.Certification.Aerodynamic_data.CD0.value', ...
         'Aircraft.Certification.Aerodynamic_data.e.value', ...
         'Aircraft.Certification.Aerodynamic_data.CL0.value', ...
         'Aircraft.Certification.Aerodynamic_data.CL_star.value', ...
         'Aircraft.Certification.Aerodynamic_data.CM0.value', ...
         'Aircraft.Certification.Aerodynamic_data.CMCL.value', ...
         'Aircraft.Certification.Aerodynamic_data.bcg.value', ...
         'Aircraft.Certification.Aerodynamic_data.CM_landing_gear', ...
         'Aircraft.Geometry.Wing.croot', ...
         'Aircraft.Geometry.Wing.ctip', ...
         'Aircraft.Geometry.Wing.twist_angle', ...   
         'Aircraft.Geometry.Aileron.S', ...
         'Aircraft.Geometry.Aileron.b', ...
         'Aircraft.Geometry.Aileron.ca', ...
         'Aircraft.Geometry.Aileron.cb', ...
         'Aircraft.Geometry.Aileron.y_iniziale', ...
         'Aircraft.Geometry.Aileron.y_finale', ...
         'Aircraft.Geometry.Aileron.eta_iniziale', ...
         'Aircraft.Geometry.Aileron.eta_finale', ...
         'Normal_Force_Curve_Slope_deg', ...
         'Aircraft.Geometry.Fuselage.length', ...
         'Aircraft.Geometry.Fuselage.diameter', ...
         'Aircraft.Certification.Aerodynamic_data.Ideal_cl', ...
         'Aircraft.Geometry.Wing.camberloc', ...
         'Aircraft.Geometry.Wing.thickchord', ...
         'Aircraft.Geometry.Wing.xle', ...
         'Aircraft.Geometry.Wing.yle', ...
         'Aircraft.Geometry.Wing.zle', ...
         'Aircraft.Geometry.Horizontal.camber', ...
         'Aircraft.Geometry.Horizontal.location_of_camber', ...
         'Aircraft.Geometry.Horizontal.thickchord', ...
         'Aircraft.Geometry.Horizontal.twistvalue', ...
         'Aircraft.Geometry.Horizontal.twistlocation', ...
         'Aircraft.Geometry.Horizontal.zero_x', ...
         'Aircraft.Geometry.Horizontal.xloc', ...
         'Aircraft.Geometry.Horizontal.zloc', ...
         'Aircraft.Geometry.Horizontal.b', ...
         'Aircraft.Geometry.Horizontal.ctip', ...
         'Aircraft.Geometry.Horizontal.croot', ...
         'Aircraft.Geometry.Horizontal.sweep', ...
         'Aircraft.Geometry.Horizontal.freccia_loc', ...
         'Aircraft.Geometry.Horizontal.freccia_sec_loc', ...
         'Aircraft.Geometry.Horizontal.dihedral', ...
         'Aircraft.Geometry.Horizontal.Movable.eta_inner', ...
         'Aircraft.Geometry.Horizontal.Movable.eta_outer', ...
         'Aircraft.Geometry.Horizontal.Movable.cf_c_inner', ...
         'Aircraft.Geometry.Horizontal.Movable.cf_c_outer', ...
         'Aircraft.Geometry.Wing.sweep', ...
         'Aircraft.Geometry.Wing.secondary_sweep_location', ...
         'Aircraft.Geometry.Wing.freccia_posizione', ...
         'Aircraft.Geometry.Wing.dihedral', ...
         'Aircraft.Certification.ISA_Condition.Sea_level', ...
         'Operative_ceiling', ...
         'Theoretical_ceiling', ... 
         'Airloads_flag', ...
         'Takeoff_power', ...
         'Takeoff_rpm', ...
         'Max_Continous_power', ...
         'Max_Continous_rpm', ...
         'Correction_factor_flag1', ...
         'Correction_factor_flag2', ...
         'Reduction_ratio', ...
         'Engine.Limit_side_load', ...
         'Engine_mount_mass', ... 
         'Engine_accessories', ...
         'Propeller_spinner', ...
         'Propeller_polarmoment', ...
         'Engine_pitch_speed', ...
         'Engine_yaw_speed', ...
         'Propeller_blade_number', ...
         'Engine_normal_load_factor', ...
         'C_h_delta', ...
         'C_h_alfa', ...
         'CL_MAX_TAKEOFF', ...
         'CL_MAX_LANDING', ...
         'S_elevator', ...
         'c_elevator', ...
         'chord_ratio_elevator', ...
         'elevator_overhang', ...
         'span_ratio_elevator', ...
         'S_hinge_elevator', ...
         'Chdeltaelevator', ...
         'Chalfaelevator', ...
         'S_vertical', ...
         'chord_vertical', ...
         'S_rudder', ...
         'chord_rudder', ...
         'chord_ratio_rudder_cf_c', ...
         'overhang_rudder', ...
         'span_ratio_rudder', ...
         'max_deflection_rudder', ...
         'Chdeltarudder', ...
         'Chalfarudder', ...
         'ca_c_inner', ...
         'ca_c_outer'...
         'c_root_horizontal', ...
         'c_tip_horizontal', ...
         'ce_c_horizontal_root', ...
         'ce_c_horizontal_tip', ...
         'c_root_vertical', ...
         'c_tip_vertical', ...
         'cr_c_root', ...
         'cr_c_tip', ...
         'eta_inner_rudder', ...
         'eta_outer_rudder', ...
         'c_kink_one', ...
         'c_kink_two', ...
         'wing_type', ...
         'vertical_flag', ...
         'panel_span1', ...
         'panel_span2', ...
         'panel_span3', ...
         'horizontal_tail_damping_factor', ...
         'airfoil_first_panel', ...
         'airfoil_sec_panel', ...
         'airfoil_third_panel', ...
         'airfoil_fourth_panel', ...
         'sweep_first', ...
         'sweep_second', ...
         'sweep_third', ...
         'dihedral_first', ...
         'dihedral_second', ...
         'dihedral_third', ...
         'twist_angle_second', ...
         'twist_angle_third', ...
         'twist_angle_fourth', ...
         'vertical_sweep', ...
         'sweep_loc_vert', ...
         'sec_sweeploc_vert', ...
         'vertical_dihedral', ...
         'vertical_twist', ...
         'ver_twist_loc', ...
         'vertical_x_leadingedge', ...
         'vertical_croot', ...
         'vertical_ctip', ...
         'vertical_span', ...
         'vertical_z_position', ...
         'flaps_eta_iniziale', ...
         'flaps_eta_finale', ...
         'flaps_y_iniziale', ...
         'flaps_y_finale', ...
         'cf_c_root_flap', ...
         'croot_flap', ...
         'ctip_flap', ...
         'cf_flap', ...
         'flap_surface_S', ...
         'vertical_yaw_angle', ...
         'verticaloads_dr', ...
         'vertical_CY_vector', ...
         'vertical_CY0', ...
         'maximum_delta_rudder', ...
         'CY_vertical_tailplane', ...
         'case3_CY_vert_tp', ...
         'vert_tp_K'};

Index=zeros(1,length(VarText));
for i=1:length(VarText)
    Index(i) = find(~cellfun(@isempty, strfind(label,VarText{i})));
end
LengthIndex = nnz(Index);
if LengthIndex==length(VarText)
    p = strcmp('Aircraft.Certification.Aircraft_Name.value',label);
    Aircraft.Certification.Aircraft_Name.value = char(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aircraft_Name.Attributes.type',label);
    Aircraft.Certification.Aircraft_Name.Attributes.type = char(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aircraft_Name.Attributes.designer',label);
    Aircraft.Certification.Aircraft_Name.Attributes.designer = char(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Regulation.value',label);
    Aircraft.Certification.Regulation.value = char(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Regulation.Attributes.Date',label);
    Aircraft.Certification.Regulation.Attributes.Date = char(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Regulation.Attributes.Amendment',label);
    Aircraft.Certification.Regulation.Attributes.Amendment = char(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax',label);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value = str2double(table2array(value(p==1,1))); 
    Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin',label);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value = str2double(table2array(value(p==1,1))); 
    Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.Attributes.unit = char(table2array(unit(p==1,1)));             
    p = strcmp('Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient',label);
    Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value = str2double(table2array(value(p==1,1))); 
    Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.Attributes.unit = char(table2array(unit(p==1,1)));          
    p = strcmp('Aircraft.Certification.Aerodynamic_data.Min_Lift_Coefficient',label);
    Aircraft.Certification.Aerodynamic_data.Min_Lift_Coefficient.value = str2double(table2array(value(p==1,1))); 
    Aircraft.Certification.Aerodynamic_data.Min_Lift_Coefficient.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('Aircraft.Certification.Regulation.SubpartC.Flightloads.Max_Continous_Power_Speed_VH',label);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Max_Continous_Power_Speed_VH.value = str2double(table2array(value(p==1,1))); 
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Max_Continous_Power_Speed_VH.Attributes.unit = char(table2array(unit(p==1,1)));      
    p = strcmp('Aircraft.Certification.Regulation.SubpartC.Flightloads.Min_Design_Cruise_Speed',label);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Min_Design_Cruise_Speed.value = str2double(table2array(value(p==1,1))); 
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Min_Design_Cruise_Speed.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Wing.mac',label);
%     Aircraft.Geometry.Wing.mac.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Wing.mac.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('Aircraft.Certification.Regulation.SubpartC.Gustloads.Gust_speed_cruise',label);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.value = str2double(table2array(value(p==1,1)));
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Certification.Regulation.SubpartC.Gustloads.Gust_speed_dive',label);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.value = str2double(table2array(value(p==1,1)));
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.Attributes.unit = char(table2array(unit(p==1,1)));   
    p = strcmp('Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope',label);
    Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value = str2double(table2array(value(p==1,1)));
    Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('Aircraft.Certification.Aerodynamic_data.alpha.value',label);
    Aircraft.Certification.Aerodynamic_data.alpha.value = char(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.CL.value',label);
    Aircraft.Certification.Aerodynamic_data.CL.value = char(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.CD.value',label);
    Aircraft.Certification.Aerodynamic_data.CD.value = char(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.CM.value',label);
    Aircraft.Certification.Aerodynamic_data.CM.value = char(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Geometry.Wing.S.value', label);
    Aircraft.Geometry.Wing.S.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Geometry.Wing.b.value', label);
    Aircraft.Geometry.Wing.b.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Geometry.Wing.AR.value', label);
    Aircraft.Geometry.Wing.AR.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Weight.I_Level.W_maxTakeOff.value', label);
    Aircraft.Weight.I_Level.W_maxTakeOff.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.xac.value', label);
%     Aircraft.Geometry.General.xac.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.yac.value', label);
%     Aircraft.Geometry.General.yac.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.zac.value', label);
%     Aircraft.Geometry.General.zac.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.XAC_nondim', label);
    Aircraft.Geometry.General.XAC_nondim.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.XCG_nondim', label);
    Aircraft.Geometry.General.XCG_nondim.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.S.value', label);
    Aircraft.Geometry.Horizontal.S.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.l.value', label);
    Aircraft.Geometry.Horizontal.l.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Geometry.Engine.Primary.Thrust_axes.value', label);
    Aircraft.Geometry.Engine.Primary.Thrust_axes.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.CD_landing_gear.value', label);
    Aircraft.Certification.Aerodynamic_data.CD_landing_gear.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.xcg.value', label);
    Aircraft.Geometry.General.xcg.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.ycg.value', label);
    Aircraft.Geometry.General.ycg.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.zcg.value', label);
    Aircraft.Geometry.General.zcg.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.CD0.value', label);
    Aircraft.Certification.Aerodynamic_data.CD0.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.e.value', label);
    Aircraft.Certification.Aerodynamic_data.e.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.CL0.value', label);
    Aircraft.Certification.Aerodynamic_data.CL0.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.CL_star.value', label);
    Aircraft.Certification.Aerodynamic_data.CL_star.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.CM0.value', label);
    Aircraft.Certification.Aerodynamic_data.CM0.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.CMCL.value', label);
    Aircraft.Certification.Aerodynamic_data.CMCL.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.bcg.value', label);
    Aircraft.Geometry.General.bcg.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.General.bcg.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('Aircraft.Certification.Aerodynamic_data.CM_landing_gear', label);
    Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value = str2double(table2array(value(p==1,1)));
    Aircraft.Certification.Aerodynamic_data.CM_landing_gear.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('Aircraft.Geometry.Wing.croot', label);
    Aircraft.Geometry.Wing.croot.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Wing.croot.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Wing.ctip', label);
    Aircraft.Geometry.Wing.ctip.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Wing.ctip.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Wing.twist_angle', label);
    Aircraft.Geometry.Wing.twist_angle_first.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Wing.twist_angle_first.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Aileron.S', label); 
    Aircraft.Geometry.Aileron.S.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Aileron.S.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Aileron.b', label);
    Aircraft.Geometry.Aileron.b.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Aileron.b.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Aileron.ca', label);
    Aircraft.Geometry.Aileron.ca.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Aileron.ca.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Aileron.cb', label);
    Aircraft.Geometry.Aileron.cb.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Aileron.ca.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Aileron.y_iniziale', label);
    Aircraft.Geometry.Aileron.y_inner.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Aileron.y_inner.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Aileron.y_finale', label);
    Aircraft.Geometry.Aileron.y_outer.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Aileron.y_outer.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Aileron.eta_iniziale', label); 
    Aircraft.Geometry.Aileron.eta_inner.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Aileron.eta_inner.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Aileron.eta_finale', label);
    Aircraft.Geometry.Aileron.eta_outer.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Aileron.eta_outer.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Normal_Force_Curve_Slope_deg', label);
    Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value = str2double(table2array(value(p==1,1)));
    Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Fuselage.length', label);   
    Aircraft.Geometry.Fuselage.length.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Fuselage.length.Attributes.unit = char(table2array(unit(p==1,1)));  
    p = strcmp('Aircraft.Geometry.Fuselage.diameter', label);
    Aircraft.Geometry.Fuselage.diameter.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Fuselage.diameter.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.Ideal_cl', label);
    Aircraft.Certification.Aerodynamic_data.Ideal_cl.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Certification.Aerodynamic_data.Ideal_cl.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Wing.camberloc', label);
    Aircraft.Geometry.Wing.camberloc.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Wing.camberloc.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Wing.thickchord', label);    
    Aircraft.Geometry.Wing.thickchord.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Wing.thickchord.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Wing.xle', label); 
    Aircraft.Geometry.Wing.xle.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Wing.xle.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Wing.yle', label); 
    Aircraft.Geometry.Wing.yle.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Wing.yle.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Wing.zle', label); 
    Aircraft.Geometry.Wing.zle.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Wing.zle.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.camber', label);
    Aircraft.Geometry.Horizontal.camber.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Horizontal.camber.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.location_of_camber', label);
    Aircraft.Geometry.Horizontal.location_of_camber.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Horizontal.location_of_camber.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.thickchord', label);
    Aircraft.Geometry.Horizontal.thickchord.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Horizontal.thickchord.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.twistvalue', label);
    Aircraft.Geometry.Horizontal.twist.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Horizontal.twist.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.twistlocation', label);
    Aircraft.Geometry.Horizontal.twistloc.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Horizontal.twistloc.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.zero_x', label);
    Aircraft.Geometry.Horizontal.xloc0.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Horizontal.xloc0.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.xloc', label);
    Aircraft.Geometry.Horizontal.xloc.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Horizontal.xloc.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.zloc', label);
    Aircraft.Geometry.Horizontal.zloc.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Horizontal.zloc.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.b', label);
    Aircraft.Geometry.Horizontal.b.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Horizontal.b.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.ctip', label);
    Aircraft.Geometry.Horizontal.ctip.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Horizontal.ctip.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.croot', label);
    Aircraft.Geometry.Horizontal.croot.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Horizontal.croot.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.sweep', label);
    Aircraft.Geometry.Horizontal.sweep.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Horizontal.sweep.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.freccia_loc', label);
    Aircraft.Geometry.Horizontal.sweeploc.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Horizontal.sweeploc.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.freccia_sec_loc', label);
    Aircraft.Geometry.Horizontal.secondary_sweep_location.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Horizontal.secondary_sweep_location.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.dihedral', label);
    Aircraft.Geometry.Horizontal.dihedral.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Horizontal.dihedral.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.Movable.eta_inner', label);
    Aircraft.Geometry.Elevator.eta_inner.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Elevator.eta_inner.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.Movable.eta_outer', label);
    Aircraft.Geometry.Elevator.eta_outer.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Elevator.eta_outer.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.Movable.cf_c_inner', label);
    Aircraft.Geometry.Elevator.cf_c_inner.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Elevator.cf_c_inner.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.Movable.cf_c_outer', label);
    Aircraft.Geometry.Elevator.cf_c_outer.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Elevator.cf_c_outer.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Wing.sweep', label);
    Aircraft.Geometry.Wing.sweep.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Wing.sweep.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Wing.secondary_sweep_location', label);
    Aircraft.Geometry.Wing.secondary_sweep_location.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Wing.secondary_sweep_location.Attributes.unit = char(table2array(unit(p==1,1)));  
    p = strcmp('Aircraft.Geometry.Wing.freccia_posizione', label);
    Aircraft.Geometry.Wing.sweep_location.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Wing.sweep_location.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Wing.dihedral', label);
    Aircraft.Geometry.Wing.dihedral.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Wing.dihedral.Attributes.unit = char(table2array(unit(p==1,1)));  
    p = strcmp('Aircraft.Certification.ISA_Condition.Sea_level', label);
    Aircraft.Certification.ISA_Condition.Sea_level.Altitude.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Certification.ISA_Condition.Sea_level.Altitude.Attributes.unit = char(table2array(unit(p==1,1)));  
    p = strcmp('Operative_ceiling', label);
    Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('Theoretical_ceiling', label);
    Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Airloads_flag', label);
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Airload_case.Attributes.case = char(table2array(value(p==1,1)));
    p = strcmp('Takeoff_power', label); 
    Aircraft.Engine.Takeoff.Power.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Engine.Takeoff.Power.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Takeoff_rpm', label);
    Aircraft.Engine.Takeoff.RPM.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Engine.Takeoff.RPM.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Max_Continous_power', label); 
    Aircraft.Engine.Max_Continous.Power.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Engine.Max_Continous.Power.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Max_Continous_rpm', label); 
    Aircraft.Engine.Max_Continous.RPM.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Engine.Max_Continous.RPM.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Correction_factor_flag1', label);
    Aircraft.Engine.Correction_factor.Attributes.flag1 = char(table2array(value(p==1, 1)));
    p = strcmp('Correction_factor_flag2', label);
    Aircraft.Engine.Correction_factor.Attributes.flag2 = str2double(table2array(value(p==1, 1)));
    p = strcmp('Reduction_ratio', label);
    Aircraft.Engine.Reduction_ratio.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Engine.Reduction_ratio.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Engine.Limit_side_load', label);
    Aircraft.Engine.Limit_side_load.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Engine.Limit_side_load.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Engine_mount_mass', label);
    Aircraft.Engine.Engine_mount_mass.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Engine.Engine_mount_mass.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Engine_accessories', label);
    Aircraft.Engine.Engine_accessories_mass.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Engine.Engine_accessories_mass.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Propeller_spinner', label);
    Aircraft.Engine.Propeller_spinner_mass.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Engine.Propeller_spinner_mass.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Propeller_polarmoment', label);
    Aircraft.Engine.Propeller_polar_moment.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Engine.Propeller_polar_moment.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Engine_pitch_speed', label);
    Aircraft.Engine.Pitch_speed.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Engine.Pitch_speed.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('Engine_yaw_speed', label);
    Aircraft.Engine.Yaw_speed.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Engine.Yaw_speed.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('Propeller_blade_number', label);
    Aircraft.Engine.Propeller_blade_number.value = str2double(table2array(value(p==1, 1)));
    p = strcmp('Engine_normal_load_factor', label);
    Aircraft.Engine.Engine_normal_load_factor.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Engine.Engine_normal_load_factor.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('C_h_delta', label);
    Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_delta_rad.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_delta_rad.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('C_h_alfa', label);
    Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_alfa_rad.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_alfa_rad.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('CL_MAX_TAKEOFF', label);
    Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_takeoff.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_takeoff.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('CL_MAX_LANDING', label);
    Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_landing.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_landing.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('S_elevator', label);
    Aircraft.Geometry.Elevator.S.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Elevator.S.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('c_elevator', label);
    Aircraft.Geometry.Elevator.chord.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Elevator.chord.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('chord_ratio_elevator', label);
    Aircraft.Geometry.Elevator.chord_ratio_ce_c.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Elevator.chord_ratio_ce_c.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('elevator_overhang', label);
    Aircraft.Geometry.Elevator.overhang.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Elevator.overhang.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('span_ratio_elevator', label);
    Aircraft.Geometry.Elevator.span_ratio.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Elevator.span_ratio.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('S_hinge_elevator', label);
    Aircraft.Geometry.Elevator.S_hinge.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Elevator.S_hinge.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Chdeltaelevator', label);
    Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_delta_rad.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_delta_rad.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Chalfaelevator', label);
    Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_alfa_rad.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_alfa_rad.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('S_vertical', label);
    Aircraft.Geometry.Vertical.S.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Vertical.S.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('chord_vertical', label);
    Aircraft.Geometry.Vertical.chord.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Vertical.chord.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('S_rudder', label);
    Aircraft.Geometry.Rudder.S.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Rudder.S.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('chord_rudder', label);
    Aircraft.Geometry.Rudder.chord.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Rudder.chord.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('chord_ratio_rudder_cf_c', label);
    Aircraft.Geometry.Rudder.chord_ratio_cf_c.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Rudder.chord_ratio_cf_c.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('overhang_rudder', label);
    Aircraft.Geometry.Rudder.overhang.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Rudder.overhang.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('span_ratio_rudder', label);
    Aircraft.Geometry.Rudder.span_ratio.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Rudder.span_ratio.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('max_deflection_rudder', label);
    Aircraft.Geometry.Rudder.max_deflection.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Rudder.max_deflection.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Chdeltarudder', label);
    Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_delta_rad.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_delta_rad.Attributes.unit = char(table2array(unit(p==1,1)));  
    p = strcmp('Chalfarudder', label);
    Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_alfa_rad.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_alfa_rad.Attributes.unit = char(table2array(unit(p==1,1)));    
    p = strcmp('ca_c_inner', label);
    Aircraft.Geometry.Aileron.ca_c_root.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Aileron.ca_c_root.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('ca_c_outer', label);
    Aircraft.Geometry.Aileron.ca_c_tip.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Aileron.ca_c_tip.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('c_root_horizontal', label);
    Aircraft.Geometry.Horizontal.croot.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Horizontal.croot.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('c_tip_horizontal', label);
    Aircraft.Geometry.Horizontal.ctip.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Horizontal.ctip.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('ce_c_horizontal_root', label);
    Aircraft.Geometry.Horizontal.ce_c_root.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Horizontal.ce_c_root.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('ce_c_horizontal_tip', label);
    Aircraft.Geometry.Horizontal.ce_c_tip.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Horizontal.ce_c_tip.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('c_root_vertical', label);
    Aircraft.Geometry.Vertical.croot.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Vertical.croot.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('c_tip_vertical', label);
    Aircraft.Geometry.Vertical.ctip.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Vertical.ctip.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('cr_c_root', label);
    Aircraft.Geometry.Rudder.cr_c_root.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Rudder.cr_c_root.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('cr_c_tip', label);
    Aircraft.Geometry.Rudder.cr_c_tip.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Rudder.cr_c_tip.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('eta_inner_rudder', label);
    Aircraft.Geometry.Rudder.eta_inner.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Rudder.eta_inner.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('eta_outer_rudder', label);
    Aircraft.Geometry.Rudder.eta_outer.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Rudder.eta_outer.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('c_kink_one', label);
    Aircraft.Geometry.Wing.chord_kink_one.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Wing.chord_kink_one.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('c_kink_two', label);
    Aircraft.Geometry.Wing.chord_kink_two.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Wing.chord_kink_two.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('wing_type', label);
    Aircraft.Geometry.Wing.type.value = char(table2array(value(p==1,1)));
    p = strcmp('vertical_flag', label);
    Aircraft.Geometry.Vertical.empennage_flag.value = char(table2array(value(p==1,1)));
    p = strcmp('panel_span1', label);
    Aircraft.Geometry.Wing.panel_span1.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Wing.panel_span1.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('panel_span2', label);
    Aircraft.Geometry.Wing.panel_span2.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Wing.panel_span2.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('panel_span3', label);
    Aircraft.Geometry.Wing.panel_span3.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Wing.panel_span3.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('horizontal_tail_damping_factor', label);
    Aircraft.Certification.Aerodynamic_data.Horizontal.damping_factor.value = str2double(table2array(value(p==1, 1)));  
    Aircraft.Certification.Aerodynamic_data.Horizontal.damping_factor.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('airfoil_first_panel', label);
    Aircraft.Certification.Aerodynamic_data.airfoil_first_panel.value = char(table2array(value(p==1, 1)));
    Aircraft.Certification.Aerodynamic_data.airfoil_first_panel.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('airfoil_sec_panel', label);
    Aircraft.Certification.Aerodynamic_data.airfoil_second_panel.value = char(table2array(value(p==1, 1)));
    Aircraft.Certification.Aerodynamic_data.airfoil_second_panel.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('airfoil_third_panel', label);
    Aircraft.Certification.Aerodynamic_data.airfoil_third_panel.value = char(table2array(value(p==1, 1)));
    Aircraft.Certification.Aerodynamic_data.airfoil_third_panel.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('airfoil_fourth_panel', label);
    Aircraft.Certification.Aerodynamic_data.airfoil_fourth_panel.value = char(table2array(value(p==1, 1)));
    Aircraft.Certification.Aerodynamic_data.airfoil_fourth_panel.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('sweep_first', label);
    Aircraft.Geometry.Wing.sweep_first.value = str2double(table2array(value(p==1, 1)));  
    Aircraft.Geometry.Wing.sweep_first.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('sweep_second', label);
    Aircraft.Geometry.Wing.sweep_second.value = str2double(table2array(value(p==1, 1)));  
    Aircraft.Geometry.Wing.sweep_second.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('sweep_third', label);
    Aircraft.Geometry.Wing.sweep_third.value = str2double(table2array(value(p==1, 1)));  
    Aircraft.Geometry.Wing.sweep_third.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('dihedral_first', label);
    Aircraft.Geometry.Wing.dihedral_first.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Wing.dihedral_first.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('dihedral_second', label);
    Aircraft.Geometry.Wing.dihedral_second.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Wing.dihedral_second.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('dihedral_third', label);
    Aircraft.Geometry.Wing.dihedral_third.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Wing.dihedral_third.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('twist_angle_second', label);
    Aircraft.Geometry.Wing.twist_angle_second.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Wing.twist_angle_second.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('twist_angle_third', label);
    Aircraft.Geometry.Wing.twist_angle_third.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Wing.twist_angle_third.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('twist_angle_fourth', label);
    Aircraft.Geometry.Wing.twist_angle_fourth.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Wing.twist_angle_fourth.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('vertical_sweep', label);
    Aircraft.Geometry.Vertical.sweep.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Vertical.sweep.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('sweep_loc_vert', label);
    Aircraft.Geometry.Vertical.sweeploc.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Vertical.sweeploc.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('sec_sweeploc_vert', label);
    Aircraft.Geometry.Vertical.secsweeploc.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Vertical.secsweeploc.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('vertical_dihedral', label);
    Aircraft.Geometry.Vertical.dihedral.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Vertical.dihedral.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('vertical_twist', label);
    Aircraft.Geometry.Vertical.twist.value       = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Vertical.twist.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('ver_twist_loc', label);
    Aircraft.Geometry.Vertical.twistloc.value       = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Vertical.twistloc.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('vertical_x_leadingedge', label);
    Aircraft.Geometry.Vertical.xle.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Vertical.xle.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('vertical_croot', label);
    Aircraft.Geometry.Vertical.croot.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Vertical.croot.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('vertical_ctip', label);
    Aircraft.Geometry.Vertical.ctip.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Vertical.ctip.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('vertical_xtip_leadingedge', label);   
    Aircraft.Geometry.Vertical.xtip_le.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Vertical.xtip_le.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('vertical_span', label);  
    Aircraft.Geometry.Vertical.b.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Vertical.b.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('vertical_z_position', label);
    Aircraft.Geometry.Vertical.zpos.value = str2double(table2array(value(p==1, 1)));   %1.0; 
    Aircraft.Geometry.Vertical.zpos.Attributes.unit = char(table2array(unit(p==1,1))); % "% of fuselage diameter";
    p = strcmp('flaps_eta_iniziale', label);
    Aircraft.Geometry.Flaps.eta_inner.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Flaps.eta_inner.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('flaps_eta_finale', label);    
    Aircraft.Geometry.Flaps.eta_outer.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Flaps.eta_outer.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('flaps_y_iniziale', label);    
    Aircraft.Geometry.Flaps.y_inner.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Flaps.y_inner.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('flaps_y_finale', label);
    Aircraft.Geometry.Flaps.y_outer.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Flaps.y_outer.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('flaps_span_b', label);
    Aircraft.Geometry.Flaps.b.value = str2double(table2array(value(p==1, 1))); 
    Aircraft.Geometry.Flaps.b.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('cf_c_root_flap', label); 
    Aircraft.Geometry.Flaps.cf_c_root.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Flaps.cf_c_root.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('cf_c_tip_flap', label); 
    Aircraft.Geometry.Flaps.cf_c_tip.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Flaps.cf_c_tip.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('croot_flap', label);
    Aircraft.Geometry.Flaps.croot.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Flaps.croot.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('ctip_flap', label);
    Aircraft.Geometry.Flaps.ctip.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Flaps.ctip.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('cf_flap', label);
    Aircraft.Geometry.Flaps.cf.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Flaps.cf.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('flap_surface_S', label);
    Aircraft.Geometry.Flaps.S.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Geometry.Flaps.S.Attributes.unit = char(table2array(unit(p==1,1))); 
    % -------------------------------------------------------------------------- 
    p = strcmp('vertical_yaw_angle', label);
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('verticaloads_dr', label);
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.dr.value = char(table2array(value(p==1, 1)));
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.dr.value = str2num(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.dr.value);
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.dr.Attributes.unit = char(table2array(unit(p==1,1)));    
    p = strcmp('vertical_CY_vector', label);
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_vector.value = char(table2array(value(p==1, 1)));  % [0.00003615686 0.01291153]';
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_vector.value = str2num(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_vector.value);
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_vector.Attributes.unit = char(table2array(unit(p==1,1))); 
    p = strcmp('vertical_CY0', label);
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_0.value = str2double(table2array(value(p==1, 1)));
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_0.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('vertical_CYdr', label);
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CYdr.value = str2double(table2array(value(p==1, 1)));   % 0.000644;
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CYdr.Attributes.unit = char(table2array(unit(p==1,1))); % "1/deg";
    p = strcmp('maximum_delta_rudder', label);
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.value = str2double(table2array(value(p==1, 1)));   % 30.0;
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.Attributes.unit = char(table2array(unit(p==1,1))); % "deg";
    % CASE 2
    p = strcmp('CY_vertical_tailplane', label);
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_VTP.value = str2double(table2array(value(p==1, 1)));   % 0.0245;
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_VTP.Attributes.unit = char(table2array(unit(p==1,1))); % "Non dimensional";
    % CASE 3
    p = strcmp('case3_CY_vert_tp', label);
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_VTP.value = str2double(table2array(value(p==1, 1)));   % 0.0233;
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_VTP.Attributes.unit = char(table2array(unit(p==1,1))); % "Non dimensional";
    p = strcmp('vert_tp_K', label);
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.value = str2double(table2array(value(p==1, 1)));   % 0.3; 
    Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.Attributes.unit = char(table2array(unit(p==1,1))); % "m"; 

end

Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.Attributes.unit = "Positive g";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.Attributes.unit = "Negative g";
Aircraft.Geometry.Wing.twist_angle.value = Aircraft.Geometry.Wing.twist_angle_first.value;
Aircraft.Geometry.Wing.twist_angle.Attributes.unit = "deg";

% 
%          'Aircraft.Geometry.Wing.mac', ...
%          'Aircraft.Certification.Aerodynamic_data.xac.value', ...
%          'Aircraft.Certification.Aerodynamic_data.yac.value', ...
%          'Aircraft.Certification.Aerodynamic_data.zac.value', ...
%% Derived values (also inputs)
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value = calcn(Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.Attributes.unit = "Positive g";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.value = calcn(Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.Attributes.unit = "Negative g";
% -----------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------

% function Aircraft = FromFileCertification_fun(Aircraft,varargin)
% %% READ INPUT FILE
% %% READING EXCEL
% % cd ../_Utilities
% fprintf('\n');
% fprintf('Input File: ');
% if nargin==1
%     filename='Aircraft.xlsx';
% else
%     filename=[varargin{1},'.xlsx'];
% end
% fprintf(filename);
% 
% %reading excel spreadsheet with readtable function with options
% opts = spreadsheetImportOptions("NumVariables", 3);
% % Specify sheet and range
% %opts.Sheet = "Sheet1";
% %opts.DataRange = "A2:C25";
% 
% % Specify column names and types
% opts.VariableNames = ["field", "value", "unit"];
% opts.VariableTypes = ["char", "char", "char"];
% 
% % Specify variable properties
% %opts = setvaropts(opts, ["field", "value", "unit"], "WhitespaceRule", "preserve");
% %opts = setvaropts(opts, ["field", "value", "unit"], "EmptyFieldRule", "auto");
% Table = readtable(filename, opts, "UseExcel", false);
% 
% value=Table(:,2);
% unit=Table(:,3);
% label=Table{:,1};
% fprintf('\n');
% % cd ../_PreDesign
% 
% VarText={'Aircraft_name', ...      % PRELIMINARY INPUTS
%          'Aircraft_type', ...
%          'Designer', ...
%          'Regulation', ...
%          'Date', ...
%          'Amendment', ...
%          'Airloads_flag', ...
%          'vertical_tail_flag', ... % PRELIMINARY INPUTS
%          'MTOM', ...               % AIRCRAFT MASS/WEIGHT
%          'EM', ...
%          'UM', ...
%          'FM', ...
%          'CM', ...
%          'Oil_mass', ...           % AIRCRAFT MASS/WEIGHT
%          'Wing_type', ...          % WING CHARACTERISTICS
%          'Wing_surface', ...
%          'Wing_span', ...
%          'Wing_dihedral', ...
%          'Wing_loading', ...
%          'Wing_xac', ...
%          'panel_span1', ...
%          'panel_span2', ...
%          'panel_span3', ...
%          'first_airfoil_section', ...
%          'second_airfoil_section', ...
%          'third_airfoil_section', ...
%          'fourth_airfoil_section', ...
%          'aspect_ratio', ...
%          'first_sweep_value', ...
%          'second_sweep_value', ...
%          'third_sweep_value', ...
%          'first_dihedral_value', ...
%          'second_dihedral_value', ...
%          'third_dihedral_value', ...
%          'first_twist_value', ...
%          'second_twist_value', ...
%          'third_twist_value', ...
%          'fourth_twist_value', ...
%          'Wing_root_chord', ...
%          'Wing_tip_chord', ...
%          'Wing_kink_one', ...
%          'Wing_kink_two', ...
%          'wing_camber', ...
%          'wing_camber_location', ...   
%          'thickness_ratio', ...
%          'secondary_sweep_location', ...
%          'sweep_location', ...
%          'wing_xle', ...
%          'wing_yle', ...
%          'wing_zle', ...               % WING CHARACTERISTICS
%          'Correction_factor_flag1', ... % ENGINE DATA
%          'Correction_factor_flag2', ...
%          'takeoff_power', ...
%          'takeoff_rpm', ...
%          'reduction_ratio', ...
%          'max_cont_power', ...
%          'max_cont_rpm', ...
%          'engine_limit_side_load', ...
%          'engine_mount_mass', ...
%          'engine_accessories', ...
%          'engine_pitch_speed', ...
%          'engine_yaw_speed', ...
%          'engine_normal_load_factor', ...
%          'thrust_axes', ...
%          'single_tank_capacity', ...
%          'fuel_capacity', ...
%          'oil_capacity', ...           
%          'engine_block_mass', ...      % ENGINE DATA 
%          'propeller_diameter', ...     % PROPELLER DATA
%          'propeller_numberof_blades', ...
%          'propeller_mass', ...
%          'prop_torque_at_takeoff_pow', ...
%          'propeller_polar_moment', ...
%          'propeller_spinner_mass', ... % PROPELLER DATA
%          'Max_lift_coeff_clean', ...   % AERODYNAMIC DATA
%          'Max_lift_coeff_takeoff', ...
%          'Max_lift_coeff_landing', ...
%          'Min_lift_coeff', ...
%          'Normal_force_coeff_slope_deg', ...
%          'Normal_force_coeff_slope', ...
%          'CD0', ...
%          'CD_landing_gear', ...
%          'oswald_efficiency', ...
%          'zero_lift_coefficient', ...
%          'clstar_coefficient', ...
%          'CM0', ...
%          'Cmlanding_gear', ...
%          'CMCL_slope', ...
%          'aileron_chdelta', ...
%          'aileron_chalfa', ...
%          'elevator_chdelta', ...
%          'elevator_chalfa', ...
%          'rudder_chdelta', ...
%          'rudder_chalfa', ...
%          'rudder_max_deflection', ... 
%          'elevator_max_deflection', ...
%          'ailerons_max_deflection', ...
%          'flaps_max_deflection', ...
%          'wing_body_aoa', ...
%          'wing_body_cl', ...
%          'wing_body_cm', ...
%          'wing_body_cd', ...           % AERODYNAMIC DATA
%          'sea_level', ...              % STANDARD ATMOSPHERE
%          'operative_ceiling', ...
%          'theoretical_ceiling', ...    % STANDARD ATMOSPHERE
%          'max_forward_cg', ...         % CENTRE OF GRAVITY
%          'max_aft_cg', ...
%          'xcg', ...
%          'ycg', ...
%          'zcg', ...                    % CENTRE OF GRAVITY
%          'wing_xac', ...               % AERODYNAMIC CENTRE
%          'ailerons_surface', ...       % AILERON CHARACTERISTICS
%          'ailerons_span', ...
%          'ailerons_y_inner', ...
%          'ailerons_y_outer', ...
%          'ailerons_eta_inner', ...
%          'ailerons_eta_outer', ...
%          'ailerons_ca', ...
%          'ailerons_cb', ...
%          'ca_c_root_ailerons', ...
%          'ca_c_tip_ailerons', ...
%          'ca_root_ailerons', ...
%          'ca_tip_ailerons', ...
%          'ailerons_moment_arm', ...   
%          'cf_ailerons', ...            % AILERON CHARACTERISTICS
%          'elevator_surface', ...       % ELEVATOR CHARACTERISTICS
%          'elevator_span', ...
%          'elevator_y_inner', ...
%          'elevator_y_outer', ...
%          'elevator_eta_inner', ...
%          'elevator_eta_outer', ...
%          'elevator_ca', ...
%          'elevator_cb', ...
%          'elevator_type', ...
%          'ce_c_root_elevator', ...
%          'ce_c_tip_elevator', ...
%          'ce_root_elevator', ...
%          'ce_tip_elevator', ...
%          'elevator_S_hinge', ...
%          'elevator_moment_arm', ...
%          'cf_elevator', ...            % ELEVATOR CHARACTERISTICS
%          'rudder_surface', ...         % RUDDER CHARACTERISTICS
%          'rudder_span', ...
%          'rudder_y_inner', ...
%          'rudder_y_outer', ...
%          'rudder_eta_inner', ...
%          'rudder_eta_outer', ...
%          'rudder_ca', ...
%          'rudder_cb', ...
%          'cr_c_root_rudder', ...
%          'cr_c_tip_rudder', ...
%          'cr_root_rudder', ...
%          'cr_tip_rudder', ...
%          'rudder_moment_arm', ...
%          'cf_rudder', ...              % RUDDER CHARACTERISTICS
%          'flaps_surface', ...          % FLAPS CHARACTERISTICS
%          'flaps_span', ...
%          'flaps_y_inner', ...
%          'flaps_y_outer', ...
%          'flaps_eta_inner', ...
%          'flaps_eta_outer', ...
%          'flaps_ca', ...
%          'flaps_cb', ...
%          'cf_c_root_flaps', ...
%          'cf_c_tip_flaps', ...
%          'cf_root_flaps', ...
%          'cf_tip_flaps', ...
%          'flaps_moment_arm', ...
%          'cf_flaps', ...
%          'flaps_type', ...             % FLAPS CHARACTERISTICS
%          'overall_length', ...         % FUSELAGE CHARACTERISTICS
%          'overall_width', ...
%          'overall_height', ...         
%          'id', ...                     
%          'type', ...                   % FUSELAGE CHARACTERISTICS
%          'ht_surface', ...             % HORIZONTAL TAILPLANE 
%          'ht_span', ...
%          'ht_croot', ...
%          'ht_ctip', ...
%          'ht_arm', ...
%          'ht_camber', ...
%          'ht_camber_location', ...
%          'ht_thickness_ratio', ...
%          'ht_twist_angle', ...
%          'ht_twistangle_location', ...
%          'ht_type', ...
%          'ht_x_location', ...
%          'ht_y_location', ...
%          'ht_z_location', ...
%          'ht_sweep', ...
%          'ht_sweep_location', ...
%          'ht_sec_sweep_location', ...
%          'ht_dihedral_angle', ...      % HORIZONTAL TAILPLANE 
%          'vt_surface', ...             % VERTICAL TAILPLANE
%          'vt_span', ...
%          'vt_croot', ...
%          'vt_ctip', ...
%          'vt_arm', ...
%          'vt_twist_angle', ...
%          'vt_twistangle_location', ...
%          'vt_dihedral_angle', ...
%          'vt_x_location', ...
%          'vt_y_location', ...
%          'vt_z_location', ...
%          'vt_x_location_perc', ...
%          'vt_y_location_perc', ...
%          'vt_z_location_perc', ...
%          'vt_yaw_angle', ...
%          'vt_sweep_angle', ...
%          'vt_sweeplocation', ...
%          'vt_sec_sweeplocation', ...   % VERTICAL TAILPLANE
%          'elevator_overhang', ...      % HINGE DATA
%          'rudder_overhang', ...        % HINGE DATA
%          'vertical_tail_delta_rudder', ... % LATERAL FORCE DATA
%          'CY_case1', ...                   
%          'CY0', ...
%          'CY_dr', ...
%          'CY_case2', ...
%          'CY_case3', ...
%          'vertical_tailplane_K', ...       % LATERAL FORCE DATA 
%          'Min_design_cruise_speed', ...    % DESIGN AIRSPEED 
%          };
%      VarText = VarText';
% % -----------------------------
% % label = label(2:end,1);
% % unit  = unit(2:end,1);
% % value = value(2:end,1);
% Debugvartext = VarText';
% % VarText = VarText';
% % -----------------------------
% Index=zeros(1,length(VarText));
% for i=1:length(VarText)
%     x = find(~cellfun(@isempty, strfind(label,VarText{i})));
%     if length(x) > 1
%         x = x(1);
%         Index(i) = x;
%     elseif length(x) == 1
%         Index(i) = x;
%     end
% %     Index(i) = find(~cellfun(@isempty, strfind(label,VarText{i})));
% %     Index(i) = find(~cellfun(@isempty, strfind(label,Debugvartext{i})));
% end
% LengthIndex = nnz(Index);
% if LengthIndex==length(VarText)
%     % ---------------------------------------------------------------------
%     p = strcmp('Aircraft_name',label);
%     Aircraft.Certification.Aircraft_Name.value = char(table2array(value(p==1,1)));
%     p = strcmp('Aircraft_type',label);
%     Aircraft.Certification.Aircraft_Name.Attributes.type = char(table2array(value(p==1,1)));
%     p = strcmp('Designer',label);
%     Aircraft.Certification.Aircraft_Name.Attributes.designer = char(table2array(value(p==1,1)));
%     p = strcmp('Regulation',label);
%     Aircraft.Certification.Regulation.value = char(table2array(value(p==1,1)));
%     p = strcmp('Date',label);
%     Aircraft.Certification.Regulation.Attributes.Date = char(table2array(value(p==1,1)));
%     p = strcmp('Amendment',label);
%     Aircraft.Certification.Regulation.Attributes.Amendment = char(table2array(value(p==1,1)));
%     p = strcmp('Airloads_flag', label);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Airload_case.Attributes.case = char(table2array(value(p==1,1)));
%     p = strcmp('vertical_tail_flag', label);
%     Aircraft.Geometry.Vertical.empennage_flag.value = char(table2array(value(p==1,1)));
%     % ---------------------------------------------------------------------
%     p = strcmp('MTOM', label);
%     Aircraft.Weight.I_Level.W_maxTakeOff.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Weight.I_Level.W_maxTakeOff.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('EM', label);
%     Aircraft.Weight.I_Level.W_EmptyMass.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Weight.I_Level.W_EmptyMass.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('UM', label);
%     Aircraft.Weight.I_Level.W_UsefulMass.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Weight.I_Level.W_UsefulMass.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('FM', label);
%     Aircraft.Weight.I_Level.W_FuelMass.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Weight.I_Level.W_FuelMass.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('CM', label);
%     Aircraft.Weight.I_Level.W_CrewMass.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Weight.I_Level.W_CrewMass.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Oil_mass', label);
%     Aircraft.Weight.I_Level.Oil_mass.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Weight.I_Level.Oil_mass.Attributes.unit = char(table2array(unit(p==1,1)));
%     % ---------------------------------------------------------------------
%     p = strcmp('Wing_type', label);
%     Aircraft.Geometry.Wing.type.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Wing.type.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Wing_surface', label);
%     Aircraft.Geometry.Wing.S.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Wing.S.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Wing_span', label);
%     Aircraft.Geometry.Wing.b.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Wing.b.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Wing_dihedral', label);
%     Aircraft.Geometry.Wing.dihedral.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Wing.dihedral.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Wing_loading', label);
%     Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Certification.Performance.I_Level.Wing_loading_SI.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Wing_xac', label);
%     Aircraft.Geometry.General.XAC_nondim.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.General.XAC_nondim.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('panel_span1', label);
%     Aircraft.Geometry.Wing.panel_span1.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Wing.panel_span1.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('panel_span2', label);
%     Aircraft.Geometry.Wing.panel_span2.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Wing.panel_span2.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('panel_span3', label);
%     Aircraft.Geometry.Wing.panel_span3.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Wing.panel_span3.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('first_airfoil_section', label);
%     Aircraft.Certification.Aerodynamic_data.airfoil_first_panel.value = char(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.airfoil_first_panel.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('second_airfoil_section', label);
%     Aircraft.Certification.Aerodynamic_data.airfoil_second_panel.value = char(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.airfoil_second_panel.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('third_airfoil_section', label);
%     Aircraft.Certification.Aerodynamic_data.airfoil_third_panel.value = char(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.airfoil_third_panel.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('fourth_airfoil_section', label);
%     Aircraft.Certification.Aerodynamic_data.airfoil_fourth_panel.value = char(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.airfoil_fourth_panel.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('aspect_ratio', label);
%     Aircraft.Geometry.Wing.AR.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Wing.AR.Attributes.unit = char(table2array(unit(p==1, 1))); 
%     p = strcmp('first_sweep_value', label);
%     Aircraft.Geometry.Wing.sweep_first.value = str2double(table2array(value(p==1, 1)));  
%     Aircraft.Geometry.Wing.sweep_first.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('second_sweep_value', label);
%     Aircraft.Geometry.Wing.sweep_second.value = str2double(table2array(value(p==1, 1)));  
%     Aircraft.Geometry.Wing.sweep_second.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('third_sweep_value', label);
%     Aircraft.Geometry.Wing.sweep_third.value = str2double(table2array(value(p==1, 1)));  
%     Aircraft.Geometry.Wing.sweep_third.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('first_dihedral_value', label);
%     Aircraft.Geometry.Wing.dihedral_first.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.dihedral_first.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('second_dihedral_value', label);
%     Aircraft.Geometry.Wing.dihedral_second.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.dihedral_second.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('third_dihedral_value', label);
%     Aircraft.Geometry.Wing.dihedral_third.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.dihedral_third.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('first_twist_value', label);
%     Aircraft.Geometry.Wing.twist_angle_first.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Wing.twist_angle_first.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('second_twist_value', label);
%     Aircraft.Geometry.Wing.twist_angle_second.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.twist_angle_second.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('third_twist_value', label);
%     Aircraft.Geometry.Wing.twist_angle_third.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.twist_angle_third.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('fourth_twist_value', label);
%     Aircraft.Geometry.Wing.twist_angle_fourth.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.twist_angle_fourth.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('twist_location', label);
%     Aircraft.Geometry.Wing.twist_location.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Wing.twist_location.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Wing_root_chord', label);
%     Aircraft.Geometry.Wing.croot.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Wing.croot.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Wing_tip_chord', label);
%     Aircraft.Geometry.Wing.ctip.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Wing.ctip.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Wing_kink_one', label);
%     Aircraft.Geometry.Wing.chord_kink_one.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Wing.chord_kink_one.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Wing_kink_two', label);
%     Aircraft.Geometry.Wing.chord_kink_two.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Wing.chord_kink_two.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('wing_camber', label);
%     Aircraft.Geometry.Wing.camber.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Wing.camber.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('wing_camber_location', label);
%     Aircraft.Geometry.Wing.camberloc.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.camberloc.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('thickness_ratio', label);
%     Aircraft.Geometry.Wing.thickchord.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.thickchord.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('secondary_sweep_location', label);
%     Aircraft.Geometry.Wing.secondary_sweep_location.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.secondary_sweep_location.Attributes.unit = char(table2array(unit(p==1,1)));  
%     p = strcmp('sweep_location', label);
%     Aircraft.Geometry.Wing.sweep_location.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.sweep_location.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('wing_xle', label); 
%     Aircraft.Geometry.Wing.xle.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.xle.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('wing_yle', label); 
%     Aircraft.Geometry.Wing.yle.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.yle.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('wing_zle', label); 
%     Aircraft.Geometry.Wing.zle.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.zle.Attributes.unit = char(table2array(unit(p==1,1)));
%     % ---------------------------------------------------------------------
%     p = strcmp('Correction_factor_flag1', label);
%     Aircraft.Engine.Correction_factor.Attributes.flag1 = char(table2array(value(p==1, 1)));
%     p = strcmp('Correction_factor_flag2', label);    
%     Aircraft.Engine.Correction_factor.Attributes.flag2 = char(table2array(value(p==1, 1)));
%     p = strcmp('takeoff_power', label); 
%     Aircraft.Engine.Takeoff.Power.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Takeoff.Power.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('takeoff_rpm', label);
%     Aircraft.Engine.Takeoff.RPM.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Takeoff.RPM.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('reduction_power', label); 
%     Aircraft.Engine.Reduction_ratio.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Reduction_ratio.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('max_cont_power', label); 
%     Aircraft.Engine.Max_Continous.Power.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Max_Continous.Power.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('max_cont_rpm', label); 
%     Aircraft.Engine.Max_Continous.RPM.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Max_Continous.RPM.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('engine_limit_side_load', label);
%     Aircraft.Engine.Limit_side_load.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Limit_side_load.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('engine_mount_mass', label);
%     Aircraft.Engine.Engine_mount_mass.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Engine_mount_mass.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('engine_accessories', label);
%     Aircraft.Engine.Engine_accessories_mass.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Engine_accessories_mass.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('engine_pitch_speed', label);
%     Aircraft.Engine.Pitch_speed.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Pitch_speed.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('engine_yaw_speed', label);
%     Aircraft.Engine.Yaw_speed.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Yaw_speed.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('engine_normal_load_factor', label);
%     Aircraft.Engine.Engine_normal_load_factor.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Engine_normal_load_factor.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('thrust_axes', label);
%     Aircraft.Geometry.Engine.Primary.Thrust_axes.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Engine.Primary.Thrust_axes.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('single_tank_capacity', label);
%     Aircraft.Engine.single_tank_capacity.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.single_tank_capacity.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('fuel_capacity', label);
%     Aircraft.Engine.fuel_capacity.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.fuel_capacity.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('oil_capacity', label);
%     Aircraft.Engine.oil_capacity.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.oil_capacity.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('engine_block_mass', label);
%     Aircraft.Engine.engine_block_mass.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.engine_block_mass.Attributes.unit = char(table2array(unit(p==1,1)));
%     % ---------------------------------------------------------------------
%     p = strcmp('propeller_diameter', label);
%     Aircraft.Engine.Propeller_polar_moment.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Engine.Propeller_polar_moment.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('propeller_numberof_blades', label);
%     Aircraft.Engine.Propeller_blade_number.value = str2double(table2array(value(p==1, 1)));
%     p = strcmp('propeller_mass', label);
%     Aircraft.Engine.Propeller_mass.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Engine.Propeller_mass.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('prop_torque_at_takeoff_pow', label);
%     Aircraft.Engine.Takeoff.Propeller_torque.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Engine.Takeoff.Propeller_torque.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('propeller_polar_moment', label);
%     Aircraft.Engine.Propeller_polar_moment.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Engine.Propeller_polar_moment.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('propeller_spinner_mass', label);
%     Aircraft.Engine.Propeller_spinner_mass.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Propeller_spinner_mass.Attributes.unit = char(table2array(unit(p==1,1)));
%     % ---------------------------------------------------------------------            
%     p = strcmp('Max_lift_coeff_clean',label);
%     Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value = str2double(table2array(value(p==1,1))); 
%     Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.Attributes.unit = char(table2array(unit(p==1,1)));          
%     p = strcmp('Max_lift_coeff_takeoff', label);
%     Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_takeoff.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_takeoff.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Max_lift_coeff_landing', label);
%     Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_landing.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_landing.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Min_lift_coeff',label);
%     Aircraft.Certification.Aerodynamic_data.Min_Lift_Coefficient.value = str2double(table2array(value(p==1,1))); 
%     Aircraft.Certification.Aerodynamic_data.Min_Lift_Coefficient.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('Normal_force_coeff_slope_deg',label);
%     Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('Normal_force_coeff_slope',label);
%     Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value = str2double(table2array(value(p==1,1))); 
%     Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('bcg', label);
%     Aircraft.Geometry.General.bcg.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.General.bcg.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('CD0', label);
%     Aircraft.Certification.Aerodynamic_data.CD0.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Certification.Aerodynamic_data.CD0.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('CD_landing_gear', label);
%     Aircraft.Certification.Aerodynamic_data.CD_landing_gear.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Certification.Aerodynamic_data.CD_landing_gear.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('oswald_efficiency', label);
%     Aircraft.Certification.Aerodynamic_data.e.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Certification.Aerodynamic_data.e.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('zero_lift_coefficient', label);
%     Aircraft.Certification.Aerodynamic_data.CL0.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Certification.Aerodynamic_data.CL0.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('clstar_coefficient', label);
%     Aircraft.Certification.Aerodynamic_data.CL_star.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Certification.Aerodynamic_data.CL_star.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('CM0', label);
%     Aircraft.Certification.Aerodynamic_data.CM0.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Certification.Aerodynamic_data.CM0.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Cmlanding_gear', label);
%     Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Certification.Aerodynamic_data.CM_landing_gear.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('CMCL_slope', label);
%     Aircraft.Certification.Aerodynamic_data.CMCL.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Certification.Aerodynamic_data.CMCL.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('aileron_chdelta', label);
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_delta_rad.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_delta_rad.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('aileron_chalfa', label);
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_alfa_rad.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_alfa_rad.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('elevator_chdelta', label);
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_delta_rad.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_delta_rad.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('elevator_chalfa', label);
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_alfa_rad.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_alfa_rad.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('rudder_chdelta', label);
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_delta_rad.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_delta_rad.Attributes.unit = char(table2array(unit(p==1,1)));  
%     p = strcmp('rudder_chalfa', label);
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_alfa_rad.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_alfa_rad.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('rudder_max_deflection', label);
%     Aircraft.Geometry.Rudder.max_deflection.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Rudder.max_deflection.Attributes.unit = char(table2array(unit(p==1, 1))); 
%     p = strcmp('elevator_max_deflection', label);
%     Aircraft.Geometry.Elevator.max_deflection.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Elevator.max_deflection.Attributes.unit = char(table2array(unit(p==1, 1)));
%     p = strcmp('ailerons_max_deflection', label);
%     Aircraft.Geometry.Aileron.max_deflection.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Aileron.max_deflection.Attributes.unit = char(table2array(unit(p==1, 1))); 
%     p = strcmp('wing_body_aoa',label);
%     Aircraft.Certification.Aerodynamic_data.alpha.value = char(table2array(value(p==1,1)));
%     Aircraft.Certification.Aerodynamic_data.alpha.Attributes.unit = char(table2array(unit(p==1, 1)));
%     p = strcmp('wing_body_cl',label);
%     Aircraft.Certification.Aerodynamic_data.CL.value = char(table2array(value(p==1,1)));
%     Aircraft.Certification.Aerodynamic_data.CL.Attributes.unit = char(table2array(unit(p==1, 1)));
%     p = strcmp('wing_body_cm',label);
%     Aircraft.Certification.Aerodynamic_data.CM.value = char(table2array(value(p==1,1)));
%     Aircraft.Certification.Aerodynamic_data.CM.Attributes.unit = char(table2array(unit(p==1, 1)));
%     p = strcmp('wing_body_cd',label);
%     Aircraft.Certification.Aerodynamic_data.CD.value = char(table2array(value(p==1,1))); 
%     Aircraft.Certification.Aerodynamic_data.CD.Attributes.unit = char(table2array(unit(p==1, 1)));
%     % ---------------------------------------------------------------------
%     p = strcmp('sea_level', label);
%     Aircraft.Certification.ISA_Condition.Sea_level.Altitude.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.ISA_Condition.Sea_level.Altitude.Attributes.unit = char(table2array(unit(p==1,1)));  
%     p = strcmp('operative_ceiling', label);
%     Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('theoretical_ceiling', label);
%     Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.Attributes.unit = char(table2array(unit(p==1,1)));
%     % ---------------------------------------------------------------------
%     p = strcmp('max_forward_cg', label);
%     Aircraft.Geometry.General.Max_forward_cg.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.General.Max_forward_cg.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('max_aft_cg', label);
%     Aircraft.Geometry.General.Max_aft_cg.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.General.Max_aft_cg.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('xcg_nondimensional', label);
%     Aircraft.Geometry.General.XCG_nondim.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.General.XCG_nondim.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('xcg', label);
%     Aircraft.Geometry.General.xcg.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.General.xcg.Attributes.unit = char(table2array(unit(p==1,1))); % Measured from the aircraft nose
%     p = strcmp('ycg', label);
%     Aircraft.Geometry.General.ycg.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.General.ycg.Attributes.unit = char(table2array(unit(p==1,1))); % Measured from the aircraft nose
%     p = strcmp('zcg', label);
%     Aircraft.Geometry.General.zcg.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.General.zcg.Attributes.unit = char(table2array(unit(p==1,1))); % Measured from the aircraft nose
%     p = strcmp('wing_xac', label);
%     Aircraft.Geometry.General.XAC_nondim.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.General.XAC_nondim.Attributes.unit = char(table2array(unit(p==1,1)));
%     % ---------------------------------------------------------------------
%     
%     p = strcmp('rudder_surface', label);
%     Aircraft.Geometry.Rudder.S.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Rudder.S.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('chord_rudder', label);
%     Aircraft.Geometry.Rudder.chord.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Rudder.chord.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('chord_ratio_rudder_cf_c', label);
%     Aircraft.Geometry.Rudder.chord_ratio_cf_c.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Rudder.chord_ratio_cf_c.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('overhang_rudder', label);
%     Aircraft.Geometry.Rudder.overhang.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Rudder.overhang.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('span_ratio_rudder', label);
%     Aircraft.Geometry.Rudder.span_ratio.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Rudder.span_ratio.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('max_deflection_rudder', label);
%     Aircraft.Geometry.Rudder.max_deflection.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Rudder.max_deflection.Attributes.unit = char(table2array(unit(p==1,1)));
%     
%     % ---------------------------------------------------------------------
%     
%     p = strcmp('ailerons_surface', label); 
%     Aircraft.Geometry.Aileron.S.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Aileron.S.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('ailerons_span', label);
%     Aircraft.Geometry.Aileron.b.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Aileron.b.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('ailerons_y_inner', label);
%     Aircraft.Geometry.Aileron.y_inner.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Aileron.y_inner.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('ailerons_y_outer', label);
%     Aircraft.Geometry.Aileron.y_outer.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Aileron.y_outer.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('ailerons_eta_inner', label); 
%     Aircraft.Geometry.Aileron.eta_inner.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Aileron.eta_inner.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('ailerons_eta_outer', label);
%     Aircraft.Geometry.Aileron.eta_outer.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Aileron.eta_outer.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('ailerons_ca', label);
%     Aircraft.Geometry.Aileron.ca.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Aileron.ca.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('ailerons_cb', label);
%     Aircraft.Geometry.Aileron.cb.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Aileron.ca.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('ca_c_root_ailerons', label);
%     Aircraft.Geometry.Aileron.ca_c_root.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Aileron.ca_c_root.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('ca_c_tip_ailerons', label);
%     Aircraft.Geometry.Aileron.ca_c_tip.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Aileron.ca_c_tip.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('ca_root_ailerons', label);
%     Aircraft.Geometry.Aileron.ca_root.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Aileron.ca_root.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('ca_tip_ailerons', label);
%     Aircraft.Geometry.Aileron.ca_tip.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Aileron.ca_tip.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('ailerons_moment_arm', label);
%     Aircraft.Geometry.Aileron.moment_arm.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Aileron.moment_arm.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('cf_ailerons', label);
%     Aircraft.Geometry.Aileron.cf.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Aileron.cf.Attributes.unit = char(table2array(unit(p==1,1))); 
%     % ---------------------------------------------------------------------
%     p = strcmp('elevator_surface', label); 
%     Aircraft.Geometry.Elevator.S.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.S.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('elevator_span', label); 
%     Aircraft.Geometry.Elevator.b.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.b.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('elevator_y_inner', label);
%     Aircraft.Geometry.Elevator.y_inner.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.y_inner.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('elevator_y_outer', label);
%     Aircraft.Geometry.Elevator.y_outer.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.y_outer.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('elevator_eta_inner', label);
%     Aircraft.Geometry.Elevator.eta_inner.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.eta_inner.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('elevator_eta_outer', label);
%     Aircraft.Geometry.Elevator.eta_outer.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.eta_outer.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('elevator_ca', label);
%     Aircraft.Geometry.Elevator.ca.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.ca.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('elevator_cb', label);
%     Aircraft.Geometry.Elevator.cb.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.cb.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('elevator_type', label);
%     Aircraft.Certification.Aerodynamic_data.Horizontal.tau.Attributes.flag = char(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.Horizontal.tau.Attributes.type = char(table2array(unit(p==1,1)));
%     p = strcmp('ce_c_root_elevator', label);
%     Aircraft.Geometry.Elevator.cf_c_inner.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.cf_c_inner.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('ce_c_tip_elevator', label);
%     Aircraft.Geometry.Elevator.cf_c_outer.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.cf_c_outer.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('ce_root_elevator', label);
%     Aircraft.Geometry.Elevator.cf_inner.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.cf_inner.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('ce_tip_elevator', label);
%     Aircraft.Geometry.Elevator.cf_outer.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.cf_outer.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('elevator_S_hinge', label);
%     Aircraft.Geometry.Elevator.S_hinge.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.S_hinge.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('elevator_moment_arm', label);
%     Aircraft.Geometry.Horizontal.l.value = str2double(table2array(value(p==1, 1)));     % tail arm in meters
%     Aircraft.Geometry.Horizontal.l.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('cf_elevator', label);
%     Aircraft.Geometry.Elevator.cf.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Elevator.cf.Attributes.unit = char(table2array(unit(p==1,1))); 
%     
%     % ---------------------------------------------------------------------
%     p = strcmp('Aircraft.Certification.Regulation.SubpartC.Flightloads.Max_Continous_Power_Speed_VH',label);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Max_Continous_Power_Speed_VH.value = str2double(table2array(value(p==1,1))); 
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Max_Continous_Power_Speed_VH.Attributes.unit = char(table2array(unit(p==1,1)));      
%     p = strcmp('Aircraft.Certification.Regulation.SubpartC.Flightloads.Min_Design_Cruise_Speed',label);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Min_Design_Cruise_Speed.value = str2double(table2array(value(p==1,1))); 
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Min_Design_Cruise_Speed.Attributes.unit = char(table2array(unit(p==1,1)));
% %     p = strcmp('Aircraft.Geometry.Wing.mac',label);
% %     Aircraft.Geometry.Wing.mac.value = str2double(table2array(value(p==1,1)));
% %     Aircraft.Geometry.Wing.mac.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('Aircraft.Certification.Regulation.SubpartC.Gustloads.Gust_speed_cruise',label);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Certification.Regulation.SubpartC.Gustloads.Gust_speed_dive',label);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.Attributes.unit = char(table2array(unit(p==1,1)));   
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope',label);
%     Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.alpha.value',label);
%     Aircraft.Certification.Aerodynamic_data.alpha.value = char(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.CL.value',label);
%     Aircraft.Certification.Aerodynamic_data.CL.value = char(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.CD.value',label);
%     Aircraft.Certification.Aerodynamic_data.CD.value = char(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.CM.value',label);
%     Aircraft.Certification.Aerodynamic_data.CM.value = char(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Wing.S.value', label);
%     Aircraft.Geometry.Wing.S.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Wing.b.value', label);
%     Aircraft.Geometry.Wing.b.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Wing.AR.value', label);
%     Aircraft.Geometry.Wing.AR.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Weight.I_Level.W_maxTakeOff.value', label);
%     Aircraft.Weight.I_Level.W_maxTakeOff.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Weight.I_Level.W_maxTakeOff.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('first_sweep_value', label);
%     Aircraft.Geometry.Wing.sweep_first.value = str2double(table2array(value(p==1, 1)));  
%     Aircraft.Geometry.Wing.sweep_first.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('second_sweep_value', label);
%     Aircraft.Geometry.Wing.sweep_second.value = str2double(table2array(value(p==1, 1)));  
%     Aircraft.Geometry.Wing.sweep_second.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('third_sweep_value', label);
%     Aircraft.Geometry.Wing.sweep_third.value = str2double(table2array(value(p==1, 1)));  
%     Aircraft.Geometry.Wing.sweep_third.Attributes.unit = char(table2array(unit(p==1,1)));
% %     p = strcmp('Aircraft.Certification.Aerodynamic_data.xac.value', label);
% %     Aircraft.Geometry.General.xac.value = str2double(table2array(value(p==1,1)));
% %     p = strcmp('Aircraft.Certification.Aerodynamic_data.yac.value', label);
% %     Aircraft.Geometry.General.yac.value = str2double(table2array(value(p==1,1)));
% %     p = strcmp('Aircraft.Certification.Aerodynamic_data.zac.value', label);
% %     Aircraft.Geometry.General.zac.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.XAC_nondim', label);
%     Aircraft.Geometry.General.XAC_nondim.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.XCG_nondim', label);
%     Aircraft.Geometry.General.XCG_nondim.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.S.value', label);
%     Aircraft.Geometry.Horizontal.S.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.l.value', label);
%     Aircraft.Geometry.Horizontal.l.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Engine.Primary.Thrust_axes.value', label);
%     Aircraft.Geometry.Engine.Primary.Thrust_axes.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.CD_landing_gear.value', label);
%     Aircraft.Certification.Aerodynamic_data.CD_landing_gear.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.xcg.value', label);
%     Aircraft.Geometry.General.xcg.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.ycg.value', label);
%     Aircraft.Geometry.General.ycg.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.zcg.value', label);
%     Aircraft.Geometry.General.zcg.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.CD0.value', label);
%     Aircraft.Certification.Aerodynamic_data.CD0.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.e.value', label);
%     Aircraft.Certification.Aerodynamic_data.e.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.CL0.value', label);
%     Aircraft.Certification.Aerodynamic_data.CL0.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.CL_star.value', label);
%     Aircraft.Certification.Aerodynamic_data.CL_star.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.CM0.value', label);
%     Aircraft.Certification.Aerodynamic_data.CM0.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.CMCL.value', label);
%     Aircraft.Certification.Aerodynamic_data.CMCL.value = str2double(table2array(value(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.bcg.value', label);
%     Aircraft.Geometry.General.bcg.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.General.bcg.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.CM_landing_gear', label);
%     Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Certification.Aerodynamic_data.CM_landing_gear.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('Aircraft.Geometry.Wing.croot', label);
%     Aircraft.Geometry.Wing.croot.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Wing.croot.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Wing.ctip', label);
%     Aircraft.Geometry.Wing.ctip.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Wing.ctip.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Wing.twist_angle', label);
%     Aircraft.Geometry.Wing.twist_angle_first.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Wing.twist_angle_first.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Aileron.S', label); 
%     Aircraft.Geometry.Aileron.S.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Aileron.S.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Aileron.b', label);
%     Aircraft.Geometry.Aileron.b.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Aileron.b.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Aileron.ca', label);
%     Aircraft.Geometry.Aileron.ca.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Aileron.ca.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Aileron.cb', label);
%     Aircraft.Geometry.Aileron.cb.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Aileron.ca.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Aileron.y_iniziale', label);
%     Aircraft.Geometry.Aileron.y_inner.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Aileron.y_inner.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Aileron.y_finale', label);
%     Aircraft.Geometry.Aileron.y_outer.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Aileron.y_outer.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Aileron.eta_iniziale', label); 
%     Aircraft.Geometry.Aileron.eta_inner.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Aileron.eta_inner.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Aileron.eta_finale', label);
%     Aircraft.Geometry.Aileron.eta_outer.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Aileron.eta_outer.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Normal_Force_Curve_Slope_deg', label);
%     Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Fuselage.length', label);   
%     Aircraft.Geometry.Fuselage.length.value = str2double(table2array(value(p==1,1)));
%     Aircraft.Geometry.Fuselage.length.Attributes.unit = char(table2array(unit(p==1,1)));  
%     p = strcmp('Aircraft.Geometry.Fuselage.diameter', label);
%     Aircraft.Geometry.Fuselage.diameter.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Fuselage.diameter.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Certification.Aerodynamic_data.Ideal_cl', label);
%     Aircraft.Certification.Aerodynamic_data.Ideal_cl.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.Ideal_cl.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Wing.camberloc', label);
%     Aircraft.Geometry.Wing.camberloc.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.camberloc.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Wing.thickchord', label);    
%     Aircraft.Geometry.Wing.thickchord.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.thickchord.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Wing.xle', label); 
%     Aircraft.Geometry.Wing.xle.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.xle.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Wing.yle', label); 
%     Aircraft.Geometry.Wing.yle.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.yle.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Wing.zle', label); 
%     Aircraft.Geometry.Wing.zle.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.zle.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.camber', label);
%     Aircraft.Geometry.Horizontal.camber.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Horizontal.camber.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.location_of_camber', label);
%     Aircraft.Geometry.Horizontal.location_of_camber.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Horizontal.location_of_camber.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.thickchord', label);
%     Aircraft.Geometry.Horizontal.thickchord.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Horizontal.thickchord.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.twistvalue', label);
%     Aircraft.Geometry.Horizontal.twist.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Horizontal.twist.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.twistlocation', label);
%     Aircraft.Geometry.Horizontal.twistloc.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Horizontal.twistloc.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.zero_x', label);
%     Aircraft.Geometry.Horizontal.xloc0.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Horizontal.xloc0.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.xloc', label);
%     Aircraft.Geometry.Horizontal.xloc.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Horizontal.xloc.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.zloc', label);
%     Aircraft.Geometry.Horizontal.zloc.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Horizontal.zloc.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.b', label);
%     Aircraft.Geometry.Horizontal.b.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Horizontal.b.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.ctip', label);
%     Aircraft.Geometry.Horizontal.ctip.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Horizontal.ctip.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.croot', label);
%     Aircraft.Geometry.Horizontal.croot.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Horizontal.croot.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.sweep', label);
%     Aircraft.Geometry.Horizontal.sweep.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Horizontal.sweep.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.freccia_loc', label);
%     Aircraft.Geometry.Horizontal.sweeploc.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Horizontal.sweeploc.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.freccia_sec_loc', label);
%     Aircraft.Geometry.Horizontal.secondary_sweep_location.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Horizontal.secondary_sweep_location.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.dihedral', label);
%     Aircraft.Geometry.Horizontal.dihedral.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Horizontal.dihedral.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.Movable.eta_inner', label);
%     Aircraft.Geometry.Elevator.eta_inner.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.eta_inner.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.Movable.eta_outer', label);
%     Aircraft.Geometry.Elevator.eta_outer.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.eta_outer.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.Movable.cf_c_inner', label);
%     Aircraft.Geometry.Elevator.cf_c_inner.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.cf_c_inner.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Horizontal.Movable.cf_c_outer', label);
%     Aircraft.Geometry.Elevator.cf_c_outer.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.cf_c_outer.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Wing.sweep', label);
%     Aircraft.Geometry.Wing.sweep.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.sweep.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Wing.secondary_sweep_location', label);
%     Aircraft.Geometry.Wing.secondary_sweep_location.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.secondary_sweep_location.Attributes.unit = char(table2array(unit(p==1,1)));  
%     p = strcmp('Aircraft.Geometry.Wing.freccia_posizione', label);
%     Aircraft.Geometry.Wing.sweep_location.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.sweep_location.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Aircraft.Geometry.Wing.dihedral', label);
%     Aircraft.Geometry.Wing.dihedral.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.dihedral.Attributes.unit = char(table2array(unit(p==1,1)));  
%     p = strcmp('Aircraft.Certification.ISA_Condition.Sea_level', label);
%     Aircraft.Certification.ISA_Condition.Sea_level.Altitude.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.ISA_Condition.Sea_level.Altitude.Attributes.unit = char(table2array(unit(p==1,1)));  
%     p = strcmp('Operative_ceiling', label);
%     Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('Theoretical_ceiling', label);
%     Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Airloads_flag', label);
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Airload_case.Attributes.case = char(table2array(value(p==1,1)));
%     p = strcmp('Takeoff_power', label); 
%     Aircraft.Engine.Takeoff.Power.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Takeoff.Power.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Takeoff_rpm', label);
%     Aircraft.Engine.Takeoff.RPM.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Takeoff.RPM.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Max_Continous_power', label); 
%     Aircraft.Engine.Max_Continous.Power.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Max_Continous.Power.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Max_Continous_rpm', label); 
%     Aircraft.Engine.Max_Continous.RPM.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Max_Continous.RPM.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Correction_factor_flag1', label);
%     Aircraft.Engine.Correction_factor.Attributes.flag1 = char(table2array(value(p==1, 1)));
%     p = strcmp('Correction_factor_flag2', label);
%     Aircraft.Engine.Correction_factor.Attributes.flag2 = str2double(table2array(value(p==1, 1)));
%     p = strcmp('Reduction_ratio', label);
%     Aircraft.Engine.Reduction_ratio.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Reduction_ratio.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Engine.Limit_side_load', label);
%     Aircraft.Engine.Limit_side_load.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Limit_side_load.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Engine_mount_mass', label);
%     Aircraft.Engine.Engine_mount_mass.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Engine_mount_mass.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Engine_accessories', label);
%     Aircraft.Engine.Engine_accessories_mass.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Engine_accessories_mass.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Propeller_spinner', label);
%     Aircraft.Engine.Propeller_spinner_mass.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Propeller_spinner_mass.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Propeller_polarmoment', label);
%     Aircraft.Engine.Propeller_polar_moment.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Engine.Propeller_polar_moment.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Engine_pitch_speed', label);
%     Aircraft.Engine.Pitch_speed.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Pitch_speed.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('Engine_yaw_speed', label);
%     Aircraft.Engine.Yaw_speed.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Yaw_speed.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('Propeller_blade_number', label);
%     Aircraft.Engine.Propeller_blade_number.value = str2double(table2array(value(p==1, 1)));
%     p = strcmp('Engine_normal_load_factor', label);
%     Aircraft.Engine.Engine_normal_load_factor.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Engine.Engine_normal_load_factor.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('C_h_delta', label);
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_delta_rad.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_delta_rad.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('C_h_alfa', label);
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_alfa_rad.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_alfa_rad.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('CL_MAX_TAKEOFF', label);
%     Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_takeoff.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_takeoff.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('CL_MAX_LANDING', label);
%     Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_landing.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_landing.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('S_elevator', label);
%     Aircraft.Geometry.Elevator.S.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.S.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('c_elevator', label);
%     Aircraft.Geometry.Elevator.chord.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.chord.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('chord_ratio_elevator', label);
%     Aircraft.Geometry.Elevator.chord_ratio_ce_c.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.chord_ratio_ce_c.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('elevator_overhang', label);
%     Aircraft.Geometry.Elevator.overhang.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.overhang.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('span_ratio_elevator', label);
%     Aircraft.Geometry.Elevator.span_ratio.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.span_ratio.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('S_hinge_elevator', label);
%     Aircraft.Geometry.Elevator.S_hinge.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Elevator.S_hinge.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Chdeltaelevator', label);
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_delta_rad.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_delta_rad.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Chalfaelevator', label);
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_alfa_rad.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Elevator.C_h_alfa_rad.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('S_vertical', label);
%     Aircraft.Geometry.Vertical.S.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Vertical.S.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('chord_vertical', label);
%     Aircraft.Geometry.Vertical.chord.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Vertical.chord.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('S_rudder', label);
%     Aircraft.Geometry.Rudder.S.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Rudder.S.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('chord_rudder', label);
%     Aircraft.Geometry.Rudder.chord.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Rudder.chord.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('chord_ratio_rudder_cf_c', label);
%     Aircraft.Geometry.Rudder.chord_ratio_cf_c.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Rudder.chord_ratio_cf_c.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('overhang_rudder', label);
%     Aircraft.Geometry.Rudder.overhang.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Rudder.overhang.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('span_ratio_rudder', label);
%     Aircraft.Geometry.Rudder.span_ratio.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Rudder.span_ratio.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('max_deflection_rudder', label);
%     Aircraft.Geometry.Rudder.max_deflection.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Rudder.max_deflection.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('Chdeltarudder', label);
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_delta_rad.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_delta_rad.Attributes.unit = char(table2array(unit(p==1,1)));  
%     p = strcmp('Chalfarudder', label);
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_alfa_rad.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Certification.Aerodynamic_data.Hinge_moments.Rudder.C_h_alfa_rad.Attributes.unit = char(table2array(unit(p==1,1)));    
%     p = strcmp('ca_c_inner', label);
%     Aircraft.Geometry.Aileron.ca_c_root.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Aileron.ca_c_root.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('ca_c_outer', label);
%     Aircraft.Geometry.Aileron.ca_c_tip.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Aileron.ca_c_tip.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('c_root_horizontal', label);
%     Aircraft.Geometry.Horizontal.croot.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Horizontal.croot.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('c_tip_horizontal', label);
%     Aircraft.Geometry.Horizontal.ctip.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Horizontal.ctip.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('ce_c_horizontal_root', label);
%     Aircraft.Geometry.Horizontal.ce_c_root.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Horizontal.ce_c_root.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('ce_c_horizontal_tip', label);
%     Aircraft.Geometry.Horizontal.ce_c_tip.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Horizontal.ce_c_tip.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('c_root_vertical', label);
%     Aircraft.Geometry.Vertical.croot.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Vertical.croot.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('c_tip_vertical', label);
%     Aircraft.Geometry.Vertical.ctip.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Vertical.ctip.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('cr_c_root', label);
%     Aircraft.Geometry.Rudder.cr_c_root.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Rudder.cr_c_root.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('cr_c_tip', label);
%     Aircraft.Geometry.Rudder.cr_c_tip.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Rudder.cr_c_tip.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('eta_inner_rudder', label);
%     Aircraft.Geometry.Rudder.eta_inner.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Rudder.eta_inner.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('eta_outer_rudder', label);
%     Aircraft.Geometry.Rudder.eta_outer.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Rudder.eta_outer.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('c_kink_one', label);
%     Aircraft.Geometry.Wing.chord_kink_one.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Wing.chord_kink_one.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('c_kink_two', label);
%     Aircraft.Geometry.Wing.chord_kink_two.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Wing.chord_kink_two.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('wing_type', label);
%     Aircraft.Geometry.Wing.type.value = char(table2array(value(p==1,1)));
%     p = strcmp('vertical_flag', label);
%     Aircraft.Geometry.Vertical.empennage_flag.value = char(table2array(value(p==1,1)));
%     p = strcmp('panel_span1', label);
%     Aircraft.Geometry.Wing.panel_span1.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Wing.panel_span1.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('panel_span2', label);
%     Aircraft.Geometry.Wing.panel_span2.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Wing.panel_span2.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('panel_span3', label);
%     Aircraft.Geometry.Wing.panel_span3.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Wing.panel_span3.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('horizontal_tail_damping_factor', label);
%     Aircraft.Certification.Aerodynamic_data.Horizontal.damping_factor.value = str2double(table2array(value(p==1, 1)));  
%     Aircraft.Certification.Aerodynamic_data.Horizontal.damping_factor.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('airfoil_first_panel', label);
%     Aircraft.Certification.Aerodynamic_data.airfoil_first_panel.value = char(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.airfoil_first_panel.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('airfoil_sec_panel', label);
%     Aircraft.Certification.Aerodynamic_data.airfoil_second_panel.value = char(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.airfoil_second_panel.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('airfoil_third_panel', label);
%     Aircraft.Certification.Aerodynamic_data.airfoil_third_panel.value = char(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.airfoil_third_panel.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('airfoil_fourth_panel', label);
%     Aircraft.Certification.Aerodynamic_data.airfoil_fourth_panel.value = char(table2array(value(p==1, 1)));
%     Aircraft.Certification.Aerodynamic_data.airfoil_fourth_panel.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('sweep_first', label);
%     Aircraft.Geometry.Wing.sweep_first.value = str2double(table2array(value(p==1, 1)));  
%     Aircraft.Geometry.Wing.sweep_first.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('sweep_second', label);
%     Aircraft.Geometry.Wing.sweep_second.value = str2double(table2array(value(p==1, 1)));  
%     Aircraft.Geometry.Wing.sweep_second.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('sweep_third', label);
%     Aircraft.Geometry.Wing.sweep_third.value = str2double(table2array(value(p==1, 1)));  
%     Aircraft.Geometry.Wing.sweep_third.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('dihedral_first', label);
%     Aircraft.Geometry.Wing.dihedral_first.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.dihedral_first.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('dihedral_second', label);
%     Aircraft.Geometry.Wing.dihedral_second.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.dihedral_second.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('dihedral_third', label);
%     Aircraft.Geometry.Wing.dihedral_third.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.dihedral_third.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('twist_angle_second', label);
%     Aircraft.Geometry.Wing.twist_angle_second.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.twist_angle_second.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('twist_angle_third', label);
%     Aircraft.Geometry.Wing.twist_angle_third.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.twist_angle_third.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('twist_angle_fourth', label);
%     Aircraft.Geometry.Wing.twist_angle_fourth.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Wing.twist_angle_fourth.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('vertical_sweep', label);
%     Aircraft.Geometry.Vertical.sweep.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Vertical.sweep.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('sweep_loc_vert', label);
%     Aircraft.Geometry.Vertical.sweeploc.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Vertical.sweeploc.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('sec_sweeploc_vert', label);
%     Aircraft.Geometry.Vertical.secsweeploc.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Vertical.secsweeploc.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('vertical_dihedral', label);
%     Aircraft.Geometry.Vertical.dihedral.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Vertical.dihedral.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('vertical_twist', label);
%     Aircraft.Geometry.Vertical.twist.value       = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Vertical.twist.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('ver_twist_loc', label);
%     Aircraft.Geometry.Vertical.twistloc.value       = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Vertical.twistloc.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('vertical_x_leadingedge', label);
%     Aircraft.Geometry.Vertical.xle.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Vertical.xle.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('vertical_croot', label);
%     Aircraft.Geometry.Vertical.croot.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Vertical.croot.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('vertical_ctip', label);
%     Aircraft.Geometry.Vertical.ctip.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Vertical.ctip.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('vertical_xtip_leadingedge', label);   
%     Aircraft.Geometry.Vertical.xtip_le.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Vertical.xtip_le.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('vertical_span', label);  
%     Aircraft.Geometry.Vertical.b.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Vertical.b.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('vertical_z_position', label);
%     Aircraft.Geometry.Vertical.zpos.value = str2double(table2array(value(p==1, 1)));   %1.0; 
%     Aircraft.Geometry.Vertical.zpos.Attributes.unit = char(table2array(unit(p==1,1))); % "% of fuselage diameter";
%     p = strcmp('flaps_eta_iniziale', label);
%     Aircraft.Geometry.Flaps.eta_inner.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Flaps.eta_inner.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('flaps_eta_finale', label);    
%     Aircraft.Geometry.Flaps.eta_outer.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Flaps.eta_outer.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('flaps_y_iniziale', label);    
%     Aircraft.Geometry.Flaps.y_inner.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Flaps.y_inner.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('flaps_y_finale', label);
%     Aircraft.Geometry.Flaps.y_outer.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Flaps.y_outer.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('flaps_span_b', label);
%     Aircraft.Geometry.Flaps.b.value = str2double(table2array(value(p==1, 1))); 
%     Aircraft.Geometry.Flaps.b.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('cf_c_root_flap', label); 
%     Aircraft.Geometry.Flaps.cf_c_root.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Flaps.cf_c_root.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('cf_c_tip_flap', label); 
%     Aircraft.Geometry.Flaps.cf_c_tip.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Flaps.cf_c_tip.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('croot_flap', label);
%     Aircraft.Geometry.Flaps.croot.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Flaps.croot.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('ctip_flap', label);
%     Aircraft.Geometry.Flaps.ctip.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Flaps.ctip.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('cf_flap', label);
%     Aircraft.Geometry.Flaps.cf.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Flaps.cf.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('flap_surface_S', label);
%     Aircraft.Geometry.Flaps.S.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Geometry.Flaps.S.Attributes.unit = char(table2array(unit(p==1,1))); 
%     % -------------------------------------------------------------------------- 
%     p = strcmp('vertical_yaw_angle', label);
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('verticaloads_dr', label);
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.dr.value = char(table2array(value(p==1, 1)));
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.dr.value = str2num(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.dr.value);
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.dr.Attributes.unit = char(table2array(unit(p==1,1)));    
%     p = strcmp('vertical_CY_vector', label);
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_vector.value = char(table2array(value(p==1, 1)));  % [0.00003615686 0.01291153]';
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_vector.value = str2num(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_vector.value);
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_vector.Attributes.unit = char(table2array(unit(p==1,1))); 
%     p = strcmp('vertical_CY0', label);
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_0.value = str2double(table2array(value(p==1, 1)));
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_0.Attributes.unit = char(table2array(unit(p==1,1)));
%     p = strcmp('vertical_CYdr', label);
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CYdr.value = str2double(table2array(value(p==1, 1)));   % 0.000644;
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CYdr.Attributes.unit = char(table2array(unit(p==1,1))); % "1/deg";
%     p = strcmp('maximum_delta_rudder', label);
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.value = str2double(table2array(value(p==1, 1)));   % 30.0;
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.Attributes.unit = char(table2array(unit(p==1,1))); % "deg";
%     % CASE 2
%     p = strcmp('CY_vertical_tailplane', label);
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_VTP.value = str2double(table2array(value(p==1, 1)));   % 0.0245;
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_VTP.Attributes.unit = char(table2array(unit(p==1,1))); % "Non dimensional";
%     % CASE 3
%     p = strcmp('case3_CY_vert_tp', label);
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_VTP.value = str2double(table2array(value(p==1, 1)));   % 0.0233;
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_VTP.Attributes.unit = char(table2array(unit(p==1,1))); % "Non dimensional";
%     p = strcmp('vert_tp_K', label);
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.value = str2double(table2array(value(p==1, 1)));   % 0.3; 
%     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.Attributes.unit = char(table2array(unit(p==1,1))); % "m"; 
% 
% end
% 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value = NaN;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.Attributes.unit = "Positive g";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.value = NaN;
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.Attributes.unit = "Negative g";
% Aircraft.Geometry.Wing.twist_angle.value = Aircraft.Geometry.Wing.twist_angle_first.value;
% Aircraft.Geometry.Wing.twist_angle.Attributes.unit = "deg";
% 
% % 
% %          'Aircraft.Geometry.Wing.mac', ...
% %          'Aircraft.Certification.Aerodynamic_data.xac.value', ...
% %          'Aircraft.Certification.Aerodynamic_data.yac.value', ...
% %          'Aircraft.Certification.Aerodynamic_data.zac.value', ...
% %% Derived values (also inputs)
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value = calcn(Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value);
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.Attributes.unit = "Positive g";
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.value = calcn(Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value);
% % Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.Attributes.unit = "Negative g";
