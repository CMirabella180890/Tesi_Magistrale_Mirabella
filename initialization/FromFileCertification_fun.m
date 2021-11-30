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
         'Aircraft.Certification.ISA_Condition.rho0', ...
         'Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient', ...
         'Aircraft.Certification.Aerodynamic_data.Min_Lift_Coefficient', ...
         'Aircraft.Certification.Regulation.SubpartC.Flightloads.Max_Continous_Power_Speed_VH', ...
         'Aircraft.Certification.Regulation.SubpartC.Flightloads.Min_Design_Cruise_Speed', ...
         'Aircraft.Geometry.Wing.mac', ...
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
         'Aircraft.Certification.Aerodynamic_data.xac.value', ...
         'Aircraft.Certification.Aerodynamic_data.yac.value', ...
         'Aircraft.Certification.Aerodynamic_data.zac.value', ...
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
         'Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value', ...
         'Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value', ...
         'Aircraft.Certification.Aerodynamic_data.CL0.value', ...
         'Aircraft.Certification.Aerodynamic_data.CL_star.value', ...
         'Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a.value', ...
         'Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b.value', ...
         'Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c.value', ...
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
         'Aircraft.Geometry.Wing.dihedral'};

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
    p = strcmp('Aircraft.Certification.ISA_Condition.rho0',label);
    Aircraft.Certification.ISA_Condition.rho0.value = str2double(table2array(value(p==1,1))); 
    Aircraft.Certification.ISA_Condition.rho0.Attributes.unit = char(table2array(unit(p==1,1)));     
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
    p = strcmp('Aircraft.Geometry.Wing.mac',label);
    Aircraft.Geometry.Wing.mac.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Wing.mac.Attributes.unit = char(table2array(unit(p==1,1))); 
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
    p = strcmp('Aircraft.Certification.Aerodynamic_data.xac.value', label);
    Aircraft.Certification.Aerodynamic_data.xac.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.yac.value', label);
    Aircraft.Certification.Aerodynamic_data.yac.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.zac.value', label);
    Aircraft.Certification.Aerodynamic_data.zac.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.XAC_nondim', label);
    Aircraft.Certification.Aerodynamic_data.XAC_nondim.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.XCG_nondim', label);
    Aircraft.Certification.Aerodynamic_data.XCG_nondim.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.S.value', label);
    Aircraft.Geometry.Horizontal.S.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Geometry.Horizontal.l.value', label);
    Aircraft.Geometry.Horizontal.l.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Geometry.Engine.Primary.Thrust_axes.value', label);
    Aircraft.Geometry.Engine.Primary.Thrust_axes.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.CD_landing_gear.value', label);
    Aircraft.Certification.Aerodynamic_data.CD_landing_gear.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.xcg.value', label);
    Aircraft.Certification.Aerodynamic_data.xcg.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.ycg.value', label);
    Aircraft.Certification.Aerodynamic_data.ycg.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.zcg.value', label);
    Aircraft.Certification.Aerodynamic_data.zcg.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.CD0.value', label);
    Aircraft.Certification.Aerodynamic_data.CD0.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.e.value', label);
    Aircraft.Certification.Aerodynamic_data.e.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value', label);
    Aircraft.Certification.Aerodynamic_data.Pol_coeff_k1.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value', label);
    Aircraft.Certification.Aerodynamic_data.Pol_coeff_k2.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.CL0.value', label);
    Aircraft.Certification.Aerodynamic_data.CL0.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.CL_star.value', label);
    Aircraft.Certification.Aerodynamic_data.CL_star.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a.value', label);
    Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a.value = str2double(table2array(value(p==1,1))); 
    p = strcmp('Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b.value', label);
    Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c.value', label);
    Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.CM0.value', label);
    Aircraft.Certification.Aerodynamic_data.CM0.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.CMCL.value', label);
    Aircraft.Certification.Aerodynamic_data.CMCL.value = str2double(table2array(value(p==1,1)));
    p = strcmp('Aircraft.Certification.Aerodynamic_data.bcg.value', label);
    Aircraft.Certification.Aerodynamic_data.bcg.value = str2double(table2array(value(p==1,1)));
    Aircraft.Certification.Aerodynamic_data.bcg.Attributes.unit = char(table2array(unit(p==1,1))); 
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
    Aircraft.Geometry.Wing.twist_angle.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Wing.twist_angle.Attributes.unit = char(table2array(unit(p==1,1)));
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
    Aircraft.Geometry.Aileron.y_iniziale.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Aileron.y_iniziale.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Aileron.y_finale', label);
    Aircraft.Geometry.Aileron.y_finale.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Aileron.y_finale.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Aileron.eta_iniziale', label); 
    Aircraft.Geometry.Aileron.eta_iniziale.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Aileron.eta_iniziale.Attributes.unit = char(table2array(unit(p==1,1)));
    p = strcmp('Aircraft.Geometry.Aileron.eta_finale', label);
    Aircraft.Geometry.Aileron.eta_finale.value = str2double(table2array(value(p==1,1)));
    Aircraft.Geometry.Aileron.eta_finale.Attributes.unit = char(table2array(unit(p==1,1)));
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
end

Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.Attributes.unit = "Positive g";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.value = NaN;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.Attributes.unit = "Negative g";
%% Derived values (also inputs)
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value = calcn(Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.Attributes.unit = "Positive g";
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.value = calcn(Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value);
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.Attributes.unit = "Negative g";
