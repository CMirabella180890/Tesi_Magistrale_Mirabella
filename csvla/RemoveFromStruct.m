%% REMOVE FIELDS INSIDE THE AIRCRAFT STRUCT VARIABLE 

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% AILERON
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% =================================================================
switch (Aircraft.Certification.Regulation.SubpartC.Flightloads.Airload_case.Attributes.case)
    case 'OPEN VSP'
        % -----------------------------------------------------------------
        field = 'y_span';
        Aircraft.Geometry.Aileron = rmfield(Aircraft.Geometry.Aileron, field);
        % -----------------------------------------------------------------
        Aircraft.Certification.Aerodynamic_data.Aileron.Max_deflection.value = Aircraft.Geometry.Aileron.Max_deflection.value;
        Aircraft.Certification.Aerodynamic_data.Aileron.Max_deflection.Attributes.unit = "deg";
        field = 'Max_deflection';
        Aircraft.Geometry.Aileron = rmfield(Aircraft.Geometry.Aileron, field);
        % -----------------------------------------------------------------
        Aircraft.Certification.Aerodynamic_data.Elevator.Max_deflection.value = Aircraft.Geometry.Elevator.max_deflection.value;
        Aircraft.Certification.Aerodynamic_data.Elevator.Max_deflection.Attributes.unit = "deg";
        Aircraft.Certification.Aerodynamic_data.Elevator.total_deflection_time.value = Aircraft.Geometry.Elevator.total_deflection_time.value;
        Aircraft.Certification.Aerodynamic_data.Elevator.total_deflection_time.Attributes.unit = "s";
        field = 'max_deflection';
        Aircraft.Geometry.Elevator = rmfield(Aircraft.Geometry.Elevator, field);
        field = 'total_deflection_time';
        Aircraft.Geometry.Elevator = rmfield(Aircraft.Geometry.Elevator, field);
        % -----------------------------------------------------------------
        Aircraft.Certification.Aerodynamic_data.Rudder.Max_deflection.value = Aircraft.Geometry.Rudder.max_deflection.value;
        Aircraft.Certification.Aerodynamic_data.Rudder.Max_deflection.Attributes.unit = "deg";
        field = 'max_deflection';
        Aircraft.Geometry.Rudder = rmfield(Aircraft.Geometry.Rudder, field);
        % -----------------------------------------------------------------
        Aircraft.Geometry.General.xcg_divided_by_mac.value = 0.25;
        Aircraft.Geometry.General.xcg_divided_by_mac.Attributes.unit = "% MAC"; % Measured from the aircraft nose
        Aircraft.Geometry.General.ycg_divided_by_mac.value = 0.0;
        Aircraft.Geometry.General.ycg_divided_by_mac.Attributes.unit = "% MAC"; % Measured from the aircraft nose
        Aircraft.Geometry.General.zcg_divided_by_mac.value = 0.0;
        Aircraft.Geometry.General.zcg_divided_by_mac.Attributes.unit = "% MAC"; % Measured from the aircraft nose
        Aircraft.Geometry.General.xac_nondim.value = Aircraft.Geometry.General.XAC_nondim.value;
        Aircraft.Geometry.General.xac_nondim.Attributes.unit = "Non dimensional";
        Aircraft.Geometry.General.xcg_nondim.value = Aircraft.Geometry.General.XCG_nondim.value;
        Aircraft.Geometry.General.xcg_nondim.Attributes.unit = "Non dimensional";
        % -----------------------------------------------------------------
        fields = {'xac', 'yac', 'zac', 'X_cg', 'XAC_nondim', 'XCG_nondim'};
        Aircraft.Geometry.General = rmfield(Aircraft.Geometry.General, fields);
        % -----------------------------------------------------------------
        Aircraft.Propeller.Propeller_spinner_mass.value = Aircraft.Engine.Propeller_spinner_mass.value;  
        Aircraft.Propeller.Propeller_spinner_mass.Attributes.unit = "kg";
        Aircraft.Propeller.Propeller_polar_moment.value = Aircraft.Engine.Propeller_polar_moment.value;  
        Aircraft.Propeller.Propeller_polar_moment.Attributes.unit = "kg*m^2";
        Aircraft.Propeller.Propeller_blade_number.value = Aircraft.Engine.Propeller_blade_number.value;  
        Aircraft.Propeller.Propeller_blade_number.Attributes.unit = "Pure number";
        % -----------------------------------------------------------------
        fields = {'Propeller_spinner_mass', 'Propeller_polar_moment', 'Propeller_blade_number'};
        Aircraft.Engine = rmfield(Aircraft.Engine, fields);
        % -----------------------------------------------------------------
    case 'SCHRENK'
        % -----------------------------------------------------------------
        Aircraft.Geometry.General.xcg_divided_by_mac.value = 0.25;
        Aircraft.Geometry.General.xcg_divided_by_mac.Attributes.unit = "% MAC"; % Measured from the aircraft nose
        Aircraft.Geometry.General.ycg_divided_by_mac.value = 0.0;
        Aircraft.Geometry.General.ycg_divided_by_mac.Attributes.unit = "% MAC"; % Measured from the aircraft nose
        Aircraft.Geometry.General.zcg_divided_by_mac.value = 0.0;
        Aircraft.Geometry.General.zcg_divided_by_mac.Attributes.unit = "% MAC"; % Measured from the aircraft nose
        Aircraft.Geometry.General.xac_nondim.value = Aircraft.Geometry.General.XAC_nondim.value;
        Aircraft.Geometry.General.xac_nondim.Attributes.unit = "Non dimensional";
        Aircraft.Geometry.General.xcg_nondim.value = Aircraft.Geometry.General.XCG_nondim.value;
        Aircraft.Geometry.General.xcg_nondim.Attributes.unit = "Non dimensional";
        % -----------------------------------------------------------------
        fields = {'xac', 'yac', 'zac', 'X_cg', 'XAC_nondim', 'XCG_nondim'};
        Aircraft.Geometry.General = rmfield(Aircraft.Geometry.General, fields);
        % -----------------------------------------------------------------
        Aircraft.Propeller.Propeller_spinner_mass.value = Aircraft.Engine.Propeller_spinner_mass.value;  
        Aircraft.Propeller.Propeller_spinner_mass.Attributes.unit = "kg";
        Aircraft.Propeller.Propeller_polar_moment.value = Aircraft.Engine.Propeller_polar_moment.value;  
        Aircraft.Propeller.Propeller_polar_moment.Attributes.unit = "kg*m^2";
        Aircraft.Propeller.Propeller_blade_number.value = Aircraft.Engine.Propeller_blade_number.value;  
        Aircraft.Propeller.Propeller_blade_number.Attributes.unit = "Pure number";
        % -----------------------------------------------------------------
        fields = {'Propeller_spinner_mass', 'Propeller_polar_moment', 'Propeller_blade_number'};
        Aircraft.Engine = rmfield(Aircraft.Engine, fields);
        % -----------------------------------------------------------------
end
