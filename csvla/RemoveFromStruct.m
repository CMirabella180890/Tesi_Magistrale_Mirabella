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
        Aircraft.Certification.Aerodynamic_data.Aileron.Max_deflection.Attributes.unit = "degrees";
        field = 'Max_deflection';
        Aircraft.Geometry.Aileron = rmfield(Aircraft.Geometry.Aileron, field);
        % -----------------------------------------------------------------
        Aircraft.Certification.Aerodynamic_data.Elevator.Max_deflection.value = Aircraft.Geometry.Elevator.max_deflection.value;
        Aircraft.Certification.Aerodynamic_data.Elevator.Max_deflection.Attributes.unit = "degrees";
        Aircraft.Certification.Aerodynamic_data.Elevator.total_deflection_time.value = Aircraft.Geometry.Elevator.total_deflection_time.value;
        Aircraft.Certification.Aerodynamic_data.Elevator.total_deflection_time.Attributes.unit = "seconds";
        field = 'max_deflection';
        Aircraft.Geometry.Elevator = rmfield(Aircraft.Geometry.Elevator, field);
        field = 'total_deflection_time';
        Aircraft.Geometry.Elevator = rmfield(Aircraft.Geometry.Elevator, field);
        % -----------------------------------------------------------------
        Aircraft.Certification.Aerodynamic_data.Rudder.Max_deflection.value = Aircraft.Geometry.Rudder.max_deflection.value;
        Aircraft.Certification.Aerodynamic_data.Rudder.Max_deflection.Attributes.unit = "degrees";
        field = 'max_deflection';
        Aircraft.Geometry.Rudder = rmfield(Aircraft.Geometry.Rudder, field);
        % -----------------------------------------------------------------
        Aircraft.Geometry.General.xcg_divided_by_mac.value = 0.25;
        Aircraft.Geometry.General.xcg_divided_by_mac.Attributes.unit = "% MAC"; % Measured from the aircraft nose
        Aircraft.Geometry.General.ycg_divided_by_mac.value = 0.0;
        Aircraft.Geometry.General.ycg_divided_by_mac.Attributes.unit = "% MAC"; % Measured from the aircraft nose
        Aircraft.Geometry.General.zcg_divided_by_mac.value = 0.0;
        Aircraft.Geometry.General.zcg_divided_by_mac.Attributes.unit = "% MAC"; % Measured from the aircraft nose
        % -----------------------------------------------------------------
        fields = {'xac', 'yac', 'zac'};
        Aircraft.Geometry.General = rmfield(Aircraft.Geometry.General, fields);
        % -----------------------------------------------------------------
    case 'SCHRENK'
        % -----------------------------------------------------------------
        Aircraft.Geometry.General.xcg_divided_by_mac.value = 0.25;
        Aircraft.Geometry.General.xcg_divided_by_mac.Attributes.unit = "% MAC"; % Measured from the aircraft nose
        Aircraft.Geometry.General.ycg_divided_by_mac.value = 0.0;
        Aircraft.Geometry.General.ycg_divided_by_mac.Attributes.unit = "% MAC"; % Measured from the aircraft nose
        Aircraft.Geometry.General.zcg_divided_by_mac.value = 0.0;
        Aircraft.Geometry.General.zcg_divided_by_mac.Attributes.unit = "% MAC"; % Measured from the aircraft nose
        % -----------------------------------------------------------------
        fields = {'xac', 'yac', 'zac'};
        Aircraft.Geometry.General = rmfield(Aircraft.Geometry.General, fields);
        % -----------------------------------------------------------------
end
