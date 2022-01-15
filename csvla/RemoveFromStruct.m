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
    case 'SCHRENK'
        % -----------------------------------------------------------------
end
