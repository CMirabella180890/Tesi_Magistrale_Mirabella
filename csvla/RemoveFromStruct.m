%% REMOVE FIELDS INSIDE THE AIRCRAFT STRUCT VARIABLE 

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% AILERON
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
field = 'y_span';
Aircraft.Geometry.Aileron = rmfield(Aircraft.Geometry.Aileron, field);
