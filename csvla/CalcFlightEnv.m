%% FLIGHT ENVELOPE
% x = calcn(obj, nmax) - from csvla.m
% This function defines a vector with load factor values between two pre-
% scribed limits. Check the class file csvla.m to have a complete
% documentation.
% Positive load factor values
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value = calcn(obj, Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value);
% Negative load factor values
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.value = calcn(obj, Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value);
% -------------------------------------------------------------------------
%% x = calcvs(obj, rho, WingLoading, MaxLiftCoeff, PositiveLoadFactors)
% This function defines a vector with stall airspeed for the chosen
% aircraft, within the precribed range of load factors. Check the
% class file csvla.m to have a complete documentation.
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_VS.value = calcvs(obj, Aircraft.Certification.ISA_Condition.rho.value, ...                            % Standard atmosphere density
                                                                                 Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value, ...                % Wing Loading in SI units 
                                                                                 Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value, ...              % Maximum Lift coefficient
                                                                                 Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value); % A vector of load factors
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_VS.Attributes.unit = "m/s";                                                                              
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_VS.value = calcvs(obj, Aircraft.Certification.ISA_Condition.rho.value, ...                             % Standard atmosphere density
                                                                                  Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value, ...                % Wing Loading in SI units 
                                                                                  Aircraft.Certification.Aerodynamic_data.Min_Lift_Coefficient.value, ...              % Maximum Lift coefficient
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.value); % A vector of load factors
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_VS.Attributes.unit = "m/s"; 
% -------------------------------------------------------------------------
%% x = calcvc(obj, WingLoading, MaxContinuousPowerSpeedVH)
% This function identifies (following CS-VLA airworthiness reg.)
% maximum cruise speed (Point C) for flight envelope calculations. 
% To have a complete documentation check the class file csvla.m
% VH design speed for max continous power: this airspeed is not available
% but must be known. From CS - VLA Airworthiness rules
Aircraft.Certification.Regulation.SubpartC.Flightloads.Cruise_Speed_VC.value = calcvc(obj, Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Cruise_Speed_VC.Attributes.unit = "m/s";
% -------------------------------------------------------------------------
%% x = calcvd(obj, MinDesignCruiseSpeed, CruiseSpeedVC)
% This function identifies (following CS-VLA airworthiness reg.)
% the maximum dive speed (Point D) for flight envelope
% calculations. To have a complete documentation check the class
% file csvla.m
Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value = calcvd(obj, Aircraft.Certification.Regulation.SubpartC.Flightloads.Min_Design_Cruise_Speed.value, ... % Min design cruise speed 
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Cruise_Speed_VC.value);            % Cruise speed from previous calculations
Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.Attributes.unit = "m/s";
% -------------------------------------------------------------------------
%% Point G definition
% Assign speed at Point G equals to the speed at Point D
Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VE.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value; % Speed at points G and D are equal 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VE.Attributes.unit = "m/s";
% -------------------------------------------------------------------------
%% Point A definition
% Assign speed at Point A (Maneuver point) equals to the maximum
% permissible positive load factor stall speed. 

% index = dsearchn(Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value);
tol = 1E-3;
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value)
    if abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value - Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value(i)) < tol
        index = i;               
    end
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_VS.value(index);
clear index
Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.Attributes.unit = "m/s";
% -------------------------------------------------------------------------
%% Point G definition
% Assign speed at Point G (Neg. Maneuver point) equals to the
% maximum permissible negative load factor stall speed
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.value)
    if abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value - Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.value(i)) < tol
        temp = i;               
    end
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_Design_manoeuvring_speed_VG.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_VS.value(temp);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_Design_manoeuvring_speed_VG.Attributes.unit = "m/s";
% -------------------------------------------------------------------------

% POSITIVE STALL SPEED 
Aircraft.Certification.Regulation.SubpartC.Flightloads.print_positive_vs.value = calcvs(obj, Aircraft.Certification.ISA_Condition.rho.value, ...         % Standard atmosphere density
                                                                                 Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value, ...   % Wing Loading in SI units 
                                                                                 Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value, ... % Maximum Lift coefficient
                                                                                 1.0);                                                                   % A vector of load factors
Aircraft.Certification.Regulation.SubpartC.Flightloads.print_positive_vs.Attributes.unit = "m/s";

% NEGATIVE STALL SPEED 
Aircraft.Certification.Regulation.SubpartC.Flightloads.print_negative_vs.value = calcvs(obj, Aircraft.Certification.ISA_Condition.rho.value, ...         % Standard atmosphere density
                                                                                 Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value, ...   % Wing Loading in SI units 
                                                                                 Aircraft.Certification.Aerodynamic_data.Max_Inverted_Lift_Coefficient.value, ... % Maximum Lift coefficient
                                                                                 1.0);                                                                   % A vector of load factors
Aircraft.Certification.Regulation.SubpartC.Flightloads.print_negative_vs.Attributes.unit = "m/s";


%% INPUT TRACKING - VN DIAGRAM 
% A possible way to track inputs for the various data will be provided
% inside the .txt file used as a log for the program.

disp(" ")
disp(" ++++ INPUT TO V - N DIAGRAM ++++");
% Horizontal tail loads increments
Data1 = [  Aircraft.Certification.Regulation.SubpartC.Flightloads.print_positive_vs.value, ...
           Aircraft.Certification.Regulation.SubpartC.Flightloads.print_negative_vs.value, ...
           Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.value, ... % Max positive value of load factors
           Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_Design_manoeuvring_speed_VG.value];                 % VG = VD on the negative side of V - n diagram
disp(" ++++++++++ DATA USED TO PLOT V - N DIAGRAM ++++++++++ ")
format = ' %6.6f          %6.6f          %6.6f          %6.6f\n';
label  = ' VS+                VS-                 VA                VG\n';
fprintf(label);
fprintf(format, Data1.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

disp(" ")
% Horizontal tail loads increments
Data1 = [  Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value, ...                       % Max positive value of load factors
           Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value, ...                       % Min (negative) value of load factors
           Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value, ...              % Max dive speed from V - n diagram
           Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VE.value];                 % VG = VD on the negative side of V - n diagram
disp(" ++++++++++ DATA USED TO PLOT V - N DIAGRAM ++++++++++ ")
format = ' %6.6f          %6.6f          %6.6f          %6.6f\n';
label  = ' nmax                nmin                 VD                VE\n';
fprintf(label);
fprintf(format, Data1.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

%% OUTPUT 
% x = V_n_diagram(npos, nneg, nmax, nmin, VSpos, VSneg, VD, VG, Reg, Aircraft_name)
% This function construct, plot and save in various format the
% flight envelope of the aircraft, following CS-VLA airworthiness
% prescription. To have a complete documentation check the class
% file csvla.m

disp(" ")
disp(" ++++ FIGURE 1 - FLIGHT ENVELOPE DIAGRAM ++++ ");
Aircraft.Certification.Regulation.SubpartC.Flightloads.V_n_diagram.value = V_n_diagram(obj, Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value, ... % Positive load factors
                                                                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.value, ...      % Negative load factors
                                                                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value, ...                       % Max positive value of load factors
                                                                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value, ...                       % Min (negative) value of load factors
                                                                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_VS.value, ...                % Positive stall speed vectors
                                                                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_VS.value, ...                % Negative stall speed vectors
                                                                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value, ...              % Max dive speed from V - n diagram
                                                                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VE.value, ...              % VE = VD on the negative side of V - n diagram
                                                                                       Aircraft.Certification.Regulation.value, ...                                                 % Chosen certification 
                                                                                       Aircraft.Certification.Aircraft_Name.value);                                                 % Selected aircraft name

% Saving figures inside correct folder
dir = pwd;
fprintf("--------------------------------------");
fprintf('\n');
fprintf('### Saving outpus inside correct Folder ###');
fprintf('\n');
SaveFolder = strcat(dir,'\Output');
fprintf('Saving Vndiagram.pdf in: ');
fprintf('\n');      
fprintf('%s\n', SaveFolder);
% SaveFolder
% Moving file inside correct folder
movefile Vndiagram.pdf Output
movefile Vndiagram.png Output
% -------------------------------------------------------------------------
%% GUST ENVELOPE 
% Vectors with airspeed values
% Two vectors with airspeed values within following ranges:
% 1st ---> [0, VC] where VC = MaxCruiseSpeed
% 2nd ---> [0, VD] where VD = MaxDiveSpeed
% indexes = 500;
indexes = 1e5;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value = linspace(0.0, Aircraft.Certification.Regulation.SubpartC.Flightloads.Cruise_Speed_VC.value, indexes)'; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.Attributes.unit = 'm/s';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_dive.value = linspace(0.0, Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value, indexes)'; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_dive.Attributes.unit = 'm/s';
% -------------------------------------------------------------------------
%% x = calcmug(obj, Wingloading, MAC, NormalForceCurveSlope, g)
% This function calculates the MASS RATIO for the selected airplane
% following the CS-VLA airworthiness prescriptions. To have a
% complete documentation check the class file csvla.m
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Mass_ratio.value = calcmug(obj, Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value, ...   % Wing loading in SI units
                                                                                Aircraft.Geometry.Wing.mac.value, ...                                        % Mean Aerodynamic Chord in meters
                                                                                Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value, ...  % Normal force curve slope (practically equal to lift curve slope
                                                                                Aircraft.Certification.ISA_Condition.rho.value, Aircraft.Constants.g.value); % Gravity acceleration g
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Mass_ratio.Attributes.unit = 'Non dim.';
% -------------------------------------------------------------------------
%% x = calckg(obj, MassRatio)
% This function calculates the GUST ALLEVIATION FACTOR for the
% selected airplane and flight conditions following the CS-VLA
% airworthiness prescriprions. To have a complete documentation
% check the class fil csvla.m
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor.value = calckg(obj, Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Mass_ratio.value);  
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor.Attributes.unit = 'Non dim.'; 
% -------------------------------------------------------------------------
%% x = calcngust(rho, NormalForceCurveSlope, GustAlleviationFact, GustSpeedCruiseVect, WingLoading, CruiseSpeed, DiveSpeed, FlagToCalc)
% This function is able to calculates in any possible case a vector
% which contains gust load factors value, following CS-VLA
% airworthiness prescription. To have a complete documentation
% check the class file csvla.m

Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.value = calcngust(obj, Aircraft.Certification.ISA_Condition.rho.value, ...                                          % Standard atmosphere density
                                                                                            Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value, ...                       % Normal force curve slope [1/rad]
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor.value, ...           % Gust alleviation factor KG
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.value, ...                 % Gust speed vector at flight speed V = VC
                                                                                            Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value, ...                             % Wing loading in SI units
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Cruise_Speed_VC.value, ...                 % Cruise speed from the V - n diagram 
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value, ...                   % Dive speed from the V - n diagram
                                                                                            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.Attributes.case(1))); % A conveniently defined case switch
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.Attributes.unit = 'g';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise.value = calcngust(obj, Aircraft.Certification.ISA_Condition.rho.value, ...                                           % Standard atmosphere density
                                                                                            Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value, ...                        % Normal force curve slope [1/rad]
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor.value, ...            % Gust alleviation factor KG
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.value, ...                  % Gust speed vector at flight speed V = VC
                                                                                            Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value, ...                              % Wing loading in SI units
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Cruise_Speed_VC.value, ...                  % Cruise speed from the V - n diagram
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value, ...                    % Dive speed from the V - n diagram 
                                                                                            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.Attributes.case(2)));  % A conveniently defined case switch
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise.Attributes.unit = 'g';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive.value = calcngust(obj, Aircraft.Certification.ISA_Condition.rho.value, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value, ...
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor.value, ...
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.value, ... 
                                                                                            Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value, ...
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Cruise_Speed_VC.value, ...
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value, ...
                                                                                            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.Attributes.case(1)));
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive.Attributes.unit = 'g';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive.value = calcngust(obj, Aircraft.Certification.ISA_Condition.rho.value, ...
                                                                                            Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value, ...
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor.value, ...
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.value, ... 
                                                                                            Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value, ...
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Cruise_Speed_VC.value, ...
                                                                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value, ...
                                                                                            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.Attributes.case(2)));
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive.Attributes.unit = 'g';     
% -------------------------------------------------------------------------

%% INPUT TRACKING - VN DIAGRAM 
% A possible way to track inputs for the various data will be provided
% inside the .txt file used as a log for the program.

%% OUTPUT
% x = Gust_envelope(npos, nneg, VSpos, VSneg, VD, nmax, nmin, VG, ...
%                   V_cruise, V_dive, ngust_plus_cruise, ngust_neg_cruise, ...
%                   ngust_plus_dive, ngust_neg_dive, Reg, Aircraft_name ) 
% This function construct, plot and save in various format the
% gust envelope of the aircraft, following CS-VLA airworthiness
% prescription. To have a complete documentation check the class
% file csvla.m 

disp(" ")
disp(" ++++ FIGURE 2 - GUST ENVELOPE ++++ "); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_envelope.value = Gust_envelope(obj, Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value, ... % Vector which contains positive load factor values
                                                                                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.value, ...    % Vector which contains negative load factor values 
                                                                                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_VS.value, ...              % Vector which contains positive stall speed values 
                                                                                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_VS.value, ...              % Vector which contains negative stall speed values 
                                                                                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value, ...            % Max dive speed from V - n diagram 
                                                                                           Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value, ...                     % Max positive load factor 
                                                                                           Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value, ...                     % Min (negative) load factor
                                                                                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VE.value, ...            % VG = VD on the V - n diagram 
                                                                                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value, ...            % Airspeed resulting from gust when flight speed is V = VC
                                                                                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_dive.value, ...              % Airspeed resulting from gust when flight speed is V = VD
                                                                                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.value, ...       % Positive load factors associated with wind gust, V = VC 
                                                                                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise.value, ...       % Negative load factors associated with wind gust, V = VC
                                                                                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive.value, ...         % Positive load factors associated with wind gust, V = VD
                                                                                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive.value, ...         % Negative load factors associated with wind gust, V = VD
                                                                                           Aircraft.Certification.Regulation.value, ...                                               % Airworthiness rules applied 
                                                                                           Aircraft.Certification.Aircraft_Name.value);                                               % Selected aircraft name

% Saving figures inside correct folder
fprintf('Saving Gustenvelope.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile Gustenvelope.pdf Output
movefile Gustenvelope.png Output 
% -----------------------------------------------------------------
%% FINAL ENVELOPE 
% Testing functions for the final envelope. In this part of the code, it is
% necessary to combine Flight Envelope and Gust envelope airspeeds and load
% factors to obtain the Final Flight Envelope. In the latter diagram, we
% can clearly understand the effect on structural design and sizing of wind
% gusts, following airworthiness rules prescription (in absence of other
% rational criteria or methodologies).
% [new_vstall new_nstall] = stall_speed_limit1(VA, vstall, vgust, nstall, ngust, nmax)
%  This function is able to obtain a stall speed vector which
%  include the gust speed lines tracked with respect the cruise
%  speed of the aircraft. This function can be used for positive
%  and negative stall speed and load factors values. 
tol = 1e-2;
twod_vect = stall_speed_limit1(obj, Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.value, ... % Design manoeuvring speed VA
                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_VS.value, ...            % Positive stall speed vector
                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value, ...          % Airspeed from gust envelope, V = VC
                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_load_factors.value, ...  % Positive load factor vector
                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.value, ...     % Positive load factor resulting from gust, V = VC
                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value, ...
                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Cruise_Speed_VC.value);                      % Max positive load factor
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value = twod_vect(:, 1); % Store output from calculations inside the struct variable AIRCRAFT
index = 1;
y = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value(end);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value)
    x = abs(y - Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value(i));
    if x < tol
        index = i;
    end
end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value; Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value(index:end)];
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.Attributes.unit = 'm/s';                                                                           
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value = twod_vect(:, 2);  % Store output from calculations inside the struct variable AIRCRAFT
y = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value(end);
index_load_factor = 1;
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.value)
    x = abs(y - Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.value(i));
    if x < tol
        index_load_factor = i;
    end
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value; Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value(index_load_factor:end)];
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.Attributes.unit = 'm/s';    
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value; Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.value(index_load_factor:end)];
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.Attributes.unit = 'g';
clear index y;
% Testing functions for the final envelope
twod_vect = stall_speed_limit1(obj, Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_Design_manoeuvring_speed_VG.value, ...
                                                                                                              Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_VS.value, ...
                                                                                                              Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value, ...
                                                                                                              Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_load_factors.value, ...
                                                                                                              Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise.value, ...
                                                                                                              Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value, ...
                                                                                                              Aircraft.Certification.Regulation.SubpartC.Flightloads.Cruise_Speed_VC.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value = twod_vect(:, 1);
index = 1;
y = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value(end);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value)
    x = abs(y - Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value(i));
    if x < tol
        index = i;
    end
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value; Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value(index:end)];
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.Attributes.unit = 'm/s';                                                                           
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.value = twod_vect(:, 2); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.value = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.value; Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise.value(index:end)]; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.Attributes.unit = 'g';
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value; Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value(index_load_factor:end)];
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.Attributes.unit = 'm/s';                                                                           

clear index;
% -------------------------------------------------------------------------
% [vAD; nAD] =  back_envelope_points(obj, vgust, ngust, VD, vstall, nstall, nmax)    
%  This function is able to obtain a stall speed vector which
%  include the gust speed lines tracked with respect the cruise
%  speed of the aircraft. This function can be used for positive
%  and negative stall speed and load factors values. Notice that the output
%  is not a single vector. 
point_env = back_envelope_points(obj, Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_dive.value, ... 
                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive.value, ...
                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value, ...
                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value, ...
                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value, ...
                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value);
% The following lines store all the calculated values inside the struct variable AIRCRAFT
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.value = point_env(1, 1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.Attributes.unit = 'm/s';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.value = point_env(2, 1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.Attributes.unit = 'g';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value = point_env(2, 2);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.Attributes.unit = 'g';
point_env = back_envelope_points(obj, Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_dive.value, ... 
                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive.value, ...
                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value, ...
                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value, ...
                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.value, ...
                                            Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.value = point_env(1, 1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.Attributes.unit = 'm/s';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.value = point_env(2, 1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.Attributes.unit = 'g';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value = point_env(2, 2);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.Attributes.unit = 'g';
% CHECK ON LOAD FACTORS AND SPEED
Check1 = abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value - Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value);
if Check1 < tol
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value;
end
%
V_check2 = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value(end) Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_dive.value(end)]';
n_check2 = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.value(end) Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive.value(end)]';
V_interp = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_cruise.value(end), Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_dive.value(end), indexes)';
N_straight_line = interp1(V_check2, n_check2, V_interp);
for i = 1:length(N_straight_line)
    Check2 = abs(N_straight_line(i) - Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value);
    if Check2 < tol
        V_final_gust = V_interp(i);
        N_final_gust = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value;
    end
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_speed.value = V_final_gust;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_speed.Attributes.unit = "m/s";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_load_factor.value = N_final_gust;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_load_factor.Attributes.unit = "m/s";


%% INPUT TO THE FINAL ENVELOPE DIAGRAM 

disp(" ")
disp(" ++++ INPUT FINAL ENVELOPE ++++");
% Horizontal tail loads increments
Data4 = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.value, ...
         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.value, ...
         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_speed.value, ...
         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_load_factor.value, ...
         Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value, ...
         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value, ...
         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.value, ...
         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.value];            % Negative load factors associated with wind gust, V = VD                 % VG = VD on the negative side of V - n diagram
disp(" ++++++++++ DATA USED TO PLOT FINAL ENVELOPE ++++++++++ ")
format = ' %6.6f    %6.6f     %6.6f    %6.6f    %6.6f    %6.6f    %6.6f    %6.6f\n';
label  = ' VC           nC           V_fg         n_fg        VD           nD          VF            nF\n';
fprintf(label);
fprintf(format, Data4.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

disp(" ")
% Horizontal tail loads increments
Data5 = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value, ...
         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value, ...
         Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.value, ...
         Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value, ...
         Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_Design_manoeuvring_speed_VG.value, ...
         Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value];            % Negative load factors associated with wind gust, V = VD                 % VG = VD on the negative side of V - n diagram
disp(" ++++++++++ DATA USED TO PLOT FINAL ENVELOPE ++++++++++ ")
format = ' %6.6f    %6.6f     %6.6f    %6.6f    %6.6f    %6.6f\n';
label  = ' VE            nE           VA           nA          VG            nG\n';
fprintf(label);
fprintf(format, Data5.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

disp(" ")
disp(" ++++ FIGURE 3 - FINAL ENVELOPE PLOT ++++ ");
% -------------------------------------------------------------------------
%  fig1 = V_n_diagram(npos, nneg, VSpos, VSneg, VD, VG, VA, VE, Reg, Aircraft_name)
%  This function plot the V - n diagram, based on the applied regulation.
%  The applied regulation is stored inside the local variable 'Reg' for
%  convenience and it is used to automatically change the output figure
%  title name. Also, the selected aircraft name is stored inside the
%  variable 'Aircraft_name' and is inserted in the diagram as plain text.
%  This might be a useful feature. 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Diagram.value = Final_envelope(obj, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value, ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.value, ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value, ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value, ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VC.value, ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.value, ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_speed.value, ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Final_gust_load_factor.value, ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value, ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value, ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_speed_VF.value, ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nF.value, ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Dive_Speed_VD.value, ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value, ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.value, ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value, ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_Design_manoeuvring_speed_VG.value, ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value, ...
                                                                                         Aircraft.Certification.Regulation.value, ... 
                                                                                         Aircraft.Certification.Aircraft_Name.value);


% Saving figures inside correct folder
fprintf('Saving Finalenvelope.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile Finalenvelope.pdf Output
movefile Finalenvelope.png Output
% -------------------------------------------------------------------------
% Aerodynamic Data useful to calculate balancing loads. It is useful to
% store data inside the struct variable AIRCRAFT. The following lines are 
% necessary to transform the imported string variable in a usable numerical
% vector. 
Aircraft.Certification.Aerodynamic_data.alpha.value = str2num(Aircraft.Certification.Aerodynamic_data.alpha.value); % A vector which contains AoA values 
Aircraft.Certification.Aerodynamic_data.CL.value    = str2num(Aircraft.Certification.Aerodynamic_data.CL.value);    % A vector which contains CL values 
Aircraft.Certification.Aerodynamic_data.CD.value    = str2num(Aircraft.Certification.Aerodynamic_data.CD.value);    % A vector which contains CD values 
Aircraft.Certification.Aerodynamic_data.CM.value    = str2num(Aircraft.Certification.Aerodynamic_data.CM.value);    % A vector which contains CM values 
% A figure with polars subplots is automatically generated from the vectors

disp(" ")
disp(" ++++ FIGURE 4 - AERODYNAMIC DATA ++++ ");
Aircraft.Certification.Aerodynamic_data.Polars = AeroPlot(Aircraft.Certification.Aerodynamic_data.alpha.value, ...
                                                          Aircraft.Certification.Aerodynamic_data.CL.value, ...
                                                          Aircraft.Certification.Aerodynamic_data.CD.value, ...
                                                          Aircraft.Certification.Aerodynamic_data.CM.value);

% Saving figures inside correct folder
fprintf('Saving Polars.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile Polars.pdf Output
movefile Polars.png Output
% -------------------------------------------------------------------------
%% Evaluation of the balancing loads 
% BalancingLoads