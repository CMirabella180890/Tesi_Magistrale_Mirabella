%% Script to evaluate shear and bending moment distr. along the main wing span

%   DESCRIPTION
%    In this script, the lift and drag distribution along the span will be 
%    used to evaluate shear and bending moment along the main wing span.
%    First, it is necessary to store the critical point, which can be
%    exctracted from the V-N diagram. 

%% LOAD A CLASS OF FUNCTIONS USEFUL TO EVALUATE ALL THE REQUIRED DATA
obj2 = auxiliaries; 

% Test speed
Aircraft.Certification.Regulation.SubpartC.Balancingloads.Test_speed.value = 34.0; 
Aircraft.Certification.Regulation.SubpartC.Balancingloads.Test_speed.Attributes.unit = "m/s";

%% DEFINE A CHORD DISTRIBUTION 
Aircraft.Geometry.Wing.y.value = linspace(0, ...
                                          Aircraft.Geometry.Wing.b.value, ...
                                          length(Aircraft.Certification.Aerodynamic_data.Cl_alongspanCL1.value));
Aircraft.Geometry.Wing.y.Attributes.unit = "m";
                                     
% Calculation of a chord distribution with a convenient, simple function
Aircraft.Geometry.Wing.taper.value = Aircraft.Geometry.Wing.ctip.value/Aircraft.Geometry.Wing.croot.value;
Aircraft.Geometry.Wing.taper.Attributes.unit = "m";

Aircraft.Certification.Regulation.SubpartC.Balancingloads.chord_distr.value = calc_chord(obj2, Aircraft.Geometry.Wing.S.value, ...
                                                                                         Aircraft.Geometry.Wing.taper.value, ...
                                                                                         Aircraft.Geometry.Wing.b.value, ...
                                                                                         Aircraft.Geometry.Wing.y.value);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.chord_distr.Attributes.unit = "m";

%% PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
Aircraft.Certification.Regulation.SubpartC.Balancingloads.cCl_distr.value = Aircraft.Certification.Aerodynamic_data.Cl_alongspanCL1.value*Aircraft.Certification.Regulation.SubpartC.Balancingloads.chord_distr.value;
Aircraft.Certification.Regulation.SubpartC.Balancingloads.cCl_distr.Attributes.unit = "m";

Aircraft.Certification.Regulation.SubpartC.Balancingloads.cCd_distr.value = Aircraft.Certification.Aerodynamic_data.Cd_alongspanCL1.value*Aircraft.Certification.Regulation.SubpartC.Balancingloads.chord_distr.value;
Aircraft.Certification.Regulation.SubpartC.Balancingloads.cCd_distr.Attributes.unit = "m";

%% ALPHA VECTOR ALONG THE MAIN WING SPAN 
cd .. 
cd utilities
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

% Calling the function alpha_calc
Aircraft.Certification.Regulation.SubpartC.Balancingloads.AoA_alongspan_rad.value = alpha_calc(obj1, ...
                                        Aircraft.Certification.Aerodynamic_data.Cl_alongspanCL1.value, ...
                                        Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                        Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                        Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value, ...
                                        Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a.value, ...
                                        Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b.value, ...
                                        Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c.value);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.AoA_alongspan_rad.Attributes.unit = "Radians"; 

% Convert to radians 
Aircraft.Certification.Regulation.SubpartC.Balancingloads.AoA_alongspan_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.Balancingloads.AoA_alongspan_rad.value);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.AoA_alongspan_deg.Attributes.unit = "Degrees";

% Changing folder 
cd .. 
cd csvla
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

%% TOTAL ANGLE OF ATTACK DATA STORAGE

% AoA_Tot = AoA + Twist_angle of the main wing
Aircraft.Certification.Regulation.SubpartC.Balancingloads.AoA_Tot_deg.value = Aircraft.Geometry.Wing.twist_angle.value + Aircraft.Certification.Regulation.SubpartC.Balancingloads.AoA_alongspan_deg.value;
Aircraft.Certification.Regulation.SubpartC.Balancingloads.AoA_Tot_deg.Attributes.unit = "Degrees"; 

% Convert in radians 
Aircraft.Certification.Regulation.SubpartC.Balancingloads.AoA_Tot_rad.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.Balancingloads.AoA_Tot_deg.value);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.AoA_Tot_rad.Attributes.unit = "Radians";

%% PROJECTION OF FORCES ALONG WING AXES - NORMAL AND AXIAL FORCES 

% Calculation of the normal force coefficient
Aircraft.Certification.Regulation.SubpartC.Balancingloads.cCz.value = calc_normal_force(obj2, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Balancingloads.AoA_Tot_rad.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Balancingloads.cCl_distr.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Balancingloads.cCd_distr.value);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.cCz.Attributes.unit = "m"; 

% Calculation of the axial force coefficient 
Aircraft.Certification.Regulation.SubpartC.Balancingloads.cCa.value = calc_axial_force(obj2, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Balancingloads.AoA_Tot_rad.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Balancingloads.cCl_distr.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Balancingloads.cCd_distr.value);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.cCa.Attributes.unit = "m"; 

% Normal force 
Aircraft.Certification.Regulation.SubpartC.Balancingloads.Normal_force.value = Aircraft.Certification.Regulation.SubpartC.Balancingloads.cCz.value*(0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Certification.Regulation.SubpartC.Balancingloads.Test_speed.value^2);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.Normal_force.Attributes.unit = "N/m";

% Axial force 
Aircraft.Certification.Regulation.SubpartC.Balancingloads.Axial_force.value = Aircraft.Certification.Regulation.SubpartC.Balancingloads.cCa.value*(0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Certification.Regulation.SubpartC.Balancingloads.Test_speed.value^2);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.Axial_force.Attributes.unit = "N/m";

%% SHEAR FORCE CALCULATION 
Aircraft.Certification.Regulation.SubpartC.Balancingloads.Shear_distr.value = calc_shear_force(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                               Aircraft.Certification.Regulation.SubpartC.Balancingloads.Normal_force.value);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.Shear_distr.Attributes.unit = "N/m";

prova_figura_Taglio = figure;
plot(flip(Aircraft.Geometry.Wing.y.value'), Aircraft.Certification.Regulation.SubpartC.Balancingloads.Shear_distr.value);
%% BENDING MOMENT CALCULATION 
Aircraft.Certification.Regulation.SubpartC.Balancingloads.Bend_mom_distr.value = calc_bend_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                               Aircraft.Certification.Regulation.SubpartC.Balancingloads.Shear_distr.value);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.Bend_mom_distr.Attributes.unit = "N";

%% PLANS FOR THE STRUCTURAL DIMENSIONING
% Tu correctly size the aerostructures of the main lifting surface 
% it is necessary to apply the procedure just developed to the 
% critical points coming from the V-N diagram. Those point represents 
% the most demanding flight conditions that our aircraft could survive. 
% Those points are all stored inside: 
% --> Aircraft.Certification.Regulation.SubpartC.Final_Envelope
% Retrieve the values and apply formulas to them. 