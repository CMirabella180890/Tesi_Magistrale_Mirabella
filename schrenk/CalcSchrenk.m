
%% CALCULATION OF THE LOAD DISTRIBUTION WITH SCHRENK'S METHOD

%% Number of elements
numb = 39; 

%% Aircraft data 
b      = Aircraft.Geometry.Wing.b.value; 
S      = Aircraft.Geometry.Wing.S.value; 
c_root = Aircraft.Geometry.Wing.croot.value;
c_tip  = Aircraft.Geometry.Wing.ctip.value;

%% Solution of the problem 
Results         = Schrenk_load_distr(b, S, c_root, c_tip, numb);
y               = Results(:, 1);
eta             = Results(:, 2); 
chord           = Results(:, 3);
elliptical_load = Results(:, 4); 
Schrenk_cCl     = Results(:, 5);
Unit_Cl         = Results(:, 6);

%% Plot results 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% FIGURE 7 - SPANWISE LIFT DISTRIBUCTION
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
disp(" ")
disp(" ++++ FIGURE 7 - SCHRENK SPANWISE LIFT DISTRIBUTION ++++ ");
figure1 = Schrenk_Cl_diagram(eta, Unit_Cl);

% EXPORT FIGURE
exportgraphics(figure1, 'Unit_Cl.pdf', 'ContentType', 'vector')
exportgraphics(figure1, 'Unit_Cl.png', 'ContentType', 'vector')

% Saving figures inside correct folder
fprintf('Saving Unit_Cl.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile Unit_Cl.pdf Output
movefile Unit_Cl.png Output 
% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% FIGURE 8 - SCHRENK'S LOADS 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
disp(" ")
disp(" ++++ FIGURE 8 - SCHRENK SPANWISE LIFT DISTRIBUTION ++++ ");
figure2 = Schrenk_plots(eta, chord, elliptical_load, Schrenk_cCl);

% EXPORT FIGURE
exportgraphics(figure1, 'SchrenkLoads.pdf', 'ContentType', 'vector')
exportgraphics(figure1, 'SchrenkLoads.png', 'ContentType', 'vector')

% Saving figures inside correct folder
fprintf('Saving SchrenkLoads.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile SchrenkLoads.pdf Output
movefile SchrenkLoads.png Output 
% -------------------------------------------------------------------------

% STORE INSIDE THE AIRCRAFT STRUCT VARIABLE
% A vector of spanwise stations
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.y.value = y; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.y.Attributes.unit = "m";
% Non dimensional spanwise stations
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.eta.value = eta; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.eta.Attributes.unit = "Non dimensional";
% Chord value at each spanwise station
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord.value = chord; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord.Attributes.unit = "m";
% Elliptical load
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.elliptical_load.value = elliptical_load; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.elliptical_load.Attributes.unit = "m";
% Schrenk's load distribution
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Schrenk_cCl.value = Schrenk_cCl; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Schrenk_cCl.Attributes.unit = "m";
% Unit Cl 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Unit_Cl.value = Unit_Cl; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Unit_Cl.Attributes.unit = "m";