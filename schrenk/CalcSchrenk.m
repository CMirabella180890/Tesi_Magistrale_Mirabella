
%% CALCULATION OF THE LOAD DISTRIBUTION WITH SCHRENK'S METHOD

% Number of elements
numb = 1e3; 

% Aircraft data 
b      = Aircraft.Geometry.Wing.b.value; 
S      = Aircraft.Geometry.Wing.S.value; 
c_root = Aircraft.Geometry.Wing.croot.value;
c_tip  = Aircraft.Geometry.Wing.ctip.value;

% Solution of the problem 
Results         = Schrenk_load_distr(b, S, c_root, c_tip, numb);
y               = Results(:, 1);
eta             = Results(:, 2); 
chord           = Results(:, 3);
elliptical_load = Results(:, 4); 
Schrenk_cCl     = Results(:, 5);
Unit_Cl         = Results(:, 6);

% Plot results 
figure1 = Schrenk_Cl_diagram(eta, Unit_Cl);
figure2 = Schrenk_plots(eta, chord, elliptical_load, Schrenk_cCl);

% STORE INSIDE THE AIRCRAFT STRUCT VARIABLE