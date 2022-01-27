%% main file to run
% 
% clc; clear all; close all;

diary logfile.txt

%% aircraft name here

ac.name = Aircraft.Certification.Aircraft_Name.value;
                              % aircraft name

str = ['### aircraft name is ' ac.name];
disp('##############################################');
disp('##############################################');
disp(' ');
disp(str);

%% Checking results directory

main_dir = pwd;                                 % current directory
res_dir  = [main_dir '\' ac.name '_results'];   % results directory

if exist(res_dir,'dir')
    rmdir(res_dir,'s')                          % removing directory
end

mkdir(res_dir)                                  % making directory

%% aircraft data

%% fuselage

ac.fuse.id       = 'fuselage';          % component name
ac.fuse.type     = 'TransportFuse';     % OpenVSP component type
ac.fuse.length = Aircraft.Geometry.Fuselage.length.value;
% ac.fuse.length   = 3.64;               % length
ac.fuse.diameter = Aircraft.Geometry.Fuselage.diameter.value;
% ac.fuse.diameter = 0.42;                % diameter
ac.fuse.utess    = 30;                  % tessellation
ac.fuse.wtess    = 15;

																		% There are no other (in the form of blank components) reference system 
% There are no other (in the form of blank components) reference system 
% than OpenVSP global one, thus, if given, position and orientation specs
% may only refers to the latter one. Since such info are not given, OpenVSP
% automatically points the nose at the origin and lays the body on the
% positive X global axis.

%% wing

% By choosing WING type, OpenVSP automatically makes a one-panel tapered
% wing. It's up to the user to add panels and shape each one of them.
% The software actually creates only - from a rear point of view - the
% right wing and then automatically the left one as its symmetrical.
% Winglets are abtained as very small panels with a big dihedral.

ac.wing.id          = 'wing';                  % component name
ac.wing.type        = 'WING';                  % OpenVSP component type


% By choosing XS_FIVE_DIGITS, OpenVSP automatically generates a symmetrical
% 5-digirs NACA airfoil. It's not the only option the software gives you:
% both type and shape of the airfoil may be changed and a custom one
% imported.

ac.wing.airfoil     = 'XS_FILE_AIRFOIL';         % OpenVSP airfoil standard type

% Define the airfoil for each section in a cell array and point to the airfoil folder
% ac.wing.affile = {'Nasa_Langley_GAW1_0417.dat', 'Nasa_Langley_GAW1_0417.dat'};
ac.wing.affile = {'e205_Selig.dat', 'e205_Selig.dat', 'e205_Selig.dat', 'e205_Selig.dat'};
%ac.wing.affile = {'E216_104.dat', 'E216_104.dat', 'E216_104.dat', 'E216_104.dat'};

ac.wing.affile = strcat('airfoil/',ac.wing.affile); 

ac.wing.idealcl     = [0.6 0.6 0.6 0.6];           % camber of each airfoil             [0.6 0.6 0.6 0.6]
ac.wing.camberloc   = [0.15 0.15 0.15 0.15];        % camber location of each airfoil     [0.15 0.15 0.15 0.15]
ac.wing.thickchord  = [0.18 0.18 0.18 0.18];        % thickness of each airfoil           [0.18 0.18 0.18 0.18]
ac.wing.twist       = [3.0 3.0 3.0 3.0];               % twist angles                   [3.0 3.0 3.0 3.0]
ac.wing.twistloc    = [0 0 0 0];                 % twist locations along the chord  [0 0 0 0]
ac.wing.xloc        = Aircraft.Geometry.Wing.xle.value;                    % X global coordinate     
ac.wing.yloc        = 0.0;
ac.wing.zloc        = Aircraft.Geometry.Wing.zle.value;                    % Z global coordinate    
ac.wing.xrot        = 0;                       % X global rotation
ac.wing.yrot        = 0;                       % Y global rotation
% -180° <= y.rot <= 180°; yrot = 180° means right wing becomes left wing
ac.wing.zrot        = 0;                         % Z global rotation
ac.wing.sym         = 2;                       % wing planar symmetry
% 0: none; 1: XY; 2: XZ; 3: YZ
% to be MODIFIED!!!
panel_span = Aircraft.Geometry.Wing.b.value/6;
%
ac.wing.span        = [panel_span panel_span panel_span];     % span(s)
ac.wing.ctip        = [Aircraft.Geometry.Wing.ctip.value Aircraft.Geometry.Wing.ctip.value Aircraft.Geometry.Wing.ctip.value];  % tip chord(s)[k1,k2,tip]
ac.wing.croot       = [Aircraft.Geometry.Wing.croot.value Aircraft.Geometry.Wing.croot.value Aircraft.Geometry.Wing.croot.value];  % root chord(s)[root,k1,k2]
% ac.wing.croot       = [ac.wing.ctip(1)];  % root chord(s)
ac.wing.sweep       = [0 0 0];                % sweep angle(s) [25 25 25]
ac.wing.sweeploc    = [0 0 0];                   % sweep location(s) along the chord
ac.wing.secsweeploc = [1 1 1];                   % secondary sweep location(s) along the chord
ac.wing.dihedral    = [0 0 0];                 % dihedral angle(s) [2 2 2]
ac.wing.utess       = [15 15 12];                 % tessellation [15 15 12]
ac.wing.wtess       = 40;
ac.wing.captess     = 9;    % Number of tessellated curves on Wing Root and Tip
ac.wing.tecluster   = 0.15;               %TE Tess Cluster Control
ac.wing.lecluster   = 0.15;               %LE Tess Cluster Control

for i = 1:numel(ac.wing.span)                     % panels surface(s)
    ac.wing.S(i)    = (ac.wing.croot(i)+...
        ac.wing.ctip(i))*ac.wing.span(i)/2*cosd(ac.wing.dihedral(i));
end

% mod PDV 10-2021
%movables
ac.wing.mov.number = 2.0;               %number of movables on wing
% If ac.wing.mov.number = 1.0;              
% only aileron or flap to be specified -> NO BOTH
% Example flap + aileron 
ac.wing.mov.type = ['f' ; 'a'];                 %1st element flap - 2nd aileron
ac.wing.mov.eta_inner = [0.10; 0.70];          % 
ac.wing.mov.eta_outer = [0.70; 1.0];
ac.wing.mov.cf_c_inner = [0.25; 0.3];
ac.wing.mov.cf_c_outer = [0.3; 0.3];
ac.wing.mov.tessellation = [10; 10];
% ac.wing.mov.type = 'a';                 %aileron type
% ac.wing.mov.eta_inner = 0.70;
% ac.wing.mov.eta_outer = 1.0;
% ac.wing.mov.cf_c_inner = 0.3;
% ac.wing.mov.cf_c_outer = 0.3;
% ac.wing.mov.tessellation = 10;

%% Mac calculator
[ac.wing.mac,ac.wing.macxloc] = FindMac(ac.wing);
															  %further check for rectangula wing with multiple panles PDV-22
if isfield(Aircraft.Geometry.Wing,'type') == 1
    if Aircraft.Geometry.Wing.type.value == 'Rectangular'
        ac.wing.mac = ac.wing.croot(1);
        ac.wing.macxloc = ac.wing.xloc;       % Rectangular or...
    end
end
% Vectors referring to wing panels must be the same size.
% Vectors referring to wing airfoils must be the same size.
% Latter ones are one element bigger than first ones.

%% horizontal tail
% The same as wing section.

ac.hor.id          = 'horizontal';
ac.hor.type        = 'WING';
ac.hor.airfoil     = 'XS_FOUR_SERIES';
ac.hor.camber      = [0 0];
ac.hor.camberloc   = [0.2 0.2];
ac.hor.thickchord  = [0.12 0.12];
ac.hor.twist       = [0 0];
ac.hor.twistloc    = [0.25 0.25];
ac.hor.xloc        = 1.49+1.638;
ac.hor.yloc        = 0.0;
ac.hor.zloc        = 0.15;
ac.hor.xrot        = 0;
ac.hor.yrot        = 0;
ac.hor.zrot        = 0;
ac.hor.sym         = 2;
ac.hor.span        = 0.748;
ac.hor.ctip        = 0.3136;
ac.hor.croot       = 0.3929;
ac.hor.sweep       = 15;
ac.hor.sweeploc    = 0;
ac.hor.secsweeploc = 1;
ac.hor.dihedral    = 0;
ac.hor.utess       = [15];
ac.hor.wtess       = 25;
%movables
%elevator
ac.hor.mov.type = 'e';          %elevator type
ac.hor.mov.eta_inner = 0.0;
ac.hor.mov.eta_outer = 1.0;
ac.hor.mov.cf_c_inner = 0.3;
ac.hor.mov.cf_c_outer = 0.3;
ac.hor.mov.tessellation = 10;

%% vertical tail

% The same as wing section.
% Vertical tail is obtained as a 90° rotated right wing.

ac.ver.id          = 'vertical';
ac.ver.type        = 'WING';
ac.ver.airfoil     = 'XS_FOUR_SERIES';
ac.ver.camber      = [0 0 0 0];
ac.ver.camberloc   = [0.2 0.2 ];
ac.ver.thickchord  = [0.1 0.1 ];
ac.ver.twist       = [0 0 ];
ac.ver.twistloc    = [0 0 ];
ac.ver.xloc        = ac.hor.xloc+ac.hor.span*tan(ac.hor.sweep/57.3);  %Vtail at the Htail tip
ac.ver.yloc        = ac.hor.span;                     % Y global coordinate
ac.ver.zloc        = 0.15;
ac.ver.xrot        = 90;
ac.ver.yrot        = 0;
ac.ver.zrot        = 0;
ac.ver.sym         = 2;
ac.ver.span        = 0.438;
ac.ver.ctip        = 0.1534725;
ac.ver.croot       = 0.3136;
ac.ver.sweep       = 20;
ac.ver.sweeploc    = 0;
ac.ver.secsweeploc = 1;
ac.ver.dihedral    = 0;
ac.ver.utess       = [15];
ac.ver.wtess       = 25;%movables
%rudder
ac.ver.mov.type = 'r';          %elevator type
ac.ver.mov.eta_inner = 0.0;
ac.ver.mov.eta_outer = 1.0;
ac.ver.mov.cf_c_inner = 0.25;
ac.ver.mov.cf_c_outer = 0.3;
ac.ver.mov.tessellation = 10;

%% pods

% Here is actually only given Y location of the pod on the - from a rear
% point of view - right wing. Its symmetricals is obtained by the very
% function making the engine.

% ac.pod.id        = 'pod';                  % component name
% ac.pod.type      = 'POD';                  % OpenVSP component type
% ac.pod.length    = 4.4;                    % length
% ac.pod.fineratio = 10;                     % fineratio
% ac.pod.xloc      = 8.14;                   % X local coordinate
% ac.pod.yloc      = 4.05;                   % Y local coordinate
% ac.pod.zloc      = 1.42;                   % Z local coordinate
% ac.pod.yrot      = 0;                      % Y global rotation
% ac.pod.utess     = 5;                      % tessellation
% ac.pod.wtess     = 5;

%% propellers

%
% 
% ac.prop.id       = 'propeller';            % component name
% ac.prop.type     = 'Disk';                 % OpenVSP component type
% ac.prop.xloc     = [-0.1 -0.1 -0.1 -0.1];            % X global coordinate
% ac.prop.yloc     = [0.4 0.7 1.0 1.4];            % Y global coordinate
% ac.prop.zloc     = [0 0 0 0];            % Z global coordinate
% ac.prop.yrot     = [0 0 0 0];                     % Y global rotation
% ac.prop.diameter = [0.2 0.2 0.2 0.4];                   % diameter
% ac.prop.utess    = 15;                     % tessellation
% ac.prop.wtess    = 15;

%% configuration cases
ac.config = 'WB';                  % configuration to be analyzed
% W    : wing
% WB   : wing-body
% WH   : wing-horizontal
% WBH  : wing-body-horizontal
% WBHV : wing-body-horizontal-vertical
% WBHVP: wing-body-horizontal-vertical-turboprop
% WBHVF: wing-body-horizontal-vertical-turbofan

ac.movables = 'none';
% ac.movables: movables configuration
%           none       = none movables
%           f       = flap
%           a       = aileron
%           fa      = flap aileron
%           fe      = flap elevator
%           fr      = flap rudder
%           fae     = flap aileron elevator
%           faer    = flap aileron elevator rudder
%           e       = elevator
%           r       = rudder

disp(' ');
disp(' ');
str = ['### aircraft configurations to be analyzed is ' ac.config];
disp(' ');
disp(str);

disp(' ');
disp(' ');
str = ['### aircraft movables to be analyzed is ' ac.movables];
disp(' ');
disp(str);

%% making vspscript

disp(' ');
disp(' ');
str = ['### creating vspscript...'];
disp(' ');
disp(str);

MakeVspscript(ac);

disp(' ');
str = ['### ... vspscript created'];
disp(' ');

%% running vspscript.exe

disp(' ');
disp(' ');
str = ['### degenerating geometry and calculating mass properties...'];
disp(' ');
disp(str);

vspscript_exe = ['vspscript.exe -script ' ac.name '.vspscript'];
dos(vspscript_exe)

disp(' ');
str = ['### ... degenerate geometry created and mass properties calculated'];
disp(' ');

%% vspaero setup

set.h                       = 0;                      % altitude
 [set.T set.a set.p set.rho] = ...                    % freestream thermodynamic
     atmosisa(set.h);                                 % characteristics
set.mu                      = ...                     % freestream dynamic viscosity
     (1.4580e-06*set.T^(1.5))/(set.T+110.4);
% set.Mach                    = 0.1;                  % freestream Mach number
set.Vinf = 19.0;
%set.Vinf                    = double(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value(1));       % freestream speed
set.Mach                    = set.Vinf/set.a;         % freestream Mach number
% set.rho                     = 1.225;                % freestream density
% set.Vinf                    = 34.0;                 % freestream speed
set.Cref                    = ac.wing.mac;            % reference chord
% set.ReCref                  = ...                   % freestream Re based
%     (set.rho*set.Vinf*set.Cref)/set.mu;             % on Cref
set.ReCref                  = 1e6;                    % freestream Re based
% on Cref
set.Sref                    = 2*sum(ac.wing.S);       % reference surface
set.Bref                    = 2*sum(ac.wing.span);    % reference span
set.X_cg                    = ac.wing.macxloc+...     % aircraft center of gravity (PDV 9.64)
    0.25*ac.wing.mac;                                 % global coordinates
set.Y_cg                    = 0;
set.Z_cg                    = 0;
set.AoA                    = 0:2:12;                 % attack angle(s)
% +++++++++++++++++++++++++++++++++++
% ANGLE OF ATTACK CORRESPONDING TO CL = 1
% set.AoA                     = 3.1862;
% +++++++++++++++++++++++++++++++++++
set.Beta                    = 0;                      % sideslip angle(s)
set.ClMax                   = -1;                     % max sectional lift coefficient
% -1: no limit
set.MaxTurningAngle         = -1;                     % ???

switch set.Beta                                       % imposing symmetry only if
    case 0                                            % set.Beta = 0
        set.Symmetry        = 'Y';                    %
    otherwise                                         % plane of symmetry:
        set.Symmetry        = 'N';                    % N: no; X: YZ; Y: XZ; Z: XY
end

set.FarDist                 = -1;                     % distance where wakes return to freestream:
% -1: VSPAERO computes the distamce
set.NumWakeNodes            = -1;                     % number of wake nodes
% -1: free
set.WakeIters               = 3;                      % number of iterations
set.NumberOfRotors          = 0;                      % number of actuator disks
set.VRTXLift                = 'Y';                    % ???
% Y: yes; N: no
set.LESuction               = 'Y';                    % leading edge suction
% Y: yes; N: no
set.KT                      = 'Y';                    % 2nd order Karman-Tsien Mach correction
% Y: yes; N: no

%% making vspaero based on given setup

disp(' ');
disp(' ');
str = ['### creating vspaero...'];
disp(' ');
disp(str);

MakeVspaero(ac.name,set);

disp(' ');
str = ['### ... vspaero created'];
disp(' ');

%% running vspaero.exe

disp(' ');
disp(' ');
str = ['### running solver...'];
disp(' ');
disp(str);
disp(' ');
disp(' ');

num_cpu         = 4;                                     % cpus number
vspaero_exe     = ['vspaero -omp ' num2str(num_cpu)...   % command
    ' ' ac.name '_DegenGeom'];
try
    ex = dos(vspaero_exe);                                % executing command
    if ex ~= 0
        warning('VSPaero failed to execute')
    end
catch
    warning('VSPaero failed to execute')
end

%% importing nd plotting results and then moving all i/o files of both vspscript.exe
% and vspaero.exe in results directory

res = ImportResults(ac);

PlotResults(set,res,ac.name,ac.config);

MoveResults(ac.name,res_dir);

diary off

%% moving remaining file in results directory

movefile('logfile.txt',res_dir);
