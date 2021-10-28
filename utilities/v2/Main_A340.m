%% main file to run

clc; clear all; close all;

diary logfile.txt

%% aircraft name here

ac.name = 'A340';                               % aircraft name

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
ac.fuse.type     = 'TransportFuse';     % OpenVSP coponent type
ac.fuse.length   = 73.84;               % length
ac.fuse.diameter = 5.64;                % diameter
ac.fuse.utess    = 30;                  % tessellation
ac.fuse.wtess    = 15;

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

ac.wing.id          = 'wing';                     % component name
ac.wing.type        = 'WING';                     % OpenVSP component type
ac.wing.airfoil     = 'XS_FOUR_SERIES';           % OpenVSP airfoil standard type
% By choosing XS_FOUR_SERIES, OpenVSP automatically generates a symmetrical
% 5-digirs NACA airfoil. It's not the only option the software gives you:
% both type and shape of the airfoil may be changed and a custom one
% imported.
ac.wing.camber      = [0.02 0.02 0.02 0.02];      % camber of each airfoil
ac.wing.camberloc   = [0.4 0.4 0.4 0.4];          % camber location of each airfoil
ac.wing.thickchord  = [0.12 0.12 0.12 0.12];      % thickness of each airfoil
ac.wing.twist       = [3 1 -3 0];                  % twist angles
ac.wing.twistloc    = [0 0 0 0];                  % twist locations along the chord
ac.wing.xloc        = 26.19;                      % X global coordinate
ac.wing.yloc        = 0.0;                        % Y global coordinate
ac.wing.zloc        = -2.15;                          % Z global coordinate
ac.wing.xrot        = 0;                          % X global rotation
ac.wing.yrot        = 0;                          % Y global rotation
% -180° <= y.rot <= 180°; yrot = 180° means right wing becomes left wing
ac.wing.sym         = 2;                          % wing planar symmetry
                                                  % 0: none; 1: XY; 2: XZ; 3: YZ
ac.wing.span        = [9.33 16.28 1.39];          % span(s)
ac.wing.ctip        = [9.35 2.74 0.51];           % tip chord(s)
ac.wing.croot       = [14.50 ac.wing.ctip(1)...
    ac.wing.ctip(2)];  % root chord(s)
%flap
ac.wing.mov.type = 'f';                 %elevator type
ac.wing.mov.eta_inner = 0.10;
ac.wing.mov.eta_outer = 0.70;
ac.wing.mov.cf_c_inner = 0.25;
ac.wing.mov.cf_c_outer = 0.3;
ac.wing.mov.tessellation = 10;

% ac.wing.mov.type = 'a';                 %aileron type
% ac.wing.mov.eta_inner = 0.70;
% ac.wing.mov.eta_outer = 1.0;
% ac.wing.mov.cf_c_inner = 0.3;
% ac.wing.mov.cf_c_outer = 0.3;
% ac.wing.mov.tessellation = 10;

for i = 1:numel(ac.wing.span)                     % panels surface(s)
    ac.wing.S(i)    = (ac.wing.croot(i)+...
        ac.wing.ctip(i))*ac.wing.span(i)/2;
end
ac.wing.sweep       = [32 32 60];                 % sweep angle(s)
ac.wing.sweeploc    = [0 0 0];                    % sweep location(s) along the chord 
ac.wing.secsweeploc = [1 1 1];                    % secondary sweep location(s) along the chord
ac.wing.dihedral    = [6.2 6.2 50];               % dihedral angle(s)
% [ac.wing.mac,ac.wing.macxloc] = FindMac(ac.wing);
ac.wing.utess       = [15 15 4];                         % tessellation
ac.wing.wtess       = 25;

[ac.wing.mac,ac.wing.macxloc] = FindMac(ac.wing);

% Vectors referring to wing panels must be the same size.
% Vectors referring to wing airfoils must be the same size.
% Latter ones are one element bigger than first ones.

%% horizontal tail

% The same as wing section.

ac.hor.id          = 'horizontal';
ac.hor.type        = 'WING';
ac.hor.airfoil     = 'XS_FOUR_SERIES';
ac.hor.camber      = [0.02 0.02];
ac.hor.camberloc   = [0.4 0.4];
ac.hor.thickchord  = [0.12 0.12];
ac.hor.twist       = [0 0];
ac.hor.twistloc    = [0 0];
ac.hor.xloc        = 64.67;
ac.hor.yloc        = 0.0;                        % Y global coordinate
ac.hor.zloc        = 1.50;
ac.hor.xrot        = 0;
ac.hor.yrot        = 0;
ac.hor.sym         = 2;
ac.hor.span        = 11.54;
ac.hor.ctip        = 2.43;
ac.hor.croot       = 5.85;
ac.hor.sweep       = 30;
ac.hor.sweeploc    = 0; 
ac.hor.secsweeploc = 1;
ac.hor.dihedral    = 5.85;
ac.hor.utess       = [15];
ac.hor.wtess       = 25;
%movables
%elevator
ac.hor.mov.type = 'e';          %elevator type
ac.hor.mov.eta_inner = 0.0;
ac.hor.mov.eta_outer = 1.0;
ac.hor.mov.cf_c_inner = 0.25;
ac.hor.mov.cf_c_outer = 0.3;
ac.hor.mov.tessellation = 10;
%% vertical tail

% The same as wing section.
% Vertical tail is obtained as a 90° rotated right wing.

ac.ver.id          = 'vertical';
ac.ver.type        = 'WING';
ac.ver.airfoil     = 'XS_FOUR_SERIES';
ac.ver.camber      = [0 0];
ac.ver.camberloc   = [0.2 0.2];
ac.ver.thickchord  = [0.1 0.1];
ac.ver.twist       = [0 0];
ac.ver.twistloc    = [0 0];
ac.ver.xloc        = 60;
ac.ver.zloc        = 0;
ac.ver.yloc        = 0.0;                        % Y global coordinate
ac.ver.xrot        = 90;
ac.ver.yrot        = 0;
ac.ver.sym         = 0;
ac.ver.span        = [11.05];
ac.ver.ctip        = [3.42];
ac.ver.croot       = 10.49;
ac.ver.sweep       = 45;
ac.ver.sweeploc    = 0; 
ac.ver.secsweeploc = 1;
ac.ver.dihedral    = 0;
ac.ver.utess       = 15;
ac.ver.wtess       = 25;
%movables
%rudder
ac.ver.mov.type = 'r';          %elevator type
ac.ver.mov.eta_inner = 0.0;
ac.ver.mov.eta_outer = 1.0;
ac.ver.mov.cf_c_inner = 0.25;
ac.ver.mov.cf_c_outer = 0.3;
ac.ver.mov.tessellation = 10;

%% nacelles

% Here are actually only given Y locations of the nacelles on the - from a
% rear point of view - right wing. Their symmetricals are obtained by the
% very function making the nacelles.

ac.nacelle.id         = 'nacelle';          % component name
ac.nacelle.type       = 'Duct';             % OpenVSP component type
ac.nacelle.width      = 2.62;               % width
ac.nacelle.height     = 2.62;               % height
ac.nacelle.chord      = 7.35;               % chord
ac.nacelle.thickchord = 0.05;               % thickchord
ac.nacelle.xloc       = [27.32 33.95];      % X global coordinate(s)
ac.nacelle.yloc       = [9.53 19.85];       % Y global coordinate(s)
ac.nacelle.zloc       = [-2 -1.5];            % Z global coordinate(s)
ac.nacelle.yrot       = [0 0];              % Y global rotation(s)    
ac.nacelle.utess      = 5;                 % tessellation
ac.nacelle.wtess      = 5;

%% configuration cases

ac.config = 'WBHV';                  % configuration to be analyzed
% W    : wing
% WB   : wing-body
% WH   : wing-horizontal
% WBH  : wing-body-horizontal
% WBHV : wing-body-horizontal-vertical
% WBHVP: wing-body-horizontal-vertical-prop
% WBHVF: wing-body-horizontal-vertical-trubofan

ac.movables = 'faer';
% ac.movables: movables configuration
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
str = ['### degenerating geometry and caltulating mass properties...'];
disp(' ');
disp(str);

vspscript_exe = ['vspscript.exe -script ' ac.name '.vspscript'];
dos(vspscript_exe)

disp(' ');
str = ['### ... degenerate geometry created and mass properties calculated'];
disp(' ');

%% vspaero setup

set.h                       = 1e5;                    % altitude
[set.T set.a set.p set.rho] = ...                     % freestream thermodynamic 
    atmosisa(set.h);                                  % characteristics
set.mu                      = ...                     % freestream dynamic viscosity
    (1.4580e-06*set.T^(1.5))/(set.T+110.4);
set.Mach                    = 0.78;                   % freestream Mach number
set.Vinf                    = set.Mach*set.a;
set.Cref                    = ac.wing.mac;            % reference chord
set.ReCref                  = ...                     % freestream Re based
    (set.rho*set.Vinf*set.Cref)/set.mu;               % on Cref
set.Sref                    = 2*sum(ac.wing.S);       % reference surface
set.Bref                    = 2*sum(ac.wing.span);    % reference span
set.X_cg                    = ac.wing.macxloc +...    % aircraft center of gravity
    0.25*ac.wing.mac;                                 % global coordinates
set.Y_cg                    = 0;                      % global coordinates
set.Z_cg                    = 0;
set.AoA                     = 0:2:10;                 % attack angle(s)
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
set.LESuction               = 'N';                    % leading edge suction
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

%% importing results and then moving all i/o files of both vspscript.exe
 % and vspaero.exe in results directory

res = ImportResults(ac);

PlotResults(set,res,ac.name,ac.config);

MoveResults(ac.name,res_dir);

diary off

%% moving remaining file in results directory

movefile('logfile.txt',res_dir);
