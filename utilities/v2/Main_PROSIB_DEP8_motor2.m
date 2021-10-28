%% main file to run

clc; clearvars; close all;

diary logfile.txt

%% aircraft name here

ac.name = 'PROSIB_DEP8_motor2_cruise_propON';      % aircraft name
propON = 1;     % 1 if rotors are to be enabled into VSPaero, 0 otherwise

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
ac.fuse.length   = 22.67;               % length
ac.fuse.diameter = 2.87;                % diameter
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

ac.wing.id          = 'wing';                  % component name
ac.wing.type        = 'WING';                  % OpenVSP component type
ac.wing.airfoil     = 'XS_FIVE_DIGIT';         % OpenVSP airfoil standard type
% By choosing XS_FIVE_DIGITS, OpenVSP automatically generates a symmetrical
% 5-digits NACA airfoil. It's not the only option the software gives you:
% both type and shape of the airfoil may be changed and a custom one
% imported.
ac.wing.idealcl     = [0.6 0.6 0.6];           % camber of each airfoil
ac.wing.camberloc   = [0.15 0.15 0.15];        % camber location of each airfoil
ac.wing.thickchord  = [0.18 0.18 0.12];        % thickness of each airfoil
ac.wing.twist       = [0 0 -2];                % twist angles
ac.wing.twistloc    = [0 0 0];                 % twist locations along the chord
ac.wing.xloc        = 8.77;                    % X global coordinate
ac.wing.yloc        = 0.0;                    % y global coordinate
ac.wing.zloc        = 1.50;                    % Z global coordinate
ac.wing.xrot        = 0;                       % X global rotation
ac.wing.yrot        = 0;                       % Y global rotation
% -180° <= y.rot <= 180°; yrot = 180° means right wing becomes left wing
ac.wing.zrot        = 0;                       % Z global rotation

ac.wing.sym         = 2;                       % wing planar symmetry
% 0: none; 1: XY; 2: XZ; 3: YZ
ac.wing.span        = [4.73 7.56];             % span(s)
ac.wing.ctip        = [2.14 1.09];             % tip chord(s)
ac.wing.croot       = [2.14 ac.wing.ctip(1)];  % root chord(s)
%movables
ac.wing.mov.number = 1.0;               %number of movables on wing
%flap
ac.wing.mov.type = 'f';                 %elevator type
ac.wing.mov.eta_inner = 0.08;
ac.wing.mov.eta_outer = 0.725;
ac.wing.mov.cf_c_inner = 0.35;
ac.wing.mov.cf_c_outer = 0.33;
ac.wing.mov.tessellation = 10;
ac.wing.mov.rotation = 0;           % flap deflection

% ac.wing.mov.type = 'a';                 %aileron type
% ac.wing.mov.eta_inner = 0.725;
% ac.wing.mov.eta_outer = 0.980;
% ac.wing.mov.cf_c_inner = 0.33;
% ac.wing.mov.cf_c_outer = 0.33;
% ac.wing.mov.tessellation = 10;


ac.wing.sweep       = [0 5.3];                % sweep angle(s)
ac.wing.sweeploc    = [0 0];                   % sweep location(s) along the chord
ac.wing.secsweeploc = [1 1];                   % secondary sweep location(s) along the chord
ac.wing.dihedral    = [0 0];                 % dihedral angle(s)
ac.wing.utess       = [40 40];                 % tessellation
ac.wing.wtess       = 60;

for i = 1:numel(ac.wing.span)                     % panels surface(s)
    ac.wing.S(i)    = (ac.wing.croot(i)+...
        ac.wing.ctip(i))*ac.wing.span(i)/2*cosd(ac.wing.dihedral(i));
end

%% Mac calculator
[ac.wing.mac,ac.wing.macxloc] = FindMac(ac.wing);

% Vectors referring to wing panels must be the same size.
% Vectors referring to wing airfoils must be the same size.
% Latter ones are one element bigger than first ones.


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
ac.ver.xloc        = 17.96;
ac.ver.yloc        = 0.0;                    % y global coordinate

ac.ver.zloc        = 0.9*ac.fuse.diameter/2;
ac.ver.xrot        = 90;
ac.ver.yrot        = 0;
ac.ver.zrot        = 0;
ac.ver.sym         = 0;
ac.ver.span        = 4.55;
ac.ver.ctip        = 2.13;
ac.ver.croot       = 3.38;
ac.ver.sweep       = 28.75;
ac.ver.sweeploc    = 0;
ac.ver.secsweeploc = 1;
ac.ver.dihedral    = 0;
ac.ver.utess       = 5;
ac.ver.wtess       = 25;
%movables
%rudder
ac.ver.mov.type = 'r';          %elevator type
ac.ver.mov.eta_inner = 0.0;
ac.ver.mov.eta_outer = 0.9;
ac.ver.mov.cf_c_inner = 0.25;
ac.ver.mov.cf_c_outer = 0.35;
ac.ver.mov.tessellation = 10;

%% horizontal tail

% The same as wing section.

ac.hor.id          = 'horizontal';
ac.hor.type        = 'WING';
ac.hor.airfoil     = 'XS_FOUR_SERIES';
ac.hor.camber      = [0 0];
ac.hor.camberloc   = [0.2 0.2];
ac.hor.thickchord  = [0.10 0.10];
ac.hor.twist       = [0 0];
ac.hor.twistloc    = [0.25 0.25];
ac.hor.xloc        = 20.33;
ac.hor.yloc        = 0.0;                    % y global coordinate

ac.hor.zloc        = ac.ver.zloc + ac.ver.span;
ac.hor.xrot        = 0;
ac.hor.yrot        = 0;
ac.hor.zrot        = 0;
ac.hor.sym         = 2;
ac.hor.span        = 3.55;
ac.hor.ctip        = 1.10;
ac.hor.croot       = 1.88;
ac.hor.sweep       = 12.2;
ac.hor.sweeploc    = 0;
ac.hor.secsweeploc = 1;
ac.hor.dihedral    = 0;
ac.hor.utess       = 15;
ac.hor.wtess       = 25;
%movables
%elevator
ac.hor.mov.type = 'e';          %elevator type
ac.hor.mov.eta_inner = 0.0;
ac.hor.mov.eta_outer = 1.0;
ac.hor.mov.cf_c_inner = 0.50;
ac.hor.mov.cf_c_outer = 0.50;
ac.hor.mov.tessellation = 10;
ac.hor.mov.rotation = 0;           % elevator deflection

%% pods

% Here is actually only given Y location of the pod on the - from a rear
% point of view - right wing. Its symmetrical is obtained by the very
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

numDEP = 8;         % total number of distributed propulsors
etaDEPin = 0.207;     % DEP inner station
etaDEPout = 0.750;    % DEP outer station
ylocDEP = linspace(etaDEPin,etaDEPout,numDEP/2) .* sum(ac.wing.span);
ac.pod.id        = 'pod';   % component name
ac.pod.type      = 'POD';   % OpenVSP component type
ac.pod.utess     = 5;       % u tessellation
ac.pod.wtess     = 5;       % u tessellation
for ii = 1:numDEP/2
    ac.pod.length(ii)    = 0.452;                    % length
    ac.pod.fineratio(ii) = 1.65;                     % fineratio
    ac.pod.xloc(ii)      = 8.685;                  % X local coordinate
    ac.pod.yloc(ii)      = ylocDEP(ii);            % Y local coordinate
    ac.pod.zloc(ii)      = 1.50;                   % Z local coordinate
    ac.pod.yrot(ii)      = 0;                      % Y global rotation
end
% TIP propeller
ac.pod.length(ii+1)    = 3.16;                   % length
ac.pod.fineratio(ii+1) = 3.46;                    % fineratio
ac.pod.xloc(ii+1)      = 8.685;                   % X local coordinate
ac.pod.yloc(ii+1)      = sum(ac.wing.span);     % Y local coordinate
ac.pod.zloc(ii+1)      = 1.50;                  % Z local coordinate
ac.pod.yrot(ii+1)      = 0;                     % Y global rotation

%% propellers

% The same as pods section;

% ac.prop.id       = 'propeller';            % component name
% ac.prop.type     = 'Disk';                 % OpenVSP component type
% ac.prop.xloc     = ac.pod.xloc;            % X global coordinate
% ac.prop.yloc     = ac.pod.yloc;            % Y global coordinate
% ac.prop.zloc     = ac.pod.zloc;            % Z global coordinate
% ac.prop.yrot     = -2;                     % Y global rotation
% ac.prop.diameter = 3.93;                   % diameter
% ac.prop.utess    = 15;                     % tessellation
% ac.prop.wtess    = 15;

ac.prop.id       = 'propeller';         % component name
ac.prop.type     = 'Disk';              % OpenVSP component type
ac.prop.utess    = 15;                  % u tessellation
ac.prop.wtess    = 15;                  % w tessellation
for ii = 1:numDEP/2
    ac.prop.xloc(ii)     = ac.pod.xloc(ii); % X global coordinate
    ac.prop.yloc(ii)     = ac.pod.yloc(ii);     % Y global coordinate
    ac.prop.zloc(ii)     = ac.pod.zloc(ii); % Z global coordinate
    ac.prop.yrot(ii)     = ac.pod.yrot(ii);  % Y global rotation
    ac.prop.diameter(ii) = 2.20;             % diameter
end
% TIP propeller
ac.prop.xloc(ii+1)     = ac.pod.xloc(ii+1);   % X global coordinate
ac.prop.yloc(ii+1)     = ac.pod.yloc(ii+1);       % Y global coordinate
ac.prop.zloc(ii+1)     = ac.pod.zloc(ii+1);   % Z global coordinate
ac.prop.yrot(ii+1)     = ac.pod.yrot(ii+1);    % Y global rotation
ac.prop.diameter(ii+1) = 3.93;                 % diameter

%% configuration cases
ac.config = 'WBHVP';                  % configuration to be analyzed
%
% ac    : aircraft
% ac.config: configuration case
%           W       = wing
%           WP      = wing propeller
%           WB      = wing body
%           WH      = wing horizontal
%           WV      = wing vertical
%           WHP      = wing horizontal propeller
%           WHF      = wing horizontal fan
%           WVP      = wing vertical propeller
%           WVF      = wing vertical fan
%           WBH     = wing body horizontal
%           WBHP     = wing body horizontal propeller
%           WBHF     = wing body horizontal fan
%           WBVP     = wing body vertical propeller
%           WBVF     = wing body vertical fan
%           WBV     = wing body vertical
%           WBHV    = wing body horizontal vertical
%           WBHVP   = wing body horizontal vertical propeller
%           WBHVF   = wing body horizontal vertical fan

ac.movables = 'fe';
% ac.movables: movables configuration (leave empty array for none)
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
str = '### creating vspscript...';
disp(' ');
disp(str);

MakeVspscript(ac);

disp(' ');
str = '### ... vspscript created';
disp(' ');
disp(str);

%% running vspscript.exe

disp(' ');
disp(' ');
str = '### degenerating geometry and calculating mass properties...';
disp(' ');
disp(str);

vspscript_exe = ['vspscript.exe -script ' ac.name '.vspscript'];
dos(vspscript_exe)

disp(' ');
str = '### ... degenerate geometry created and mass properties calculated';
disp(' ');
disp(str);

%% vspaero setup

set.h                       = 20000*.3048;          % altitude in meters
[set.T, set.a, set.p, set.rho] = ...                % freestream thermodynamic
    atmosisa(set.h);                                % characteristics
set.mu                      = ...                   % freestream dynamic viscosity
    (1.4580e-06*set.T^(1.5))/(set.T+110.4);
set.Mach                    = 0.40;                   % freestream Mach number
set.Vinf                    = set.Mach*set.a;       % freestream speed
set.Cref                    = ac.wing.mac;            % reference chord
set.ReCref                  = ...                   % freestream Re based
    (set.rho*set.Vinf*set.Cref)/set.mu;             % on Cref
set.Sref                    = 2*sum(ac.wing.S);       % reference surface
set.Bref                    = 2*sum(ac.wing.span);    % reference span
set.X_cg                    = ac.wing.macxloc+...     % aircraft center of gravity (PDV 9.64)
    0.0*ac.wing.mac;                                 % global coordinates
set.Y_cg                    = 0;
set.Z_cg                    = 0.25*ac.fuse.diameter;

set.AoA                     = -2:2:6;                % attack angle(s)
set.Beta                    = 0;                      % sideslip angle(s)

set.ClMax                   = -1;                     % max sectional lift coefficient
% -1: no limit
set.MaxTurningAngle         = -1;                     % ???

switch set.Beta                                       % imposing symmetry only if
    case 0                                            % set.Beta = 0
        set.Symmetry        = 'Y';                    %
    otherwise                                         % plane of symmetry:
        set.Symmetry        = 'NO';                   % NO: no; X: YZ; Y: XZ; Z: XY
end

set.FarDist                 = -1;                     % distance where wakes return to freestream:
% -1: VSPAERO computes the distance
set.NumWakeNodes            = -1;                     % number of wake nodes
% -1: free
set.WakeIters               = 3;                      % number of iterations

if propON == 1                                        % enable all rotors
    set.NumberOfRotors = numel(ac.prop.yloc);         % set geometric and aerodyanmic data arrays
    set.xloc = ac.prop.xloc;                          % prop global X location
    set.yloc = ac.prop.yloc;                          % prop global Y location
    set.zloc = ac.prop.zloc;                          % prop global Z location
    set.xdir = cosd(ac.prop.yrot);                    % prop X unit normal
    set.ydir = zeros(1,set.NumberOfRotors);           % prop Y unit normal
    set.zdir = sind(ac.prop.yrot);                    % prop Z unit normal
    set.dia = ac.prop.diameter;                       % prop diameter
    set.hub = 0.1.*set.dia;                           % prop hub
    
    % set RPM to zero to disable that rotor into VSPaero
    set.rpm = [0 * ones(1,set.NumberOfRotors-1), ... % DEP props RPM
        -900];                                          % TIP props RPM
    set.ct = [0.15 * ones(1,set.NumberOfRotors-1), ...   % DEP props CT
        0.25];                                           % TIP props CT
    set.cp = [0.134 * ones(1,set.NumberOfRotors-1), ...   % DEP props CP
        0.548];                                          % TIP props CP
    
else
    set.NumberOfRotors          = 0;                  % prop-off condition
end

set.movables            = ac.movables;                % enable analysis with movables
set.wingID = ac.wing.id;
set.horID = ac.hor.id;
set.flapRot             = ac.wing.mov.rotation;
set.elevRot             = ac.hor.mov.rotation;

set.VRTXLift                = 'Y';                    % ???
% Y: yes; N: no
set.LESuction               = 'Y';                    % leading edge suction
% Y: yes; N: no
set.KT                      = 'Y';                    % 2nd order Karman-Tsien Mach correction
% Y: yes; N: no

%% making vspaero based on given setup

disp(' ');
disp(' ');
str = '### creating vspaero...';
disp(' ');
disp(str);

MakeVspaero(ac.name,set);

disp(' ');
str = '### ... vspaero created';
disp(' ');
disp(str);

%% running vspaero.exe

disp(' ');
disp(' ');
str = '### running solver...';
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

%% importing and plotting results and then moving all i/o files of both vspscript.exe
% and vspaero.exe in results directory

res = ImportResults(ac);
save([ac.name,'_res.mat'],'res')

PlotResults(set,res,ac.name,ac.config);

MoveResults(ac.name,res_dir);
movefile([ac.name,'_res.mat'],res_dir)

%% Stability analysis about alfa = 0° and beta = 0°
if ~exist('stab','dir')
    mkdir('stab')
end
copyfile([res_dir,'\',ac.name '_DegenGeom.csv'],pwd)

commandwindow % focus on the command window, leaving the plots behind

set.AoA = 0;
set.Beta = 0;
set.Symmetry = 'NO'; % symmetry must be set off on stability calculation

disp(' ');
disp(' ');
str = '### creating vspaero for stability analysis about alfa=0 and beta=0...';
disp(' ');
disp(str);

MakeVspaero(ac.name,set);

disp(' ');
str = '### ... vspaero created';
disp(' ');
disp(str);

disp(' ');
disp(' ');
str = '### running solver for stability analysis...';
disp(' ');
disp(str);
disp(' ');
disp(' ');

vspaero_exe     = ['vspaero -omp ' num2str(num_cpu)...   % command
    ' -stab ' ac.name '_DegenGeom'];
try
    ex = dos(vspaero_exe);                                % executing command
    if ex ~= 0
        warning('VSPaero failed to perform stability calculations')
    end
catch
    warning('VSPaero failed to perform stability calculations')
end

movefile([ac.name '_DegenGeom.*'],'stab');
movefile('stab',res_dir)

%% moving remaining file in results directory
diary off
movefile('logfile.txt',res_dir);
