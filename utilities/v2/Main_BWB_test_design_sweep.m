%% main file to run
clc; clearvars; close all;

%% Creating diary...
diary logfile.txt
    
     
%% Combinations       
% Valori delle singole variabili assunte nel design sweep
AR_inner = [0.5, 2];
AR_outer = [3.5, 5.5];
taper_inner = [0.2, 0.35];
taper_outer = [0.15, 0.5];
C_root = 40;  % C_tip_i = C_root_i * taper_i
sweep_inner = [61, 64];
sweep_outer = [37, 41];
thick_root = 0.18;
thick_kink = [0.10, 0.12];
thick_tip = [0.08, 0.10];
twist_kink = [1, 2];
twist_tip = [-4, -5];


%% Aircraft name...
dir_name = 'BWB_test';

%% Design sweep...
% Create progress bar
progBar = waitbar(0,'Running loop...','Name','Progress bar',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(progBar,'canceling',0);
set(progBar.Children(1), 'String', 'Abort')

% Pre-allocate arrays
steps = length(AR_inner) * length(AR_outer) * length(taper_inner) * ...
    length(taper_outer) * length(C_root) * length(sweep_inner) * ...
    length(sweep_outer) * length(thick_root) * length(thick_kink) * ...
    length(thick_tip) * length(twist_kink) * length(twist_tip);
stepTimeHistory = NaN(steps,1);


file_index = 1;
for i_twist_tip=1:length(twist_tip)
    for i_twist_kink=1:length(twist_kink)
        for i_thick_tip=1:length(thick_tip)
            for i_thick_kink=1:length(thick_kink)
                for i_sweep_outer=1:length(sweep_outer)
                    for i_sweep_inner=1:length(sweep_inner)
                        for i_taper_outer=1:length(taper_outer) 
                            for i_taper_inner=1:length(taper_inner) 
                                for i_AR_outer=1:length(AR_outer)
                                    for i_AR_inner=1:length(AR_inner)
            
                                tic;
                                        
                                % Check for clicked canceling button
                                if getappdata(progBar,'canceling')
                                    waitbar(file_index/steps,progBar,sprintf('Aborting...'))
                                    break
                                end
        
                                %% Variable assignment...
                                c_root_inner = C_root; % root velivolo
                                c_tip_inner = c_root_inner * taper_inner(i_taper_inner); % tip primo pannello = root secondo pannello
                                c_root_outer = c_tip_inner;
                                c_tip_outer = c_root_outer * taper_outer(i_taper_outer);
                                b_inner = AR_inner(i_AR_inner) * c_root_inner * (1 + taper_inner(i_taper_inner))/2;
                                b_outer = AR_outer(i_AR_outer) * c_root_outer * (1 + taper_outer(i_taper_outer))/2;
            
                                %% Enumerating aircraft...
                                Code = "";
                                if file_index < 1000
                                  Code = Code + "0";
                                end
                                if file_index < 100
                                  Code = Code + "0";
                                end
                                if file_index < 10
                                  Code = Code + "0";
                                end
                                Code = Code + file_index;
        
                                ac.name = strcat(dir_name,'_',num2str(file_index));
                                file_index = file_index + 1;
                                disp('##############################################');
                                disp('##############################################');
                                disp(strcat("Writing file ", ac.name, "..."));
                                propON = 0;     % 1 if rotors are to be enabled into VSPaero, 0 otherwise

                                str = ['### aircraft name is ' ac.name];
                                disp('##############################################');
                                disp('##############################################');
                                disp(' ');
                                disp(str);

        %% Checking results directory
        main_dir = [pwd '\' dir_name '_results'];                                 % current directory 
        res_dir = [strcat(main_dir,'\', strcat(dir_name,'_',Code))];
       
       
        if exist(res_dir,'dir')
            rmdir(res_dir,'s')                          % removing directory
        end
        mkdir(res_dir)                                  % making directory

        
        %% Writing legend...
        Table = table(  Code, AR_inner(i_AR_inner), AR_outer(i_AR_outer), taper_inner(i_taper_inner), taper_outer(i_taper_outer), sweep_inner(i_sweep_inner), sweep_outer(i_sweep_outer), thick_kink(i_thick_kink), thick_tip(i_thick_tip), twist_kink(i_twist_kink), twist_tip(i_twist_tip), ...
                        'VariableNames', {'Aircraft #', 'AR_inner', 'AR_outer', 'taper_inner', 'taper_outer', 'sweep_inner', 'sweep_outer', 'thick_kink', 'thick_tip', 'twist_kink', 'twist_tip'});
        if file_index == 2
            writetable(Table, 'Legend.xlsx', 'Range',strcat("A",num2str(file_index-1)));
        else
            writetable(Table, 'Legend.xlsx', 'Range',strcat("A",num2str(file_index)),'WriteVariableNames', false);
        end
       
       
              
        %% aircraft data

        % By choosing WING type, OpenVSP automatically makes a one-panel tapered
        % wing. It's up to the user to add panels and shape each one of them.
        % The software actually creates only - from a rear point of view - the
        % right wing and then automatically the left one as its symmetrical.
        % Winglets are abtained as very small panels with a big dihedral.

        ac.wing.id          = 'wing';                  % component name
        ac.wing.type        = 'WING';                  % OpenVSP component type

        % By choosing XS_FIVE_DIGITS, OpenVSP automatically generates a symmetrical
        % 5-digits NACA airfoil. It's not the only option the software gives you:
        % both type and shape of the airfoil may be changed and a custom one
        % imported.

        % ac.wing.idealcl     = [0.4 0.4 0.4];           % camber of each airfoil
        % ac.wing.camber      = [0.02 0.02 0.02];      % camber of each airfoil

        ac.wing.airfoil     = 'XS_FILE_AIRFOIL';         % OpenVSP airfoil standard type

        % Define the airfoil for each section in a cell array and point to the airfoil folder
        ac.wing.affile = {'Liebeck_BWB_airfoil.dat',...
            'NASA-SC-0612.txt', 'NASA-SC-0612.txt'};
        ac.wing.affile = strcat('airfoil/',ac.wing.affile); 

        ac.wing.thickchord  = [thick_root, thick_kink(i_thick_kink), thick_tip(i_thick_tip)];      % thickness of each airfoil
        ac.wing.twist       = [0, twist_kink(i_twist_kink), twist_tip(i_twist_tip)];                  % twist angles
        ac.wing.twistloc    = [0.25, 0.25, 0.25];                  % twist locations along the chord
        ac.wing.xloc        = 0.0;                        % X global coordinate (usually from Nose of the fuselage)
        ac.wing.zloc        = 0;                          % Z global coordinate
        ac.wing.xrot        = 0;                          % X global rotation
        ac.wing.yrot        = 0;                          % Y global rotation
        % -180° <= y.rot <= 180°; yrot = 180° means right wing becomes left wing
        ac.wing.zrot        = 0;                         % Z global rotation
        ac.wing.sym         = 2;                          % wing planar symmetry
                                                          % 0: none; 1: XY; 2: XZ; 3: YZ
        ac.wing.span        = [b_inner, b_outer];          % span(s)
        ac.wing.ctip        = [c_tip_inner, c_tip_outer];           % tip chord(s)
        ac.wing.croot       = [c_root_inner, c_root_outer];  % root chord(s)
        ac.wing.sweep       = [sweep_inner(i_sweep_inner), sweep_outer(i_sweep_outer)];  % sweep angle(s)
        ac.wing.sweeploc    = [0 0];                    % sweep location(s) along the chord 
        ac.wing.secsweeploc = [1 1];                    % secondary sweep location(s) along the chord
        ac.wing.dihedral    = [6.0 6.0];               % dihedral angle(deg.)
        % [ac.wing.mac,ac.wing.macxloc] = FindMac(ac.wing);
        ac.wing.utess       = [15 15];                         % tessellation
        ac.wing.wtess       = 25;

        % BLENDING PANELS
        % XxYyMatch is an enum used to identify the type of wing blending between XSecs
        %    BLEND_FREE = 0, 
        %    BLEND_ANGLES = 1, 
        %    BLEND_MATCH_IN_LE_TRAP = 2, 
        %    BLEND_MATCH_IN_TE_TRAP = 3, 
        %    BLEND_MATCH_OUT_LE_TRAP = 4, 
        %    BLEND_MATCH_OUT_TE_TRAP = 5, 
        %    BLEND_MATCH_IN_ANGLES = 6, 
        %    BLEND_MATCH_LE_ANGLES = 7, 
        %    BLEND_NUM_TYPES = 8 
        % where
        % BLEND_FREE 	
        % Free blending
        % 
        % BLEND_ANGLES 	
        % Blend based on angles (sweep & dihedral)
        % 
        % BLEND_MATCH_IN_LE_TRAP 	
        % Match inboard leading edge trapezoid
        % 
        % BLEND_MATCH_IN_TE_TRAP 	
        % Match inboard trailing edge trapezoid
        % 
        % BLEND_MATCH_OUT_LE_TRAP 	
        % Match outboard leading edge trapezoid
        % 
        % BLEND_MATCH_OUT_TE_TRAP 	
        % Match outboard trailing edge trapezoid
        % 
        % BLEND_MATCH_IN_ANGLES 	
        % Match inboard angles
        % 
        % BLEND_MATCH_LE_ANGLES 	
        % Match leading edge angles
        % 
        % BLEND_NUM_TYPES 	
        % Number of blending types

        ac.wing.blending = 1; % blending flag
        % the following arrays size must match section data size

        ac.wing.blend.InLEMode = [NaN, 4, 0]; % leading edge inboard blending type
        ac.wing.blend.InLESweep = [NaN, ac.wing.sweep]; % leading edge inboard sweep angle
        ac.wing.blend.InLEDihedral = [NaN, ac.wing.dihedral]; % leading edge inboard sweep angle
        ac.wing.blend.InLEStrength = [NaN, 0.4, 0.4]; % leading edge inboard sweep angle

        ac.wing.blend.OutLEMode = [1, 0, NaN];  % leading edge outboard blending type
        ac.wing.blend.OutLESweep = [0, ac.wing.sweep]; % leading edge inboard sweep angle
        ac.wing.blend.OutLEDihedral = [0, ac.wing.dihedral]; % leading edge inboard sweep angle
        ac.wing.blend.OutLEStrength = [0.4, 0.4, 0.4]; % leading edge inboard sweep angle

        ac.wing.blend.InTEMode = [NaN, 5, 0]; % trailing edge inboard blending type
        ac.wing.blend.InTESweep = [NaN, ac.wing.sweep]; % leading edge inboard sweep angle
        ac.wing.blend.InTEDihedral = [NaN, ac.wing.dihedral]; % leading edge inboard sweep angle
        ac.wing.blend.InTEStrength = [NaN, 0.4, 0.4]; % leading edge inboard sweep angle

        ac.wing.blend.OutTEMode = [1, 0, NaN];  % trailing edge outboard blending type
        ac.wing.blend.OutTESweep = [0, ac.wing.sweep]; % leading edge inboard sweep angle
        ac.wing.blend.OutTEDihedral = [0, ac.wing.dihedral]; % leading edge inboard sweep angle
        ac.wing.blend.OutTEStrength = [0.4, 0.4, 0.4]; % leading edge inboard sweep angle



        %movables
        ac.wing.mov.number = 0;               % number of movables on wing
        %flap
        % ac.wing.mov.type = 'f';                 % elevon treated as flap type
        % ac.wing.mov.eta_inner = [0.1, 1.01*ac.wing.span(1)/sum(ac.wing.span)];
        % ac.wing.mov.eta_outer = [0.99*ac.wing.span(1), ...
        %     0.99*ac.wing.span(1)+0.50*ac.wing.span(2)] / sum(ac.wing.span);
        % ac.wing.mov.cf_c_inner = [0.1, 0.2];
        % ac.wing.mov.cf_c_outer = [0.2, 0.25];
        % ac.wing.mov.tessellation = [10, 15];
        % ac.wing.mov.rotation = [0, 0];           % movable deflection

        % ac.wing.mov.type = 'a';                 %aileron type
        % ac.wing.mov.eta_inner = 0.70;
        % ac.wing.mov.eta_outer = 1.0;
        % ac.wing.mov.cf_c_inner = 0.3;
        % ac.wing.mov.cf_c_outer = 0.3;
        % ac.wing.mov.tessellation = 10;

        for i = 1:numel(ac.wing.span)                     % panels surface(s)
            ac.wing.S(i)    = (ac.wing.croot(i)+...
                ac.wing.ctip(i))*ac.wing.span(i)/2*cosd(ac.wing.dihedral(i));
        end

        %% Mac calculator
        [ac.wing.mac,ac.wing.macxloc] = FindMac(ac.wing);

        % Vectors referring to wing panels must be the same size.
        % Vectors referring to wing airfoils must be the same size.
        % Latter ones are one element bigger than first ones.


        %% Vertical tail

        % The same as wing section.
        % Vertical tail is obtained as a 90° rotated right wing.

        ac.ver.id          = 'vertical';
        ac.ver.type        = 'WING';
        ac.ver.airfoil     = 'XS_FOUR_SERIES';
        ac.ver.camber      = [0 0 0];
        ac.ver.camberloc   = [0.2 0.2 0.2];
        ac.ver.thickchord  = [0.1 0.1 0.1];
        ac.ver.twist       = [0 0 0];
        ac.ver.twistloc    = [0 0 0];
        ac.ver.xloc        = 20.0;
        ac.ver.zloc        = -4.0;
        ac.ver.xrot        = 90;
        ac.ver.yrot        = 0;
        ac.ver.zrot        = 0; % Probably OpenVSP rotates with a 321 sequence
        ac.ver.sym         = 0;
        ac.ver.span        = [2.5 6.0];
        ac.ver.ctip        = [5.75 2.30];
        ac.ver.croot       = [18.82 ac.ver.ctip(1)];
        ac.ver.sweep       = [80 45];
        ac.ver.sweeploc    = [0 0]; 
        ac.ver.secsweeploc = [1 0];
        ac.ver.dihedral    = [0 0];
        ac.ver.utess       = [12 15];
        ac.ver.wtess       = 25;
        % %movables
        % %rudder
        % ac.ver.mov.type = 'r';          %elevator type
        % ac.ver.mov.eta_inner = 0.0;
        % ac.ver.mov.eta_outer = 1.0;
        % ac.ver.mov.cf_c_inner = 0.25;
        % ac.ver.mov.cf_c_outer = 0.3;
        % ac.ver.mov.tessellation = 10;

        %% nacelles

        % Here are actually only given Y locations of the nacelles on the - from a
        % rear point of view - right wing. Their symmetricals are obtained by the
        % very function making the nacelles.
        % 
        ac.nacelle.id         = 'nacelle';          % component name
        ac.nacelle.type       = 'Duct';             % OpenVSP component type
        ac.nacelle.width      = 2.62;               % width
        ac.nacelle.height     = 2.62;               % height
        ac.nacelle.chord      = 7.35;               % chord
        ac.nacelle.thickchord = 0.05;               % thickchord
        ac.nacelle.xloc       = [28.0];             % X global coordinate(s)
        ac.nacelle.yloc       = [2.0];              % Y global coordinate(s)
        ac.nacelle.zloc       = [2];                % Z global coordinate(s)
        ac.nacelle.yrot       = [0];                % Y global rotation(s)    
        ac.nacelle.utess      = 5;                  % tessellation
        ac.nacelle.wtess      = 5;

        % %% pods
        % 
        % % Here is actually only given Y location of the pod on the - from a rear
        % % point of view - right wing. Its symmetrical is obtained by the very
        % % function making the engine.
        % 
        % % ac.pod.id        = 'pod';                  % component name
        % % ac.pod.type      = 'POD';                  % OpenVSP component type
        % % ac.pod.length    = 4.4;                    % length
        % % ac.pod.fineratio = 10;                     % fineratio
        % % ac.pod.xloc      = 8.14;                   % X local coordinate
        % % ac.pod.yloc      = 4.05;                   % Y local coordinate
        % % ac.pod.zloc      = 1.42;                   % Z local coordinate
        % % ac.pod.yrot      = 0;                      % Y global rotation
        % % ac.pod.utess     = 5;                      % tessellation
        % % ac.pod.wtess     = 5;
        % 
        % numDEP = 8;         % total number of distributed propulsors
        % etaDEPin = 0.207;     % DEP inner station
        % etaDEPout = 0.750;    % DEP outer station
        % ylocDEP = linspace(etaDEPin,etaDEPout,numDEP/2) .* sum(ac.wing.span);
        % ac.pod.id        = 'pod';   % component name
        % ac.pod.type      = 'POD';   % OpenVSP component type
        % ac.pod.utess     = 5;       % u tessellation
        % ac.pod.wtess     = 5;       % u tessellation
        % for ii = 1:numDEP/2
        %     ac.pod.length(ii)    = 0.452;                    % length
        %     ac.pod.fineratio(ii) = 1.65;                     % fineratio
        %     ac.pod.xloc(ii)      = 8.685;                  % X local coordinate
        %     ac.pod.yloc(ii)      = ylocDEP(ii);            % Y local coordinate
        %     ac.pod.zloc(ii)      = 1.50;                   % Z local coordinate
        %     ac.pod.yrot(ii)      = 0;                      % Y global rotation
        % end
        % % TIP propeller
        % ac.pod.length(ii+1)    = 3.16;                   % length
        % ac.pod.fineratio(ii+1) = 3.46;                    % fineratio
        % ac.pod.xloc(ii+1)      = 8.685;                   % X local coordinate
        % ac.pod.yloc(ii+1)      = sum(ac.wing.span);     % Y local coordinate
        % ac.pod.zloc(ii+1)      = 1.50;                  % Z local coordinate
        % ac.pod.yrot(ii+1)      = 0;                     % Y global rotation

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

        % ac.prop.id       = 'propeller';         % component name
        % ac.prop.type     = 'Disk';              % OpenVSP component type
        % ac.prop.utess    = 15;                  % u tessellation
        % ac.prop.wtess    = 15;                  % w tessellation
        % for ii = 1:numDEP/2
        %     ac.prop.xloc(ii)     = ac.pod.xloc(ii); % X global coordinate
        %     ac.prop.yloc(ii)     = ac.pod.yloc(ii);     % Y global coordinate
        %     ac.prop.zloc(ii)     = ac.pod.zloc(ii); % Z global coordinate
        %     ac.prop.yrot(ii)     = ac.pod.yrot(ii);  % Y global rotation
        %     ac.prop.diameter(ii) = 2.20;             % diameter
        % end
        % % TIP propeller
        % ac.prop.xloc(ii+1)     = ac.pod.xloc(ii+1);   % X global coordinate
        % ac.prop.yloc(ii+1)     = ac.pod.yloc(ii+1);       % Y global coordinate
        % ac.prop.zloc(ii+1)     = ac.pod.zloc(ii+1);   % Z global coordinate
        % ac.prop.yrot(ii+1)     = ac.pod.yrot(ii+1);    % Y global rotation
        % ac.prop.diameter(ii+1) = 3.93;                 % diameter

        %% configuration cases
        ac.config = 'W';                  % configuration to be analyzed
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

        ac.movables = '';
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

        % set.h                       = 10000;          % altitude in meters
        % [set.T, set.a, set.p, set.rho] = ...                % freestream thermodynamic
        %     atmosisa(set.h);                                % characteristics
        % set.mu                      = ...                   % freestream dynamic viscosity
        %     (1.4580e-06*set.T^(1.5))/(set.T+110.4);
        % set.Mach                    = 0.82;                   % freestream Mach number
        % set.Vinf                    = set.Mach*set.a;       % freestream speed
        % % set.rho                     = 0.653;                  % freestream density
        % % set.Vinf                    = 128.6;                 % freestream speed
        % set.Cref                    = ac.wing.mac;            % reference chord
        % set.ReCref                  = ...                   % freestream Re based
        %     (set.rho*set.Vinf*set.Cref)/set.mu;             % on Cref
        % % set.ReCref                  = 1e7;                    % freestream Re based on Cref
        % set.Sref                    = 2*sum(ac.wing.S);       % reference surface
        % set.Bref                    = 2*sum(ac.wing.span);    % reference span
        % set.X_cg                    = ac.wing.macxloc+...     % aircraft center of gravity (PDV 9.64)
        %     0.35*ac.wing.mac;                                 % global coordinates
        % set.Y_cg                    = 0;
        % set.Z_cg                    = 0;
        % 
        % set.AoA                     = -2:2:6;                % attack angle(s)
        % set.Beta                    = 0;                      % sideslip angle(s)
        % 
        % set.ClMax                   = -1;                     % max sectional lift coefficient
        % % -1: no limit
        % set.MaxTurningAngle         = -1;                     % ???
        % 
        % switch set.Beta                                       % imposing symmetry only if
        %     case 0                                            % set.Beta = 0
        %         set.Symmetry        = 'Y';                    %
        %     otherwise                                         % plane of symmetry:
        %         set.Symmetry        = 'NO';                   % NO: no; X: YZ; Y: XZ; Z: XY
        % end
        % 
        % set.FarDist                 = -1;                     % distance where wakes return to freestream:
        % % -1: VSPAERO computes the distance
        % set.NumWakeNodes            = -1;                     % number of wake nodes
        % % -1: free
        % set.WakeIters               = 3;                      % number of iterations
        % 
        % if propON == 1                                        % enable all rotors
        %     set.NumberOfRotors = numel(ac.prop.yloc);         % set geometric and aerodyanmic data arrays
        %     set.xloc = ac.prop.xloc;                          % prop global X location
        %     set.yloc = ac.prop.yloc;                          % prop global Y location
        %     set.zloc = ac.prop.zloc;                          % prop global Z location
        %     set.xdir = cosd(ac.prop.yrot);                    % prop X unit normal
        %     set.ydir = zeros(1,set.NumberOfRotors);           % prop Y unit normal
        %     set.zdir = sind(ac.prop.yrot);                    % prop Z unit normal
        %     set.dia = ac.prop.diameter;                       % prop diameter
        %     set.hub = 0.1.*set.dia;                           % prop hub
        %     
        %     % set RPM to zero to disable that rotor into VSPaero
        %     set.rpm = [0 * ones(1,set.NumberOfRotors-1), ... % DEP props RPM
        %         -990];                                          % TIP props RPM
        %     set.ct = [0.21 * ones(1,set.NumberOfRotors-1), ...   % DEP props CT
        %         0.35];                                           % TIP props CT
        %     set.cp = [0.190 * ones(1,set.NumberOfRotors-1), ...   % DEP props CP
        %         0.768];                                          % TIP props CP
        %     
        % else
        %     set.NumberOfRotors          = 0;                  % prop-off condition
        % end
        % 
        % set.movables            = ac.movables;                % enable analysis with movables
        % set.wingID = ac.wing.id;
        % % set.horID = ac.hor.id;
        % % set.flapRot             = ac.wing.mov.rotation;
        % % set.elevRot             = ac.hor.mov.rotation;
        % 
        % set.VRTXLift                = 'Y';                    % ???
        % % Y: yes; N: no
        % set.LESuction               = 'Y';                    % leading edge suction
        % % Y: yes; N: no
        % set.KT                      = 'Y';                    % 2nd order Karman-Tsien Mach correction
        % % Y: yes; N: no
        % 
        % %% making vspaero based on given setup
        % 
        % disp(' ');
        % disp(' ');
        % str = '### creating vspaero...';
        % disp(' ');
        % disp(str);
        % 
        % MakeVspaero(ac.name,set);
        % 
        % disp(' ');
        % str = '### ... vspaero created';
        % disp(' ');
        % disp(str);
        % 
        % %% running vspaero.exe
        % 
        % disp(' ');
        % disp(' ');
        % str = '### running solver...';
        % disp(' ');
        % disp(str);
        % disp(' ');
        % disp(' ');
        % 
        % num_cpu         = 4;                                     % cpus number
        % vspaero_exe     = ['vspaero -omp ' num2str(num_cpu)...   % command
        %     ' ' ac.name '_DegenGeom'];
        % try
        %     ex = dos(vspaero_exe);                                % executing command
        %     if ex ~= 0
        %         warning('VSPaero failed to execute')
        %     end
        % catch
        %     warning('VSPaero failed to execute')
        % end
        % 
        % %% importing and plotting results and then moving all i/o files of both vspscript.exe
        % % and vspaero.exe in results directory
        % 
        % res = ImportResults(ac);
        % save([ac.name,'_res.mat'],'res')
        % 
        % PlotResults(set,res,ac.name,ac.config);
        % 
         MoveResults(ac.name,res_dir);
        % movefile([ac.name,'_res.mat'],res_dir)
        % 
        % %% Stability analysis about alfa = 0° and beta = 0°
        % if ~exist('stab','dir')
        %     mkdir('stab')
        % end
        % copyfile([res_dir,'\',ac.name '_DegenGeom.csv'],pwd)
        % 
        % commandwindow % focus on the command window, leaving the plots behind
        % 
        % set.AoA = 0;
        % set.Beta = 0;
        % set.Symmetry = 'NO'; % symmetry must be set off on stability calculation
        % 
        % disp(' ');
        % disp(' ');
        % str = '### creating vspaero for stability analysis about alfa=0 and beta=0...';
        % disp(' ');
        % disp(str);
        % 
        % MakeVspaero(ac.name,set);
        % 
        % disp(' ');
        % str = '### ... vspaero created';
        % disp(' ');
        % disp(str);
        % 
        % disp(' ');
        % disp(' ');
        % str = '### running solver for stability analysis...';
        % disp(' ');
        % disp(str);
        % disp(' ');
        % disp(' ');
        % 
        % vspaero_exe     = ['vspaero -omp ' num2str(num_cpu)...   % command
        %     ' -stab ' ac.name '_DegenGeom'];
        % try
        %     ex = dos(vspaero_exe);                                % executing command
        %     if ex ~= 0
        %         warning('VSPaero failed to perform stability calculations')
        %     end
        % catch
        %     warning('VSPaero failed to perform stability calculations')
        % end
        % 
        % movefile([ac.name '_DegenGeom.*'],'stab');
        % movefile('stab',res_dir)

        fprintf('\n\n\n\n');
        
        % Update progress bar
        stepTimeHistory(file_index) = toc;
        avgStepTime = mean(stepTimeHistory,'omitnan');
        timeRemaining = avgStepTime * (steps-file_index);
        estimTime = datestr(seconds(timeRemaining),'HH:MM:SS');
        myString = [num2str(file_index), ' of ', num2str(steps), ...
            ' completed - ', estimTime, ' remaining'];
        waitbar(file_index/steps,progBar,sprintf('%s',myString))
        
                               
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
delete(progBar)

if file_index < steps
   disp('-------------------------------------------------')
   disp(['LOOP EXECUTION ABORTED AT CONFIGURATION no. ',num2str(file_index)])
   disp('-------------------------------------------------') 
end

disp('END')

%% Moving remaining file in results directory
diary off
movefile('logfile.txt',main_dir);
movefile('Legend.xlsx',main_dir);

