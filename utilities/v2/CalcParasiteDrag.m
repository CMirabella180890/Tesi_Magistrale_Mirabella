function CalcParasiteDrag(ac,airspeed,altitude)
% Write VSPscript file to calculate parasite drag
% Call this function in MakeVspscript function after making the geometric parts
% The parasite drag calculation is independent from VSPAERO, hence the
% reference quantities are not taken from the SET structure defined from
% the main script.
% Please, input everything in SI units.
%
% CalcParasiteDrag(ac) calculates the parasite drag at an airspeed
% of 100 m/s and sea level
%
% CalcParasiteDrag(ac,airspeed) calculates the parasite drag at an
% airspeed defined by the user and sea level
%
% CalcParasiteDrag(ac,airspeed,altitude) calculates the parasite
% drag with both airspeed and altitude defined by the user
%
% WARNING: the initial geometry part is the same of the MAKEVSPSCRIPT
% function. Therefore, any edit in that function should be also made here

if nargin == 1
   airspeed = 100;
   altitude = 0;
elseif nargin == 2
    altitude = 0;
end

disp('### calculating parasite drag...')

% Makes .vspscript file.
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
%
% ac.movables: movables configuration
%           none       = none
%           f       = flap
%           a       = aileron
%           fa      = flap aileron
%           fe      = flap elevator
%           fr      = flap rudder
%           fae     = flap aileron elevator
%           faer    = flap aileron elevator rudder
%           e       = elevator
%           r       = rudder


fileid = fopen([ac.name '_ParasiteDrag.vspscript'],'w+');    % opening vspscript file

fprintf(fileid,'void main()\r\n');              % starting vspscript file

fprintf(fileid,'{\r\n\r\n');

switch ac.config                                   % writing vspecript file
    case 'W'                                    % depending on configuration
        MakeWing(fileid,ac.wing);
    case 'WP'                                    % depending on configuration
        MakeWing(fileid,ac.wing);
        for i = 1:2*length(ac.prop.yloc)
            MakeProp(fileid,ac.prop,i);
        end               
    case 'WB'
        MakeWing(fileid,ac.wing);
        MakeFuselage(fileid,ac.fuse);
    case 'WH'
        MakeWing(fileid,ac.wing);
        MakeWing(fileid,ac.hor);
    case 'WV'
        MakeWing(fileid,ac.wing);
        MakeWing(fileid,ac.ver);
    case 'WHP'
        MakeWing(fileid,ac.wing);
        MakeWing(fileid,ac.hor);
        for i = 1:2*numel(ac.pod.yloc)
            MakePod(fileid,ac.pod,i);
            MakeProp(fileid,ac.prop,i);
        end
    case 'WHF'
        MakeWing(fileid,ac.wing);
        MakeWing(fileid,ac.hor);
        for i = 1:2*numel(ac.nacelle.yloc)
            MakeNacelle(fileid,ac.nacelle,i);
        end
        
    case 'WVP'
        MakeWing(fileid,ac.wing);
        MakeWing(fileid,ac.ver);
        for i = 1:2*numel(ac.pod.yloc)
            MakePod(fileid,ac.pod,i);
            MakeProp(fileid,ac.prop,i);
        end
    case 'WVF'
        MakeWing(fileid,ac.wing);
        MakeWing(fileid,ac.ver);
        for i = 1:2*numel(ac.nacelle.yloc)
            MakeNacelle(fileid,ac.nacelle,i);
        end            
    case 'WBH'
        MakeWing(fileid,ac.wing);
        MakeWing(fileid,ac.hor);
        MakeFuselage(fileid,ac.fuse);
    case 'WBV'
        MakeWing(fileid,ac.wing);
        MakeWing(fileid,ac.ver);
        MakeFuselage(fileid,ac.fuse);
    case 'WBHP'
        MakeWing(fileid,ac.wing);
        MakeWing(fileid,ac.hor);
        MakeFuselage(fileid,ac.fuse);
        for i = 1:2*numel(ac.pod.yloc)
            MakePod(fileid,ac.pod,i);
            MakeProp(fileid,ac.prop,i);
        end
    case 'WBHF'
        MakeWing(fileid,ac.wing);
        MakeWing(fileid,ac.hor);
        MakeFuselage(fileid,ac.fuse);
        for i = 1:2*numel(ac.nacelle.yloc)
            MakeNacelle(fileid,ac.nacelle,i);
        end
        
    case 'WBVP'
        MakeWing(fileid,ac.wing);
        MakeWing(fileid,ac.ver);
        MakeFuselage(fileid,ac.fuse);
        for i = 1:2*numel(ac.pod.yloc)
            MakePod(fileid,ac.pod,i);
            MakeProp(fileid,ac.prop,i);
        end
    case 'WBVF'
        MakeWing(fileid,ac.wing);
        MakeWing(fileid,ac.ver);
        MakeFuselage(fileid,ac.fuse);
        for i = 1:2*numel(ac.nacelle.yloc)
            MakeNacelle(fileid,ac.nacelle,i);
        end
    case 'WBHV'
        MakeWing(fileid,ac.wing);
        MakeWing(fileid,ac.hor);
        MakeWing(fileid,ac.ver);
        MakeFuselage(fileid,ac.fuse);
    case 'WBHVP'
        MakeWing(fileid,ac.wing);
        MakeWing(fileid,ac.hor);
        MakeWing(fileid,ac.ver);
        MakeFuselage(fileid,ac.fuse);
        for i = 1:2*numel(ac.pod.yloc)
            MakePod(fileid,ac.pod,i);
            MakeProp(fileid,ac.prop,i);
        end
    case 'WBHVF'
        MakeWing(fileid,ac.wing);
        MakeWing(fileid,ac.hor);
        MakeWing(fileid,ac.ver);
        MakeFuselage(fileid,ac.fuse);
        for i = 1:2*numel(ac.nacelle.yloc)
            MakeNacelle(fileid,ac.nacelle,i);
        end
end

% MOVABLES
% movables: movables configuration
%           none       = none
%           f       = flap
%           a       = aileron
%           fa      = flap aileron
%           fe      = flap elevator
%           fr      = flap rudder
%           fae     = flap aileron elevator
%           faer    = flap aileron elevator rudder
%           e       = elevator
%           r       = rudder

switch ac.movables
    
    case 'f'                            %flap on wing
        MakeMovables(fileid,ac.wing)
        
    case 'fa'                           %flap and aileron on wing
        MakeMovables(fileid,ac.wing)
        %MakeMovables(fileid,ac.wing)
        
    case 'fe'                          %flap on wing, elevator
        MakeMovables(fileid,ac.wing)
        MakeMovables(fileid,ac.hor)
        
    case 'fr'                          %flap on wing, rudder
        MakeMovables(fileid,ac.wing.f)
        MakeMovables(fileid,ac.ver.r)
        
    case 'fae'                          %flap and aileron on wing, elevator
        MakeMovables(fileid,ac.wing)
        MakeMovables(fileid,ac.hor)
        
    case 'far'                          %flap on wing, elevator
        MakeMovables(fileid,ac.wing)
        MakeMovables(fileid,ac.ver)
        
    case 'fer'                          %flap on wing, elevator
        MakeMovables(fileid,ac.wing)
        MakeMovables(fileid,ac.hor)
        MakeMovables(fileid,ac.ver)
        
    case 'faer'                       %flap and aileron on wing, elevator, rudder
        MakeMovables(fileid,ac.wing)
        %MakeMovables(fileid,ac.wing)
        MakeMovables(fileid,ac.hor)
        MakeMovables(fileid,ac.ver)
        
    case 'e'                          %flap and aileron on wing, elevator
        MakeMovables(fileid,ac.hor)
        
    case 'r'                          %flap and aileron on wing, elevator
        MakeMovables(fileid,ac.ver)
end
% END MOVABLES

fprintf(fileid,'// Parasite Drag Settings \r\n');    % comment
fprintf(fileid,'string PDcon = FindContainer( "ParasiteDragSettings", 0 ); \r\n');
fprintf(fileid,'string groupname = "ParasiteDrag"; \r\n\r\n');

fprintf(fileid,'array<string> fileNameInput = GetStringAnalysisInput( "ParasiteDrag", "FileName" ); \r\n');
fprintf(fileid,['fileNameInput[0] = "', ac.name, '_ParasiteDrag"; \r\n']);
fprintf(fileid,'SetStringAnalysisInput( "ParasiteDrag", "FileName", fileNameInput ); \r\n\r\n');

fprintf(fileid,'array<int> vinfUnitInput = GetIntAnalysisInput( "ParasiteDrag", "VelocityUnit" ); \r\n');
fprintf(fileid,'vinfUnitInput[0] = V_UNIT_M_S; \r\n');
fprintf(fileid,'SetIntAnalysisInput( "ParasiteDrag", "VelocityUnit", vinfUnitInput ); \r\n\r\n');

fprintf(fileid,'array<double> vinfFCinput = GetDoubleAnalysisInput( "ParasiteDrag", "Vinf" ); \r\n');
fprintf(fileid,['vinfFCinput[0] = ', num2str(airspeed),'; \r\n']);
fprintf(fileid,'SetDoubleAnalysisInput( "ParasiteDrag", "Vinf", vinfFCinput ); \r\n\r\n');

fprintf(fileid,'array<int> altUnitInput = GetIntAnalysisInput( "ParasiteDrag", "LengthUnit" ); \r\n');
fprintf(fileid,'altUnitInput[0] = LEN_M; \r\n');
fprintf(fileid,'SetIntAnalysisInput( "ParasiteDrag", "LengthUnit", altUnitInput ); \r\n\r\n');

% I cannot assign SI units to altitude, therefore I convert the SI values in ft
fprintf(fileid,'array<double> altInput = GetDoubleAnalysisInput( "ParasiteDrag", "Altitude" ); \r\n');
fprintf(fileid,['altInput[0] = ', num2str(altitude/.3048),'; \r\n']);
fprintf(fileid,'SetDoubleAnalysisInput( "ParasiteDrag", "Altitude", altInput ); \r\n\r\n');

fprintf(fileid,'array<double> srefInput = GetDoubleAnalysisInput( "ParasiteDrag", "Sref" ); \r\n');
fprintf(fileid,['srefInput[0] = ', num2str(2*sum(ac.wing.S)),'; \r\n']); % 2 half-wing planforms
fprintf(fileid,'SetDoubleAnalysisInput( "ParasiteDrag", "Sref", srefInput ); \r\n\r\n');

fprintf(fileid,'array<int> turbcfEqnChoiceInput = GetIntAnalysisInput( "ParasiteDrag", "TurbCfEqnChoice" ); \r\n');
fprintf(fileid,'turbcfEqnChoiceInput[0] = CF_TURB_IMPLICIT_KARMAN_SCHOENHERR; \r\n');
fprintf(fileid,'SetIntAnalysisInput( "ParasiteDrag", "TurbCfEqnChoice", turbcfEqnChoiceInput ); \r\n\r\n');

fprintf(fileid,'string FFEqnChoice; \r\n');
if contains(ac.config,'W')
    fprintf(fileid,['SetParmVal( ', ac.wing.id, ',  "FFWingEqnType", "ParasiteDragProps", FF_W_DATCOM ); \r\n']);
end
if contains(ac.config,'B')
    fprintf(fileid,['SetParmVal( ', ac.fuse.id, ',  "FFBodyEqnType", "ParasiteDragProps", FF_B_SCHEMENSKY_BODY ); \r\n']);
end
if contains(ac.config,'H')
    fprintf(fileid,['SetParmVal( ', ac.hor.id, ',  "FFWingEqnType", "ParasiteDragProps", FF_W_DATCOM ); \r\n']);
end
if contains(ac.config,'V')
    fprintf(fileid,['SetParmVal( ', ac.ver.id, ',  "FFWingEqnType", "ParasiteDragProps", FF_W_JENKINSON_TAIL ); \r\n']);
end
if contains(ac.config,'P')
    fprintf(fileid,['SetParmVal( ', ac.pod.id, ',  "FFBodyEqnType", "ParasiteDragProps", FF_B_SCHEMENSKY_NACELLE ); \r\n']);
end
fprintf(fileid,'\r\n');

% TODO: include user-defined components form factors
% TODO: include user-defined excrescences

fprintf(fileid,'string ridpd = ExecAnalysis( "ParasiteDrag" ); \r\n');
fprintf(fileid,'array<double>@ dat = GetDoubleResults( ridpd, "Total_CD_Total", 0 ); \r\n');
%  fprintf(fileid,'Print("CD0: " + dat[0]);');
fprintf(fileid,'Print(dat[0]);'); % print only a scalar (total CD)
fprintf(fileid,'\r\n\r\n');


fprintf(fileid,'}');                                 % ends vspscript file
fclose(fileid);                                      % closes vspscript file

end