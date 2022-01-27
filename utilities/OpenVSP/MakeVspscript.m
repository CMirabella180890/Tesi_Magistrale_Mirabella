function MakeVspscript(ac)

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


fileid = fopen([ac.name '.vspscript'],'w+');    % opening vspscript file

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

%% MOVABLES
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
%         for i=1:ac.wing.mov.number
%         MakeMovables(fileid,ac.wing)
%         end
        MakeMovables(fileid,ac.wing)
        
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
%% END MOVABLES


%DegenGeom function
MakeDegenGeom(fileid,ac.name);

% Calculate parasite drag NOW PERFORMED WITH ANOTHER SCRIPT FILE
% CalcParasiteDrag(fileid,ac)

% Saving .vsp3 model
MakeVSP3(fileid,ac.name);

% mass properties
% accuracy = 50;
% ComputeMassProps(fileid,ac.name,accuracy);

fprintf(fileid,'}');                                 % ends vspscript file
fclose(fileid);                                      % closes vspscript file

end
