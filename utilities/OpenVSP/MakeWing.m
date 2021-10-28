function MakeWing(fileid,wing)

% Makes wings, horizontal tails and vertical tails.
%
% fileid:file id you're writing in
% wing  : wing data

nsec = numel(wing.twist);

fprintf(fileid,['// ' wing.id '\r\n\r\n']);    % commenting

% creating one-panel wing
fprintf(fileid,'%s\r\n\r\n',['string ' wing.id ' = AddGeom("' wing.type '");']);

fprintf(fileid,'%s\r\n\r\n',['SetGeomName(' wing.id ',"' wing.id '");']);


% specifying plane of symmetry, if any
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' wing.id ',"Sym_Planar_Flag","Sym",' ...
    num2str(wing.sym) ');']);

fprintf(fileid,'Update(); \r\n\r\n');               % updating after geometry changes

% setting location and orientation
fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"X_Rel_Location","XForm",' ...
    num2str(wing.xloc) ');']);
fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"Z_Rel_Location","XForm",' ...
    num2str(wing.zloc) ');']);

% assign zero to global rotations if these are undefined
if ~isfield(wing,'xrot')
    wing.xrot = 0;
end
if ~isfield(wing,'yrot')
    wing.yrot = 0;
end
if ~isfield(wing,'zrot')
    wing.zrot = 0;
end

fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"X_Rel_Rotation","XForm",' ...
    num2str(wing.xrot) ');']);
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' wing.id ',"Y_Rel_Rotation","XForm",' ...
    num2str(wing.yrot) ');']);
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' wing.id ',"Z_Rel_Rotation","XForm",' ...
    num2str(wing.zrot) ');']);

fprintf(fileid,'Update(); \r\n\r\n');               % updating after geometry changes

%% panels

% adding panels, if necessary
if numel(wing.span)>1
    for i = 1:numel(wing.span)-1
        fprintf(fileid,'%s\r\n',['InsertXSec(' wing.id ',1,' wing.airfoil ');']);
    end
    fprintf(fileid,'\r\n');
    fprintf(fileid,'Update(); \r\n\r\n');           % updating after geometry changes
end

for i = 1:numel(wing.span)
    ShapePanel(fileid,wing,i);                 % shaping each panel
end

%% airfoils

% changing airfoils type
if strcmp(wing.airfoil,'XS_FOUR_SERIES') == 0
    for i = 1:numel(wing.airfoil)
        fprintf(fileid,'%s\r\n',['ChangeXSecShape(GetXSecSurf(' wing.id...
            ',1),' num2str(i-1) ',' wing.airfoil ');']);
    end
    
    fprintf(fileid,'\r\n');
    fprintf(fileid,'Update(); \r\n\r\n');           % updating after geometry changes
end

% changing airfoil characteristics
for i = 1:nsec
    ShapeAirfoil(fileid,wing,i);               % shaping each panel
end

%% setting tessellation
% chordwise direction x
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' wing.id ',"Tess_W","Shape",' ...
    num2str(wing.wtess) ');']);

% if cap tesselation, l.e. or t.e. clusters are undefined, provide default values
if ~isfield(wing,'captess')
    wing.captess = 9;
end
if ~isfield(wing,'lecluster')
    wing.lecluster = 0.15;
end
if ~isfield(wing,'tecluster')
    wing.tecluster = 0.15;
end

fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' wing.id ',"CapUMinTess","EndCap",' ...
    num2str(wing.captess) ');']);
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' wing.id ',"LECluster","WingGeom",' ...
    num2str(wing.lecluster) ');']);
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' wing.id ',"TECluster","WingGeom",' ...
    num2str(wing.tecluster) ');']);

fprintf(fileid,'Update(); \r\n\r\n');               % updating after geometry changes

%% Blending panels
if isfield(wing,'blending')
    if wing.blending == 1
        
        fprintf(fileid,'%s\r\n\r\n',['// Blending ', wing.id, ' sections']);
        
        for i = 1:nsec
            
            if i > 1 % inboard panels blending
                % leading edges
                fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"InLEMode","' ['XSec_',num2str(i-1)] '",' num2str(wing.blend.InLEMode(i)) ');']);
                if wing.blend.InLEMode(i) ~= 0 % mode is not FREE
                    fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"InLEStrength","' ['XSec_',num2str(i-1)] '",' num2str(wing.blend.InLEStrength(i)) ');']);
                    if wing.blend.InLEMode(i) == 1 % mode is ANGLES
                        fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"InLESweep","' ['XSec_',num2str(i-1)] '",' num2str(wing.blend.InLESweep(i)) ');']);
                        fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"InLEDihedral","' ['XSec_',num2str(i-1)] '",' num2str(wing.blend.InLEDihedral(i)) ');']);
                    end
                end
                % trailing edges
                fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"InTEMode","' ['XSec_',num2str(i-1)] '",' num2str(wing.blend.InTEMode(i)) ');']);
                if wing.blend.InTEMode(i) ~= 0 % mode is not FREE
                    if wing.blend.InTEMode(i) == 1 % mode is ANGLES
                        fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"InTESweep","' ['XSec_',num2str(i-1)] '",' num2str(wing.blend.InTESweep(i)) ');']);
                        fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"InTEDihedral","' ['XSec_',num2str(i-1)] '",' num2str(wing.blend.InTEDihedral(i)) ');']);
                    end
                    fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"InTEStrength","' ['XSec_',num2str(i-1)] '",' num2str(wing.blend.InTEStrength(i)) ');']);
                end
            end
            
            if i < nsec % outboard panels blending
                % leading edges
                fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"OutLEMode","' ['XSec_',num2str(i-1)] '",' num2str(wing.blend.OutLEMode(i)) ');']);
                if wing.blend.OutLEMode(i) ~= 0 % mode is not FREE
                    if wing.blend.OutLEMode(i) == 1 % mode is ANGLES
                        fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"OutLESweep","' ['XSec_',num2str(i-1)] '",' num2str(wing.blend.OutLESweep(i)) ');']);
                        fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"OutLEDihedral","' ['XSec_',num2str(i-1)] '",' num2str(wing.blend.OutLEDihedral(i)) ');']);
                    end
                    fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"OutLEStrength","' ['XSec_',num2str(i-1)] '",' num2str(wing.blend.OutLEStrength(i)) ');']);
                end
                % trailing edges
                fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"OutTEMode","' ['XSec_',num2str(i-1)] '",' num2str(wing.blend.OutTEMode(i)) ');']);
                if wing.blend.OutTEMode(i) ~= 0 % mode is not FREE
                    if wing.blend.OutTEMode(i) == 1 % mode is ANGLES
                        fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"OutTESweep","' ['XSec_',num2str(i-1)] '",' num2str(wing.blend.OutTESweep(i)) ');']);
                        fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"OutTEDihedral","' ['XSec_',num2str(i-1)] '",' num2str(wing.blend.OutTEDihedral(i)) ');']);
                    end
                    fprintf(fileid,'%s\r\n',['SetParmVal(' wing.id ',"OutTEStrength","' ['XSec_',num2str(i-1)] '",' num2str(wing.blend.OutTEStrength(i)) ');']);
                end
            end
            
            fprintf(fileid,'\r\n');
        end
        
        fprintf(fileid,'Update(); \r\n\r\n');  % updating after geometry changes
    end
end

end
