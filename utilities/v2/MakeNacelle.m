function MakeNacelle(fileid,nacelle,i)

% Makes nacelles.
%
% file   : VSPscript file id you're writing in
% nacelle: nacelles data
% i      : index of the particular nacelle you're creating; left wing has
%          odd index; lower the index, closer the nacelle to the fuselage

fprintf(fileid,['// ' nacelle.id ' ' num2str(i) '\r\n\r\n']);    % commenting

% creating the nacelle
fprintf(fileid,'%s\r\n\r\n',['string ' nacelle.id '_' num2str(i)...
    ' = AddGeom("' nacelle.type '");']);

% setting width, height, length and thikness/chord
fprintf(fileid,'%s\r\n',['SetParmVal(' nacelle.id '_' num2str(i)...
    ',"Width","Design",' num2str(nacelle.width) ');']);
fprintf(fileid,'%s\r\n',['SetParmVal(' nacelle.id '_' num2str(i)...
    ',"Height","Design",' num2str(nacelle.height) ');']);
fprintf(fileid,'%s\r\n',['SetParmVal(' nacelle.id '_' num2str(i)...
    ',"Chord","Design",' num2str(nacelle.chord) ');']);
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' nacelle.id '_' num2str(i)...
    ',"ThickChord","Design",' num2str(nacelle.thickchord) ');']);

fprintf(fileid,'Update(); \r\n\r\n');     % updating after geometry changes

% setting location and orientation
fprintf(fileid,'%s\r\n',['SetParmVal(' nacelle.id '_' num2str(i)...
    ',"X_Rel_Location","XForm",' num2str(nacelle.xloc(round(i/2))) ');']);
fprintf(fileid,'%s\r\n',['SetParmVal(' nacelle.id '_' num2str(i)...
    ',"Y_Rel_Location","XForm",' num2str(((-1)^i)*nacelle.yloc(round(i/2))) ');']);
fprintf(fileid,'%s\r\n',['SetParmVal(' nacelle.id '_' num2str(i)...
    ',"Z_Rel_Location","XForm",' num2str(nacelle.zloc(round(i/2))) ');']);
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' nacelle.id '_' num2str(i)...
    ',"Y_Rel_Rotation","XForm",' num2str(nacelle.yrot(round(i/2))) ');']);

fprintf(fileid,'Update(); \r\n\r\n');     % updating after geometry changes

% setting tesselation
fprintf(fileid,'%s\r\n',['SetParmVal(' nacelle.id '_' num2str(i)...
    ',"Tess_U","Shape",' num2str(nacelle.utess) ');']);
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' nacelle.id '_' num2str(i)...
    ',"Tess_W","Shape",' num2str(nacelle.wtess) ');']);

fprintf(fileid,'Update(); \r\n\r\n');      % updating after geometry changes

end