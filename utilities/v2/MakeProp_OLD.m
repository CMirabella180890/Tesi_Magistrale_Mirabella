function MakeProp(fileid,prop,i)

% Makes propellers.
%
% fileid: id file you're writing in
% prop  : pods data
% i     : index of the particular propeller you're creating; left wing has
%         odd index; lower the index, closer the propeller to the fuselage


fprintf(fileid,['// ' prop.id ' ' num2str(i) '\r\n\r\n']);    % commenting

% creating the propellers
fprintf(fileid,'%s\r\n\r\n',['string ' prop.id '_' num2str(i)...
    ' = AddGeom("' prop.type '");']);

% setting location and orientation
fprintf(fileid,'%s\r\n',['SetParmVal(' prop.id '_' num2str(i)...
    ',"X_Rel_Location","XForm",' num2str(prop.xloc(round(i/2))) ');']);
fprintf(fileid,'%s\r\n',['SetParmVal(' prop.id '_' num2str(i)...
    ',"Y_Rel_Location","XForm",' num2str(((-1)^i)*prop.yloc(round(i/2))) ');']);
fprintf(fileid,'%s\r\n',['SetParmVal(' prop.id '_' num2str(i)...
    ',"Z_Rel_Location","XForm",' num2str(prop.zloc(round(i/2))) ');']);
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' prop.id '_' num2str(i)...
    ',"Y_Rel_Rotation","XForm",' num2str(prop.yrot(round(i/2))) ');']);

fprintf(fileid,'Update();\r\n\r\n');       % updating after geometry changes

% setting length
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' prop.id '_' num2str(i)...
    ',"Diameter","Design",' num2str(prop.diameter) ');']);

fprintf(fileid,'Update();\r\n\r\n');       % updating after geometry changes

%% setting tesselation
fprintf(fileid,'%s\r\n',['SetParmVal(' prop.id '_' num2str(i)...
    ',"Tess_U","Shape",' num2str(prop.utess) ');']);
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' prop.id '_' num2str(i)...
    ',"Tess_W","Shape",' num2str(prop.wtess) ');']);

fprintf(fileid,'Update(); \r\n\r\n');      % updating after geometry changes

end
