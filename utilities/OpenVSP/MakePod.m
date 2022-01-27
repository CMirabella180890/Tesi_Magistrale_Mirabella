function MakePod(fileid,pod,i)

% Makes pods.
%
% file: id file you're writing in
% pod : pods data
% i   : index of the particular pod you're creating; left wing has odd
%       index; lower the index, closer the pod to the fuselage

fprintf(fileid,['// ' pod.id ' ' num2str(i) '\r\n\r\n']);    % commenting

% creating the pod
fprintf(fileid,'%s\r\n\r\n',['string ' pod.id '_' num2str(i)...
    ' = AddGeom("' pod.type '");']);

%shade
fprintf(fileid,'%s\r\n\r\n',['SetGeomDrawType(' wing.id ', GEOM_DRAW_SHADE );']);
                                                                           % Make component appear as shaded

% setting width, height, length and thikness/chord
fprintf(fileid,'%s\r\n',['SetParmVal(' pod.id '_' num2str(i)...
    ',"Length","Design",' num2str(pod.length(round(i/2))) ');']);
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' pod.id '_' num2str(i)...
    ',"FineRatio","Design",' num2str(pod.fineratio(round(i/2))) ');']);

fprintf(fileid,'Update(); \r\n\r\n');       % updating after geometry changes

% setting location and orientation
fprintf(fileid,'%s\r\n',['SetParmVal(' pod.id '_' num2str(i)...
    ',"X_Rel_Location","XForm",' num2str(pod.xloc(round(i/2))) ');']);
fprintf(fileid,'%s\r\n',['SetParmVal(' pod.id '_' num2str(i)...
    ',"Y_Rel_Location","XForm",' num2str(((-1)^i)*pod.yloc(round(i/2))) ');']);
fprintf(fileid,'%s\r\n',['SetParmVal(' pod.id '_' num2str(i)...
    ',"Z_Rel_Location","XForm",' num2str(pod.zloc(round(i/2))) ');']);
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' pod.id '_' num2str(i)...
    ',"Y_Rel_Rotation","XForm",' num2str(pod.yrot(round(i/2))) ');']);

fprintf(fileid,'Update(); \r\n\r\n');       % updating after geometry changes

% setting tesselation
fprintf(fileid,'%s\r\n',['SetParmVal(' pod.id '_' num2str(i)...
    ',"Tess_U","Shape",' num2str(pod.utess) ');']);
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' pod.id '_' num2str(i)...
    ',"Tess_W","Shape",' num2str(pod.wtess) ');']);

fprintf(fileid,'Update(); \r\n\r\n');       % updating after geometry changes

end
