function MakeFuselage(fileid,fuselage)

% Makes fuselages.
%
% fileid   : id file you're writing in
% fuselage : fuselage data

fprintf(fileid,'// fuselage \r\n\r\n');     % commenting

% creating the fuselage
fprintf(fileid,'%s\r\n\r\n',['string ' fuselage.id ' = AddGeom("' fuselage.type '");']);
%shade
fprintf(fileid,'%s\r\n\r\n',['SetGeomDrawType(' fuselage.id ', GEOM_DRAW_SHADE );']);                                                                       % Make component appear as shaded

% setting length and diameter
fprintf(fileid,'%s\r\n',['SetParmVal(' fuselage.id ',"Length","Design",' ...
    num2str(fuselage.length) ');']);
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' fuselage.id ',"Diameter","Design",' ...
    num2str(fuselage.diameter) ');']);

fprintf(fileid,'Update(); \r\n\r\n');       % updating after geometry changes

%% setting tessellation
fprintf(fileid,'%s\r\n',['SetParmVal(' fuselage.id ',"Tess_U","Shape",' ...
    num2str(fuselage.utess) ');']);
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' fuselage.id ',"Tess_W","Shape",' ...
    num2str(fuselage.wtess) ');']);

fprintf(fileid,'Update(); \r\n\r\n');       % updating after geometry changes

end
