function MakeDegenGeom(fileid,name)

% Integrates the vsprcipt with the lines vspscript.exe needs for
% creating the degenrate geometry.
%
% fileid: id file you're writing in
% name  : aircraft name

% setting name and extension of output file
fprintf(fileid,'%s\r\n',['SetComputationFileName(DEGEN_GEOM_CSV_TYPE,"'...
    name '_DegenGeom.csv");']);  

% always degenerating all created components
fprintf(fileid,'ComputeDegenGeom(SET_ALL,DEGEN_GEOM_CSV_TYPE);\r\n\r\n');

%matlab degen geom
fprintf(fileid,'%s\r\n',['SetComputationFileName(DEGEN_GEOM_M_TYPE,"'...
    name '_DegenGeom.m");']);  

% always degenerating all created components
fprintf(fileid,'ComputeDegenGeom(SET_ALL,DEGEN_GEOM_M_TYPE);\r\n\r\n');

end