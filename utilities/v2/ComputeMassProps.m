function ComputeMassProps(fileid,name,accuracy)

% Computes mass properties.
%
% fileid  : if file you're wrinting in
% name    : aircraft name
% accuracy: accuracy of the computation

% setting name and extensione of output file
fprintf(fileid,'%s\r\n',['SetComputationFileName(MASS_PROP_TXT_TYPE,"'...
    name '_MassProps.txt");']);  

% always computing all created components mass properties and setting accuracy
fprintf(fileid,'%s\r\n\r\n',['ComputeMassProps(SET_ALL,'...
    num2str(accuracy) ');']);

end

