function MakeVSP3(fileid,name)

% Creates the .vsp3 fiel
%
% fileid: id file you're writing in
% name  : aircraft name
% Example:
%    //==== Save Vehicle to File ====//
%    Print( string( "--->Saving VSP model\n" ) );
%    string fname = "TestSaveLoad.vsp3";
%    WriteVSPFile( fname, SET_ALL );
%

% adding a comment
fprintf(fileid,'//==== Save Vehicle to File ====//\r\n');

% setting filename and extension    
fprintf(fileid,'%s\r\n',['string fname = "'...
    name '.vsp3";']);    

% writing all on file .vsp3
fprintf(fileid,'WriteVSPFile( fname, SET_ALL );\r\n\r\n');

end