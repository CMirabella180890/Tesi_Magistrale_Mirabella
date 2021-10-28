function ShapePanel(fileid,component,i)

% Shapes panels.
%
% fileid   : id file you're writing in
% component: structure containing the particular panel you want to model
% i        : index of the particular panel you're modeling. i=1 is the root
%            panel

fprintf(fileid,['// ' component.id ' panel ' num2str(i) '\r\n\r\n']);

% choosing shape driver parameters, from which the other depends
fprintf(fileid,'%s\r\n\r\n',['SetDriverGroup(' component.id ',' num2str(i) ...
    ',SPAN_WSECT_DRIVER,ROOTC_WSECT_DRIVER,TIPC_WSECT_DRIVER);']);

% modeling the shape of the panel
fprintf(fileid,'%s\r\n',['SetParmVal(' component.id ',"Root_Chord","XSec_'...
    num2str(i) '",' num2str(component.croot(i)) ');']);
fprintf(fileid,'%s\r\n',['SetParmVal(' component.id ',"Tip_Chord","XSec_'...
    num2str(i) '",' num2str(component.ctip(i)) ');']);
fprintf(fileid,'%s\r\n',['SetParmVal(' component.id ',"Span","XSec_'...
    num2str(i) '",' num2str(component.span(i)) ');']);
fprintf(fileid,'%s\r\n',['SetParmVal(' component.id ',"Sweep","XSec_'...
    num2str(i) '",' num2str(component.sweep(i)) ');']);
fprintf(fileid,'%s\r\n',['SetParmVal(' component.id ',"Sweep_Location","XSec_'...
    num2str(i) '",' num2str(component.sweeploc(i)) ');']);
fprintf(fileid,'%s\r\n',['SetParmVal(' component.id ',"Sec_Sweep_Location","XSec_'...
    num2str(i) '",' num2str(component.secsweeploc(i)) ');']);
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' component.id ',"Dihedral","XSec_'...
    num2str(i) '",' num2str(component.dihedral(i)) ');']);
%tessellation spanwise y direction
	fprintf(fileid,'%s\r\n',['SetParmVal(' component.id ',"SectTess_U","XSec_' num2str(i) '",' ...
    num2str(component.utess(i)) ');']);
% fprintf(fileid,'%s\r\n',['SetParmVal(' component.id ',"Twist","XSec_'...
%     num2str(i-1) '",' num2str(component.twist(i)) ');']);
% fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' component.id ',"Twist_Location","XSec_'...
%     num2str(i-1) '",' num2str(component.twistloc(i)) ');']);

fprintf(fileid,'Update(); \r\n\r');       % updating after geometry changes

end