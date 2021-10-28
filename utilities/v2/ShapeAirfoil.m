function ShapeAirfoil(fileid,component,i)

% Shapes airfoils. It allows you to model only NACA 4 and 5 digit ones.
%
% file     : id file you're writing in
% component: structure containing the particular airfoil you want to model
% i        : index of the particular airfoil you're modeling. i=1 is the root
%            airfoil

fprintf(fileid,['// ' component.id ' airfoil ' num2str(i) '\r\n\r\n']);

% setting camber, camber location and thickness
switch component.airfoil
    case 'XS_FOUR_SERIES'
        fprintf(fileid,'%s\r\n',['SetParmVal(' component.id ',"Camber","XSecCurve_'...
            num2str(i-1) '",' num2str(component.camber(i)) ');']);
        fprintf(fileid,'%s\r\n',['SetParmVal(' component.id ',"CamberLoc","XSecCurve_'...
            num2str(i-1) '",' num2str(component.camberloc(i)) ');']);
        fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' component.id ',"ThickChord","XSecCurve_'...
            num2str(i-1) '",' num2str(component.thickchord(i)) ');']);
        
    case 'XS_FIVE_DIGIT'
        fprintf(fileid,'%s\r\n',['SetParmVal(' component.id ',"IdealCl","XSecCurve_'...
            num2str(i-1) '",' num2str(component.idealcl(i)) ');']);
        fprintf(fileid,'%s\r\n',['SetParmVal(' component.id ',"CamberLoc","XSecCurve_'...
            num2str(i-1) '",' num2str(component.camberloc(i)) ');']);
        fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' component.id ',"ThickChord","XSecCurve_'...
            num2str(i-1) '",' num2str(component.thickchord(i)) ');']);
        
    case 'XS_FILE_AIRFOIL'
        fprintf(fileid,'%s\r\n',['ReadFileAirfoil(GetXSec(GetXSecSurf(' component.id...
            ',1),' num2str(i-1) '),"' component.affile{i} '");']);
        if isfield(component,'thickchord') % eventually scale thickness
            fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' component.id ',"ThickChord","XSecCurve_'...
                num2str(i-1) '",' num2str(component.thickchord(i)) ');']);
        end
end

fprintf(fileid,'%s\r\n',['SetParmVal(' component.id ',"Twist","XSec_'...
    num2str(i-1) '",' num2str(component.twist(i)) ');']);
fprintf(fileid,'%s\r\n\r\n',['SetParmVal(' component.id ',"Twist_Location","XSec_'...
    num2str(i-1) '",' num2str(component.twistloc(i)) ');']);

fprintf(fileid,'Update(); \r\n\r');       % updating after geometry changes

end