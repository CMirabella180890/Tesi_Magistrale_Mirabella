function MakeMovables(fileid,wing)
%%%makemovables
% This function adds in fileid (vspscript), a movable on a given wing

% Main Input parameters:
% - npanel = total number of wing panels 
% - eta inner = inner control surface pos y/(b/2)
% - eta outer = outer control surface pos y/(b/2)
% - cf/c_inner = inner chord ratio cf/clocal
% - cf/c_outer = outer chord ratio cf/clocal

% Output parameters
% - writing vspscript

%Convetion is:
% eta_vps = 1/(npanel+2) + [1- 2/(npanel+2)*eta_input]

% Cycle for more than one movables per WING component
for w=1:length((wing.mov.type))
    
    switch wing.mov.type(w) %pdv
        case 'f'
            subID = 'flap';
        case 'a'
            subID = 'aileron';
        case 'e'
            subID = 'elevator';
        case 'r'
            subID = 'rudder';
    end
    
   
    
    %% Input from aircraft structure
    npanel = length(wing.span);
    eta_inner = wing.mov.eta_inner(w);
    eta_outer = wing.mov.eta_outer(w);
    cf_c_inner = wing.mov.cf_c_inner(w);
    cf_c_outer = wing.mov.cf_c_outer(w);
    tessellation = wing.mov.tessellation(w);
    
    %Dimensional wing semi span
    semispan = sum(wing.span);
    y_inner = eta_inner*semispan;
    y_outer = eta_outer*semispan;
    
    %%% to be checked PDV
    xvec = zeros (1,npanel+1);
    yvec = zeros (1,npanel+1);
    
    add = 0;
    for i = 1 : npanel
        add = add + wing.span(i);
        xvec(i+1) = add;
    end
    
    for i=1:npanel+1
        yvec(i) = (i)/(npanel+2);
    end
    
    %interpolate on new vector according OPENVsp
    eta_vsp_inner = interp1(xvec,yvec,y_inner);
    eta_vsp_outer = interp1(xvec,yvec,y_outer);
    
    
    
    % Writing in fileid movable
    %TYPE - %OpenVSP Available Option
    % 1-line
    % 2-rectangle
    % 3-ellipse
    % 4-control surface
    % check id to do!!!
    
    %% linear sub_surface
    % %     Print( string("Add SS_Line\n") );
    % %     string ss_line_id = AddSubSurf( wid, SS_LINE );                 // Add Sub Suface Line
    % %     SetParmVal( wid, "Const_Line_Value", "SubSurface_1", 0.4 );     // Change Loc
    
    %% rectangle sub_surface
    % %     Print( string("Add SS_Rectangle\n") );
    % %     string ss_rec_id = AddSubSurf( wid, SS_RECTANGLE );              // Add Sub Surface Rectangle
    % %     SetParmVal( wid, "Center_W", "SS_Rectangle_1", 0.0 );           // Change Loc/Size
    % %     SetParmVal( wid, "U_Length", "SS_Rectangle_1", 0.1 );
    % %     SetParmVal( wid, "W_Length", "SS_Rectangle_1", 0.1 );
    % %
    % %     Update();
    
    
    % Adding a movable "Control Surface" Type
    fprintf(fileid,'%s\r\n',['string ' subID ' = AddSubSurf( ' wing.id ' , SS_CONTROL );']);
    
    fprintf(fileid,'%s\r\n\r\n',['SetSubSurfName(' wing.id ',' subID ', "' subID '" );']);
    
    % eta inner and eta outer
    fprintf(fileid,'%s\r\n',['SetParmVal( ' wing.id ' , "UStart", "SS_Control_' num2str(w) '", ' num2str(eta_vsp_inner) ');']);
    fprintf(fileid,'%s\r\n',['SetParmVal( ' wing.id ' , "UStart", "SS_Control_' num2str(w) '", ' num2str(eta_vsp_inner) ');']);
    
    fprintf(fileid,'%s\r\n',['SetParmVal( ' wing.id ' , "UEnd", "SS_Control_' num2str(w) '", ' num2str(eta_vsp_outer) ');']);
    %SS_Control_1 "_1" must be indexed based on number of control surfaces before and after
    
    % cf/c: at the moment is constant for overall control surface spanwise
    fprintf(fileid,'%s\r\n',['SetParmVal( ' wing.id ' , "Length_C_Start", "SS_Control_' num2str(w) '", ' num2str(cf_c_inner) ');']);
    
    % tasselation on edges
    fprintf(fileid,'%s\r\n',['SetParmVal( ' wing.id ' , "Tess_Num", "SS_Control_' num2str(w) '", ' num2str(tessellation) ');']);
    
    %Update control surface
    fprintf(fileid,'Update(); \r\n\r\n');               % updating after geometry changes
    
end

end