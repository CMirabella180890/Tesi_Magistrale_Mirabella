function WingOut(Aircraft)
%% Wing Geometry Output function
% This function creates figures of 3-View of the fuselage and 3D Solid PLOT
% of the Aircraft based on here imported airfoils
%
% 1) Auxiliary variables are firstly defined
% 2) 2D plot -> 3-Views of aircraft drafting
% 3) 3D plot -> 3D PLOT of aircraft structure with imported airfoils

comp = 'Wing';

if isfield(Aircraft.Geometry, "Wing") == 1
    disp ([comp, 'exists'])
    disp ('-----------------------------------------------------------------')
    disp ('-----------------------------------------------------------------')
    disp ([comp, 'Geometry utilities end run'])
    
    %% Auxiliary declaration
    %fuselage
    df = Aircraft.Geometry.Fuselage.diameter.value;         % fuselage diameter (m)
    %     dbase = 0.1*df;                                         %fuselage base diameter
    lf = Aircraft.Geometry.Fuselage.length.value;           % fuselage lenght (m)
    %     Aircraft.Geometry.Fuselage.kcockpit.value = 1.0;
    %     Aircraft.Geometry.Fuselage.ktail.value = 1.5;
    %wing
    zpos = Aircraft.Geometry.Wing.zle.value;                                    % wing root chord zeta position (unit)
    xle = Aircraft.Geometry.Wing.xle.value/lf;                                  % wing tip leading edge %of fuselage lenght
    %sweep_le = 0; %TO BE EDITED!!!
    %sweep_le
    sweep_le = Aircraft.Geometry.Wing.sweep_second.value;
    xtip_le = xle + ...
        (Aircraft.Geometry.Wing.b.value/2*...
        tan(sweep_le/57.3))/lf;                       % wing tip leading edge %of fuselage lenght
    Aircraft.Geometry.Wing.ypos.value = 0.0;                                    % %of wing semispan
    %kink 1
    ckink1 = Aircraft.Geometry.Wing.chord_kink_one.value;                                   % kink 1 chord in meters
    ykink1 = Aircraft.Geometry.Wing.panel_span1.value * Aircraft.Geometry.Wing.b.value/2;   % kink 1 wing semispan in meters
    %kink 2
    ckink2 = Aircraft.Geometry.Wing.chord_kink_two.value;                                                                               % kink 2 chord in meters
    ykink2 = (Aircraft.Geometry.Wing.panel_span1.value + Aircraft.Geometry.Wing.panel_span2.value) * Aircraft.Geometry.Wing.b.value/2;  % kink 2 wing semispan in meters
    % Importing airfoils sections for wing, horizontal, vertical...TO BE
    % MODIFIED FROM INPUT data structure!!!
    root_coord_w = importdata('_Airfoil\E216_104.txt'); % Root Airfoil coordiante read from file (x/c & z/c) IRON_Root.txt
    tip_coord_w = importdata('_Airfoil\E216_104.txt'); % Root Airfoil coordiante read from file (x/c & z/c) IRON_Tip.txt
    
    %% TOP-VIEW
    figure('Name',[comp '-Top-View'],'NumberTitle','off');
    hold on
    plot([0 0],xle*lf + ...
        [0, Aircraft.Geometry.Wing.croot.value],'r',"LineWidth",1.5)
    plot([Aircraft.Geometry.Wing.b.value/2, Aircraft.Geometry.Wing.b.value/2],...
        xtip_le*lf +...
        [0, Aircraft.Geometry.Wing.ctip.value],'r',"LineWidth",1.5)
    plot([-Aircraft.Geometry.Wing.b.value/2, -Aircraft.Geometry.Wing.b.value/2],...
        xtip_le*lf + ...
        [0, Aircraft.Geometry.Wing.ctip.value],'r',"LineWidth",1.5)
    plot([0, Aircraft.Geometry.Wing.b.value/2],...
        [xle*lf, xtip_le*lf],...
        'r',"LineWidth",1.5)
    plot([0 -Aircraft.Geometry.Wing.b.value/2],...
        [xle*lf, xtip_le*lf],...
        'r',"LineWidth",1.5)
    plot([0, Aircraft.Geometry.Wing.b.value/2],...
        [xle*lf+Aircraft.Geometry.Wing.croot.value, xtip_le*lf+Aircraft.Geometry.Wing.ctip.value],...
        'r',"LineWidth",1.5)
    plot([0, -Aircraft.Geometry.Wing.b.value/2],...
        [xle*lf+Aircraft.Geometry.Wing.croot.value, xtip_le*lf+Aircraft.Geometry.Wing.ctip.value],...
        'r',"LineWidth",1.5)
    %kink1
    plot([ykink1 ykink1],xle*lf + ykink1*tan(sweep_le/57.3) + ...
        [0, ckink1],'r',"LineWidth",1.5)
    plot([-ykink1 -ykink1],xle*lf + ykink1*tan(sweep_le/57.3) + ...
        [0, ckink1],'r',"LineWidth",1.5)
    %kink2
    plot([ykink2 ykink2],xle*lf + ykink2*tan(sweep_le/57.3) + ...
        [0, ckink2],'r',"LineWidth",1.5)
    plot([-ykink2 -ykink2],xle*lf + ykink2*tan(sweep_le/57.3) + ...
        [0, ckink1],'r',"LineWidth",1.5)
    
    %AILERON %TO BE MODIFIED for kinks WING!!!!
    if isfield(Aircraft.Geometry, 'Aileron') == 1
        %linear chord law
        c0 = Aircraft.Geometry.Wing.croot.value;          % root chord (unit)
        c1 = Aircraft.Geometry.Wing.ctip.value;           % tip chord (unit)
        span0 = Aircraft.Geometry.Wing.yle.value;
        span1 = Aircraft.Geometry.Wing.b.value/2;
        [c_0,c_y] = chord_linear(c0 ,c1, span0, span1);
        
        y_inner = Aircraft.Geometry.Aileron.eta_inner.value * span1;
        y_outer = Aircraft.Geometry.Aileron.eta_outer.value * span1;
        
        c_inner = c_0 - c_y * y_inner;                      %htail chord at inner station (unit)
        c_outer = c_0 - c_y * y_outer;                      %htail chord at inner station (unit)
        
        cf_inner = c_inner * Aircraft.Geometry.Aileron.ca_c_root.value;   %elevetor chord at inner station (unit)
        cf_outer = c_outer* Aircraft.Geometry.Aileron.ca_c_tip.value;    %elevetor chord at inner station (unit)
        
        xle_inner = Aircraft.Geometry.Wing.xle.value +  y_inner * tan(sweep_le/57.3) +...
            c_inner - cf_inner;
        
        xle_outer = Aircraft.Geometry.Wing.xle.value +  y_outer * tan(sweep_le/57.3) +...
            c_outer - cf_outer;
        
        plot([y_inner y_inner],...
            [xle_inner xle_inner+cf_inner],'r',"LineWidth",1.5)
        plot([y_outer y_outer],...
            [xle_outer xle_outer+cf_outer],'r',"LineWidth",1.5)
        plot([y_inner y_outer],...
            [xle_inner xle_outer],'r',"LineWidth",1.5)
        fill ([y_inner y_outer y_outer y_inner],[xle_inner xle_outer xle_outer+cf_outer xle_inner+cf_inner],'r')
        %
        plot([-y_inner -y_inner],...
            [xle_inner xle_inner+cf_inner],'r',"LineWidth",1.5)
        plot([-y_outer -y_outer],...
            [xle_outer xle_outer+cf_outer],'r',"LineWidth",1.5)
        plot([-y_inner -y_outer],...
            [xle_inner xle_outer],'r',"LineWidth",1.5)
        %
        fill ([-y_inner -y_outer -y_outer -y_inner],[xle_inner xle_outer xle_outer+cf_outer xle_inner+cf_inner],'r')
    end
    
    %FLAPS %TO BE MODIFIED for kinks WING!!!!
    if isfield(Aircraft.Geometry, 'Flaps') == 1
        %linear chord law
        c0 = Aircraft.Geometry.Wing.croot.value;          % root chord (unit)
        c1 = Aircraft.Geometry.Wing.ctip.value;           % tip chord (unit)
        span0 = Aircraft.Geometry.Wing.yle.value;
        span1 = Aircraft.Geometry.Wing.b.value/2;
        [c_0,c_y] = chord_linear(c0 ,c1, span0, span1);
        
        y_inner = Aircraft.Geometry.Flaps.eta_inner.value * span1;
        y_outer = Aircraft.Geometry.Flaps.eta_outer.value * span1;
        
        c_inner = c_0 - c_y * y_inner;                      %htail chord at inner station (unit)
        c_outer = c_0 - c_y * y_outer;                      %htail chord at inner station (unit)
        
        cf_inner = c_inner * Aircraft.Geometry.Flaps.cf_c_root.value;   %elevetor chord at inner station (unit)
        cf_outer = c_outer* Aircraft.Geometry.Flaps.cf_c_tip.value;    %elevetor chord at inner station (unit)
        
        xle_inner = Aircraft.Geometry.Wing.xle.value +  y_inner * tan(sweep_le/57.3) +...
            c_inner - cf_inner;
        
        xle_outer = Aircraft.Geometry.Wing.xle.value +  y_outer * tan(sweep_le/57.3) +...
            c_outer - cf_outer;
        
        plot([y_inner y_inner],...
            [xle_inner xle_inner+cf_inner],'r',"LineWidth",1.5)
        plot([y_outer y_outer],...
            [xle_outer xle_outer+cf_outer],'r',"LineWidth",1.5)
        plot([y_inner y_outer],...
            [xle_inner xle_outer],'r',"LineWidth",1.5)
        fill ([y_inner y_outer y_outer y_inner],[xle_inner xle_outer xle_outer+cf_outer xle_inner+cf_inner],'r')
        %
        plot([-y_inner -y_inner],...
            [xle_inner xle_inner+cf_inner],'r',"LineWidth",1.5)
        plot([-y_outer -y_outer],...
            [xle_outer xle_outer+cf_outer],'r',"LineWidth",1.5)
        plot([-y_inner -y_outer],...
            [xle_inner xle_outer],'r',"LineWidth",1.5)
        %
        fill ([-y_inner -y_outer -y_outer -y_inner],[xle_inner xle_outer xle_outer+cf_outer xle_inner+cf_inner],'r')
    end
    hold off
    
    
    grid on
    title([comp 'Top-View'])
    xlabel({'y (m)'})
    ylabel({'x (m)'})
    axis equal
    
    saveas(gcf, [comp '-Top-View.fig'])
    saveas(gcf, [comp '-Top-View.png'])
    
    
    %% SIDE-VIEW
    figure('Name',[comp '-Side-View'],'NumberTitle','off');
    
    plot([xle*lf, xle*lf+Aircraft.Geometry.Wing.croot.value]...
        ,[zpos, zpos],'r',"LineWidth",1.5)                              %root chord (unit)
    hold on

    %Aifoils
    root_coord = root_coord_w;
    tip_coord = tip_coord_w;
    
    %root airfoil
    %x
    root_coord_3D(:,1) = root_coord(:,1).*Aircraft.Geometry.Wing.croot.value +...
        xle.*lf;
    %y
    root_coord_3D(:,2) = root_coord(:,1).*Aircraft.Geometry.Wing.ypos.value*Aircraft.Geometry.Wing.b.value;
    %z
    root_coord_3D(:,3) =  root_coord(:,2).*Aircraft.Geometry.Wing.croot.value +...
        Aircraft.Geometry.Wing.zle.value;
        
    %tip airfoil
    %x
    tip_coord_3D(:,1) = tip_coord(:,1).*Aircraft.Geometry.Wing.ctip.value +...
        xtip_le.*lf;
    %y
    tip_coord_3D(:,2) = Aircraft.Geometry.Wing.b.value/2;
    %z
    tip_coord_3D(:,3) =  tip_coord(:,2).*Aircraft.Geometry.Wing.ctip.value +...
        Aircraft.Geometry.Wing.zle.value+...
        Aircraft.Geometry.Wing.b.value/2*tan(Aircraft.Geometry.Wing.dihedral.value/57.3);
    
     plot(root_coord_3D(:,1),root_coord_3D(:,3),'r',"LineWidth",2)
     plot(tip_coord_3D(:,1),tip_coord_3D(:,3),'r',"LineWidth",2)
    
    hold off
    
    grid on
    title([comp 'Side-View'])
    xlabel({'y (m)'})
    ylabel({'x (m)'})
    axis equal
    
    saveas(gcf, [comp '-Side-View.fig'])
    saveas(gcf, [comp '-Side-View.png'])
    
    
    
    
    %% FRONT-VIEW
    
    figure('Name',[comp '-Front-View'],'NumberTitle','off');
    hold on
%wing
plot([0, Aircraft.Geometry.Wing.b.value/2], [zpos, ...
    zpos+...
    Aircraft.Geometry.Wing.b.value/2*tan(Aircraft.Geometry.Wing.dihedral.value/57.3)],'r',"LineWidth",1.5);
plot([0, -Aircraft.Geometry.Wing.b.value/2], [zpos, ...
    zpos+...
    Aircraft.Geometry.Wing.b.value/2*tan(Aircraft.Geometry.Wing.dihedral.value/57.3)],'r',"LineWidth",1.5);
    
    
    grid on
    title([comp 'Front-View'])
    xlabel({'y (m)'})
    ylabel({'x (m)'})
    axis equal
    
    saveas(gcf, [comp '-Front-View.fig'])
    saveas(gcf, [comp '-Front-View.png'])
    
    
    %% 3D-VIEW
    
    figure('Name',[comp '3D'],'NumberTitle','off');
    
    XXw = [root_coord_3D(:,1), tip_coord_3D(:,1)]';
    YYw = [root_coord_3D(:,2), tip_coord_3D(:,2)]';
    ZZw = [root_coord_3D(:,3), tip_coord_3D(:,3)]';
    
    surf(XXw,YYw,ZZw,'FaceColor','red','EdgeColor','none')
    hold on
    surf(XXw,-YYw,ZZw,'FaceColor','red','EdgeColor','none')
    daspect([1.0 , 1.0, 1.0])
    camlight left; lighting phong

    grid on
    title([comp '3D'])
    xlabel({'y (m)'})
    ylabel({'x (m)'})
    axis equal
       
    saveas(gcf, [comp '3D.fig'])
    saveas(gcf, [comp '3D.png'])      
    
    disp ('-----------------------------------------------------------------')
    disp ('-----------------------------------------------------------------')
    disp ([comp, 'Geometry utilities end run'])

        %% moving results
    movefile([comp '-Top-View.fig'], Aircraft.res_dir);
    movefile([comp '-Top-View.png'], Aircraft.res_dir);
    movefile([comp '-Side-View.fig'], Aircraft.res_dir);
    movefile([comp '-Side-View.png'], Aircraft.res_dir);
    movefile([comp '-Front-View.fig'], Aircraft.res_dir);
    movefile([comp '-Front-View.png'], Aircraft.res_dir);
    movefile([comp '3D.fig'], Aircraft.res_dir);
    movefile([comp '3D.png'], Aircraft.res_dir);
    
else
    disp (strcat(comp,' component does not exist!'))
    disp ('-----------------------------------------------------------------')
    disp ('-----------------------------------------------------------------')
    disp ([comp, 'Geometry utilities end run'])
    return 
    
    
end 