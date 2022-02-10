%TEMPLATE
function HorizontalOut(Aircraft)

%% TEMPLATE Geometry Output function
% This function creates figures of 3-View of the fuselage and 3D Solid PLOT
% of the Aircraft based on here imported airfoils
%
% 1) Auxiliary variables are firstly defined
% 2) 2D plot -> 3-Views of aircraft drafting
% 3) 3D plot -> 3D PLOT of aircraft structure with imported airfoils

comp = 'Horizontal';
cfg = 'lw_tt_bm';                           %aircraft configuration TO BE EDITED from outside

if isfield(Aircraft.Geometry, "Horizontal") == 1
    disp ([comp, 'exists'])
    disp ('-----------------------------------------------------------------')
    disp ('-----------------------------------------------------------------')
    disp ([comp, 'Geometry utilities end run'])
    
    %% Auxiliary declaration
    %fuselage
    df = Aircraft.Geometry.Fuselage.diameter.value;         % fuselage diameter (m)
    lf = Aircraft.Geometry.Fuselage.length.value;           % fuselage lenght (m)
    
    %horizontal
    hzpos = (Aircraft.Geometry.Horizontal.zloc.value - df/2*Aircraft.Geometry.Vertical.zpos.value  ) / Aircraft.Geometry.Vertical.b.value;   %horizontal root chord zeta position %of vtail span (poitive upward)
    Aircraft.Geometry.Horizontal.xle.value = Aircraft.Geometry.Horizontal.xloc.value/lf;      % % of lf
    xtip_le_h = Aircraft.Geometry.Horizontal.xle.value + ...
        Aircraft.Geometry.Horizontal.b.value/2*tan(Aircraft.Geometry.Horizontal.sweep.value/57.3)/lf;                                   % % of lf
    Aircraft.Geometry.Horizontal.xtip_le.value = xtip_le_h; % % of lf
    Aircraft.Geometry.Horizontal.ypos.value = Aircraft.Geometry.Horizontal.yloc.value/Aircraft.Geometry.Wing.b.value/2; % % of wing span
    
    %airfoils
    root_coord_h = importdata('_Airfoil\Simmetric_Root.txt'); % Root Airfoil coordiante read from file (x/c & z/c)
    tip_coord_h = importdata('_Airfoil\Simmetric_Tip.txt'); % Root Airfoil coordiante read from file (x/c & z/c)
    
    %% TOP-VIEW
    figure('Name',[comp '-Top-View'],'NumberTitle','off');
    hold on
    %horizontal
    plot([0 0], Aircraft.Geometry.Horizontal.xle.value*lf + ...
        [0, Aircraft.Geometry.Horizontal.croot.value],'g',"LineWidth",1.5)
    plot([Aircraft.Geometry.Horizontal.b.value/2, Aircraft.Geometry.Horizontal.b.value/2],...
        xtip_le_h*lf + ...
        [0, Aircraft.Geometry.Horizontal.ctip.value],'g',"LineWidth",1.5)
    plot([-Aircraft.Geometry.Horizontal.b.value/2, -Aircraft.Geometry.Horizontal.b.value/2],...
        xtip_le_h*lf + ...
        [0, Aircraft.Geometry.Horizontal.ctip.value],'g',"LineWidth",1.5)
    plot([0, Aircraft.Geometry.Horizontal.b.value/2],...
        [Aircraft.Geometry.Horizontal.xle.value*lf, xtip_le_h*lf],...
        'g',"LineWidth",1.5)
    plot([0, -Aircraft.Geometry.Horizontal.b.value/2],...
        [Aircraft.Geometry.Horizontal.xle.value*lf, xtip_le_h*lf],...
        'g',"LineWidth",1.5)
    plot([0, Aircraft.Geometry.Horizontal.b.value/2],...
        [Aircraft.Geometry.Horizontal.xle.value*lf+Aircraft.Geometry.Horizontal.croot.value,...
        Aircraft.Geometry.Horizontal.xtip_le.value*lf+Aircraft.Geometry.Horizontal.ctip.value],...
        'g',"LineWidth",1.5)
    plot([0, -Aircraft.Geometry.Horizontal.b.value/2],...
        [Aircraft.Geometry.Horizontal.xle.value*lf+Aircraft.Geometry.Horizontal.croot.value,...
        Aircraft.Geometry.Horizontal.xtip_le.value*lf+Aircraft.Geometry.Horizontal.ctip.value],...
        'g',"LineWidth",1.5)
    
    %ELEVATOR
    if isfield(Aircraft.Geometry, 'Elevator') == 1
        %linear chord law
        c0 = Aircraft.Geometry.Horizontal.croot.value;          % root chord (unit)
        c1 = Aircraft.Geometry.Horizontal.ctip.value;           % tip chord (unit)
        span0 = Aircraft.Geometry.Horizontal.yloc.value;
        span1 = Aircraft.Geometry.Horizontal.b.value/2;
        [c_0,c_y] = chord_linear(c0 ,c1, span0, span1);
        %
        y_inner = Aircraft.Geometry.Elevator.eta_inner.value * span1;
        y_outer = Aircraft.Geometry.Elevator.eta_outer.value * span1;
        %
        c_inner = c_0 - c_y * y_inner;                      %htail chord at inner station (unit)
        c_outer = c_0 - c_y * y_outer;                      %htail chord at inner station (unit)
        %
        ce_inner = c_inner * Aircraft.Geometry.Elevator.cf_c_inner.value;   %elevetor chord at inner station (unit)
        ce_outer = c_outer* Aircraft.Geometry.Elevator.cf_c_outer.value;    %elevetor chord at inner station (unit)
        %
        xle_inner = Aircraft.Geometry.Horizontal.xle.value*lf +  y_inner * tan(Aircraft.Geometry.Horizontal.sweep.value/57.3) +...
            c_inner - ce_inner;
        xle_outer = Aircraft.Geometry.Horizontal.xle.value*lf +  y_outer * tan(Aircraft.Geometry.Horizontal.sweep.value/57.3) +...
            c_outer - ce_outer;
        %
        plot([y_inner y_inner],...
            [xle_inner xle_inner+ce_inner],'g',"LineWidth",1.5)
        plot([y_outer y_outer],...
            [xle_outer xle_outer+ce_outer],'g',"LineWidth",1.5)
        plot([y_inner y_outer],...
            [xle_inner xle_outer],'g',"LineWidth",1.5)
        fill ([y_inner y_outer y_outer y_inner],[xle_inner xle_outer xle_outer+ce_outer xle_inner+ce_inner],'g')
        %
        plot([-y_inner -y_inner],...
            [xle_inner xle_inner+ce_inner],'g',"LineWidth",1.5)
        plot([-y_outer -y_outer],...
            [xle_outer xle_outer+ce_outer],'g',"LineWidth",1.5)
        plot([-y_inner -y_outer],...
            [xle_inner xle_outer],'g',"LineWidth",1.5)
        fill ([-y_inner -y_outer -y_outer -y_inner],[xle_inner xle_outer xle_outer+ce_outer xle_inner+ce_inner],'g')
        %
    end
    
    grid on
    title([comp 'Top-View'])
    xlabel({'y (m)'})
    ylabel({'x (m)'})
    axis equal
    
    saveas(gcf, [comp '-Top-View.fig'])
    saveas(gcf, [comp '-Top-View.png'])
    
    
    %% SIDE-VIEW
    figure('Name',[comp '-Side-View'],'NumberTitle','off');
    hold on
    
    %chord
    plot([Aircraft.Geometry.Horizontal.xle.value*lf,...
        Aircraft.Geometry.Horizontal.xle.value*lf+Aircraft.Geometry.Horizontal.croot.value]...
        ,[df/2+hzpos*Aircraft.Geometry.Vertical.b.value, ...
        df/2+hzpos*Aircraft.Geometry.Vertical.b.value],'g',"LineWidth",1.5)
    
    %airfoils
    root_coord = root_coord_h;
    tip_coord = tip_coord_h;
    
    %root airfoil
    root_coord_3D(:,1) = root_coord(:,1).*Aircraft.Geometry.Horizontal.croot.value +...
        Aircraft.Geometry.Horizontal.xle.value.*lf;
    
    root_coord_3D(:,2) = root_coord(:,1).*Aircraft.Geometry.Horizontal.ypos.value*Aircraft.Geometry.Wing.b.value;
    
    if cfg(4:5)=='ct'
        root_coord_3D(:,3) =  root_coord(:,2).*Aircraft.Geometry.Horizontal.croot.value +...
            hzpos*(Aircraft.Geometry.Vertical.b.value);
    else
        root_coord_3D(:,3) =  root_coord(:,2).*Aircraft.Geometry.Horizontal.croot.value +...
            df/2+hzpos*(Aircraft.Geometry.Vertical.b.value);
    end
    
    %tip airfoil
    tip_coord_3D(:,1) = tip_coord(:,1).*Aircraft.Geometry.Horizontal.ctip.value +...
        Aircraft.Geometry.Horizontal.xtip_le.value.*lf;
    
    tip_coord_3D(:,2) = Aircraft.Geometry.Horizontal.b.value/2;
    
    if cfg(4:5)=='ct'
        tip_coord_3D(:,3) =  tip_coord(:,2).*Aircraft.Geometry.Horizontal.ctip.value +...
            hzpos*(Aircraft.Geometry.Vertical.b.value)...
            +Aircraft.Geometry.Horizontal.b.value/2*tan(Aircraft.Geometry.Horizontal.dihedral.value/57.3);
    else
        tip_coord_3D(:,3) =  tip_coord(:,2).*Aircraft.Geometry.Horizontal.croot.value +...
            df/2+hzpos*(Aircraft.Geometry.Vertical.b.value);
    end
    
    plot(root_coord_3D(:,1),root_coord_3D(:,3),'g',"LineWidth",2)
    plot(tip_coord_3D(:,1),tip_coord_3D(:,3),'g',"LineWidth",2)
    
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
    
    %horizontal
    if cfg(4:5) == 'ct'
        plot([0, Aircraft.Geometry.Horizontal.b.value/2], [hzpos*Aircraft.Geometry.Vertical.b.value, ...
            hzpos*Aircraft.Geometry.Vertical.b.value+...
            Aircraft.Geometry.Horizontal.b.value/2*tan(Aircraft.Geometry.Horizontal.dihedral.value/57.3)],'g',"LineWidth",1.5);
        plot([0, -Aircraft.Geometry.Horizontal.b.value/2], [hzpos*Aircraft.Geometry.Vertical.b.value, ...
            hzpos*Aircraft.Geometry.Vertical.b.value+...
            Aircraft.Geometry.Horizontal.b.value/2*tan(Aircraft.Geometry.Horizontal.dihedral.value/57.3)],'g',"LineWidth",1.5);
    else
        plot([0, Aircraft.Geometry.Horizontal.b.value/2], [df/2+...
            hzpos*Aircraft.Geometry.Vertical.b.value, ...
            df/2+hzpos*Aircraft.Geometry.Vertical.b.value+...
            Aircraft.Geometry.Horizontal.b.value/2*tan(Aircraft.Geometry.Horizontal.dihedral.value/57.3)],'g',"LineWidth",1.5);
        plot([0, -Aircraft.Geometry.Horizontal.b.value/2], [df/2+...
            hzpos*Aircraft.Geometry.Vertical.b.value, ...
            df/2+hzpos*Aircraft.Geometry.Vertical.b.value+...
            Aircraft.Geometry.Horizontal.b.value/2*tan(Aircraft.Geometry.Horizontal.dihedral.value/57.3)],'g',"LineWidth",1.5);
    end
    
    grid on
    title([comp 'Front-View'])
    xlabel({'y (m)'})
    ylabel({'x (m)'})
    axis equal
    
    saveas(gcf, [comp '-Front-View.fig'])
    saveas(gcf, [comp '-Front-View.png'])
    
    
    %% 3D
    figure('Name',[comp '3D'],'NumberTitle','off');
    hold on

    XXh = [root_coord_3D(:,1), tip_coord_3D(:,1)]';
    YYh = [root_coord_3D(:,2), tip_coord_3D(:,2)]';
    ZZh = [root_coord_3D(:,3), tip_coord_3D(:,3)]';
    
    surf(XXh,YYh,ZZh,'FaceColor','green','EdgeColor','none')
    surf(XXh,-YYh,ZZh,'FaceColor','gree','EdgeColor','none')
    daspect([1.0 , 1.0, 1.0])
    camlight left; lighting phong

    grid on
    title([comp '3D'])
    xlabel({'y (m)'})
    ylabel({'x (m)'})
    axis equal
    
    saveas(gcf, [comp '3D.fig'])
    saveas(gcf, [comp '3D.png'])
    %%
    
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
    
    
    
end