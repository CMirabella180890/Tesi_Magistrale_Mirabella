%TEMPLATE
function TemplateOut(Aircraft)

%% TEMPLATE Geometry Output function
% This function creates figures of 3-View of the fuselage and 3D Solid PLOT
% of the Aircraft based on here imported airfoils
%
% 1) Auxiliary variables are firstly defined
% 2) 2D plot -> 3-Views of aircraft drafting
% 3) 3D plot -> 3D PLOT of aircraft structure with imported airfoils

comp = 'TEMPLATE ';

if isfield(Aircraft.Geometry, "Template") == 1
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
    xtip_le = xle + ...
        (Aircraft.Geometry.Wing.b.value/2*...
        tan(Aircraft.Geometry.Wing.sweep.value/57.3))/lf;                       % wing tip leading edge %of fuselage lenght
    Aircraft.Geometry.Wing.ypos.value = 0.0;                                    % %of wing semispan
    %kink 1
    ckink1 = Aircraft.Geometry.Wing.chord_kink_one.value;                                   % kink 1 chord in meters
    ykink1 = Aircraft.Geometry.Wing.panel_span1.value * Aircraft.Geometry.Wing.b.value/2;   % kink 1 wing semispan in meters
    %kink 2
    ckink2 = Aircraft.Geometry.Wing.chord_kink_two.value;                                                                               % kink 2 chord in meters
    ykink2 = (Aircraft.Geometry.Wing.panel_span1.value + Aircraft.Geometry.Wing.panel_span2.value) * Aircraft.Geometry.Wing.b.value/2;  % kink 2 wing semispan in meters
    %sweep_le
    sweep_le = Aircraft.Geometry.Wing.sweep.value;
    % Importing airfoils sections for wing, horizontal, vertical...TO BE
    % MODIFIED FROM INPUT data structure!!!
    root_coord_w = importdata('_Airfoil\E216_104.txt'); % Root Airfoil coordiante read from file (x/c & z/c) IRON_Root.txt
    tip_coord_w = importdata('_Airfoil\E216_104.txt'); % Root Airfoil coordiante read from file (x/c & z/c) IRON_Tip.txt

    %% TOP-VIEW
    figure('Name',[comp '-Top-View'],'NumberTitle','off');
    
        grid on
    title([comp 'Top-View'])
    xlabel({'y (m)'})
    ylabel({'x (m)'})
    axis equal
    
    saveas(gcf, [comp '-Top-View.fig'])
    saveas(gcf, [comp '-Top-View.png'])   
    
    
    %% SIDE-VIEW
    figure('Name',[comp '-Side-View'],'NumberTitle','off');
    
    grid on
    title([comp 'Side-View'])
    xlabel({'y (m)'})
    ylabel({'x (m)'})
    axis equal
    
    saveas(gcf, [comp '-Side-View.fig'])
    saveas(gcf, [comp '-Side-View.png'])
    
    
    %% 
    figure('Name',[comp '-Front-View'],'NumberTitle','off');
    
    grid on
    title([comp 'Front-View'])
    xlabel({'y (m)'})
    ylabel({'x (m)'})
    axis equal
    
    saveas(gcf, [comp '-Front-View.fig'])
    saveas(gcf, [comp '-Front-View.png'])
    
    
    %% 3D
    figure('Name',[comp '3D'],'NumberTitle','off');

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