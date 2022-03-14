%TEMPLATE
function EngineOut(Aircraft)

%% TEMPLATE Geometry Output function
% This function creates figures of 3-View of the fuselage and 3D Solid PLOT
% of the Aircraft based on here imported airfoils
%
% 1) Auxiliary variables are firstly defined
% 2) 2D plot -> 3-Views of aircraft drafting
% 3) 3D plot -> 3D PLOT of aircraft structure with imported airfoils

comp = 'Engine';

if isfield(Aircraft.Geometry, "Engine") == 1
    disp ([comp, 'exists'])
    disp ('-----------------------------------------------------------------')
    disp ('-----------------------------------------------------------------')
    disp ([comp, 'Geometry utilities end run'])
    
    %% Auxiliary declaration
    %Tlar and Configurations
    cfg = 'lw_tt_bm';                           %aircraft configuration TO BE EDITED from outside
    Aircraft.TLAR.Type_engine.value = "TP";     %aircraft configuration TO BE EDITED from outside
    Aircraft.Geometry.Engine.prop_config.value = 'tractor'; % 'tractor' or 'pusher' to be provided from TLARs
    prop_config = Aircraft.Geometry.Engine.prop_config.value;                     
    
    %fuselage
    df = Aircraft.Geometry.Fuselage.diameter.value;         % fuselage diameter (m)
    lf = Aircraft.Geometry.Fuselage.length.value;           % fuselage lenght (m)
    %wing
    zpos = Aircraft.Geometry.Wing.zle.value;                                    % wing root chord zeta position (unit)
    xle = Aircraft.Geometry.Wing.xle.value/lf;                                  % wing tip leading edge %of fuselage lenght
    xtip_le = xle + ...
        (Aircraft.Geometry.Wing.b.value/2*...
        tan(Aircraft.Geometry.Wing.sweep_second.value/57.3))/lf;                       % wing tip leading edge %of fuselage lenght
    Aircraft.Geometry.Wing.ypos.value = 0.0;                                    % %of wing semispan
% %     %engine
% %     Aircraft.Geometry.Engine.Primary.xpos.value = 0.90;             % % lf
% %     Aircraft.Geometry.Engine.Primary.lf.value = 0.5;                % engine lenght m
% %     Aircraft.Geometry.Engine.Primary.ypos.value = 0.0;              % %of wing semispan
% %     Aircraft.Geometry.Engine.Primary.df.value = 0.1;                % m
% %     Aircraft.Geometry.Engine.Primary.propdiam.value = 0.6;          %prop diameter in meters
% %     Aircraft.Geometry.Engine.Primary.zpos.value = 1.0;              %of fus df
    
%     Aircraft.Geometry.Engine.Primary.xpos.value = 0.0;             % % lf
%     Aircraft.Geometry.Engine.Primary.lf.value = 1.2;                % engine lenght m
%     Aircraft.Geometry.Engine.Primary.ypos.value = 0.0;              % %of wing semispan
%     Aircraft.Geometry.Engine.Primary.df.value = 0.8;                % m
%     Aircraft.Geometry.Engine.Primary.propdiam.value = 1.2;          %prop diameter in meters
%     Aircraft.Geometry.Engine.Primary.zpos.value = 0.0;              %of fus df
    
 ezpos = Aircraft.Geometry.Engine.Primary.zpos.value;       % % df engine zeta position

    
    %% TOP-VIEW
    figure('Name',[comp '-Top-View'],'NumberTitle','off');
    hold on
    %engine
    if cfg(end-1:end)== 'bm'
        xmin = Aircraft.Geometry.Engine.Primary.xpos.value * lf;
        xmax = Aircraft.Geometry.Engine.Primary.xpos.value * lf + Aircraft.Geometry.Engine.Primary.lf.value;
        yengine = Aircraft.Geometry.Engine.Primary.ypos.value * Aircraft.Geometry.Wing.b.value/2;
        yinner = yengine - Aircraft.Geometry.Engine.Primary.df.value/2;
        youter = yengine + Aircraft.Geometry.Engine.Primary.df.value/2;
    else
        xmin = Aircraft.Geometry.Engine.Primary.xpos.value * lf;
        xmax = Aircraft.Geometry.Engine.Primary.xpos.value * lf + Aircraft.Geometry.Engine.Primary.lf.value;
        yengine = Aircraft.Geometry.Engine.Primary.ypos.value * Aircraft.Geometry.Wing.b.value/2;
        yinner = yengine - Aircraft.Geometry.Engine.Primary.df.value/2;
        youter = yengine + Aircraft.Geometry.Engine.Primary.df.value/2;
    end
    
    plot([yengine yengine],[xmin xmax],'k--',"LineWidth",1)
    plot([-yengine -yengine],[xmin xmax],'k--',"LineWidth",1)
    plot([yinner yinner],[xmin xmax],'k',"LineWidth",2)
    plot([-yinner -yinner],[xmin xmax],'k',"LineWidth",2)
    plot([youter youter],[xmin xmax],'k',"LineWidth",2)
    plot([-youter -youter],[xmin xmax],'k',"LineWidth",2)
    plot([yinner youter],[xmin xmin],'k',"LineWidth",2)
    plot([yinner youter],[xmax xmax],'k',"LineWidth",2)
    plot([-yinner -youter],[xmin xmin],'k',"LineWidth",2)
    plot([-yinner -youter],[xmax xmax],'k',"LineWidth",2)
    
    % propeller (only for 'TP' aircraft)
    if  Aircraft.TLAR.Type_engine.value == ("TP")
        if strcmp(prop_config,'tractor')==1
            if cfg(end-1:end) == 'bm'
                dprop = Aircraft.Geometry.Engine.Primary.propdiam.value;    %prop diameter in meters
                plot([yengine-dprop/2 yengine+dprop/2],[xmin xmin],'k',"LineWidth",1)
                plot([-yengine-dprop/2 -yengine+dprop/2],[xmin xmin],'k',"LineWidth",1)
            else
                dprop = Aircraft.Geometry.Engine.Primary.propdiam.value;    %prop diameter in meters
                plot([yengine-dprop/2 yengine+dprop/2],[xmin xmin],'k',"LineWidth",1)
                plot([-yengine-dprop/2 -yengine+dprop/2],[xmin xmin],'k',"LineWidth",1)
            end
        else
            if cfg(end-1:end) == 'bm'
                dprop = Aircraft.Geometry.Engine.Primary.propdiam.value;    %prop diameter in meters
                plot([yengine-dprop/2 yengine+dprop/2],[xmax xmax],'k',"LineWidth",1)
                plot([-yengine-dprop/2 -yengine+dprop/2],[xmax xmax],'k',"LineWidth",1)
            else
                dprop = Aircraft.Geometry.Engine.Primary.propdiam.value;    %prop diameter in meters
                plot([yengine-dprop/2 yengine+dprop/2],[xmax xmax],'k',"LineWidth",1)
                plot([-yengine-dprop/2 -yengine+dprop/2],[xmax xmax],'k',"LineWidth",1)
            end
            
        end
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
    %ENGINE
    zmin = ezpos*df/2 - Aircraft.Geometry.Engine.Primary.df.value/2;
    zmax = ezpos*df/2 + Aircraft.Geometry.Engine.Primary.df.value/2;
    
    %PROPELLER
    if Aircraft.TLAR.Type_engine.value == "TP"
        if strcmp(prop_config,'tractor')==1
            if cfg(end-1:end) == "wm"
                zeta_engine = Aircraft.Geometry.Engine.Primary.zpos.value*df/2;
                plot ([xmin xmin],[zeta_engine+dprop/2 zeta_engine-dprop/2],'k',"LineWidth",1)
            else
                zeta_engine = Aircraft.Geometry.Engine.Primary.zpos.value*df/2;
                plot ([xmin xmin],[zeta_engine+dprop/2 zeta_engine-dprop/2],'k',"LineWidth",1)
                %         zmin = zeta_engine-dprop/2;
                %         zmax = zeta_engine+dprop/2;
            end
        else
            if cfg(end-1:end) == "wm"
                zeta_engine = Aircraft.Geometry.Engine.Primary.zpos.value*df/2;
                plot ([xmax xmax],[zeta_engine+dprop/2 zeta_engine-dprop/2],'k',"LineWidth",1)
            else
                zeta_engine = Aircraft.Geometry.Engine.Primary.zpos.value*df/2;
                plot ([xmax xmax],[zeta_engine+dprop/2 zeta_engine-dprop/2],'k',"LineWidth",1)
                %         zmin = zeta_engine-dprop/2;
                %         zmax = zeta_engine+dprop/2;
            end
        end
    end
    
    plot([xmin xmax],[zmin zmin],'k',"LineWidth",2)
    plot([xmin xmax],[zmax zmax],'k',"LineWidth",2)
    plot([xmin xmin],[zmin zmax],'k',"LineWidth",2)
    plot([xmax xmax],[zmin zmax],'k',"LineWidth",2)
    
    grid on
    title([comp 'Side-View'])
    xlabel({'y (m)'})
    ylabel({'x (m)'})
    axis equal
    
    saveas(gcf, [comp '-Side-View.fig'])
    saveas(gcf, [comp '-Side-View.png'])
    
    
    %%
    figure('Name',[comp '-Front-View'],'NumberTitle','off');
    hold on
    %engine
    teta = linspace (0,360,100);
    y = Aircraft.Geometry.Engine.Primary.df.value/2 .* cos(teta/57.3);
    z = Aircraft.Geometry.Engine.Primary.df.value/2 .* sin(teta/57.3);
    plot(y+yengine, z+Aircraft.Geometry.Engine.Primary.zpos.value*df/2,'k',"LineWidth",1.5)
    plot(y-yengine, z+Aircraft.Geometry.Engine.Primary.zpos.value*df/2,'k',"LineWidth",1.5)
    
    %propeller
    if Aircraft.TLAR.Type_engine.value == "TP"
        zeta_engine = Aircraft.Geometry.Engine.Primary.zpos.value*df/2;
        plot ([yengine-dprop/2 yengine+dprop/2], [zeta_engine zeta_engine],'k',"LineWidth",1)
        plot ([-yengine-dprop/2 -yengine+dprop/2], [zeta_engine zeta_engine],'k',"LineWidth",1)
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
    points = 100;

    %ENGINE
    if isfield(Aircraft.Geometry.Engine,'Secondary') && Aircraft.TLAR.Number_engine.value==2
        phi=[Aircraft.Powertrain.Secondary.Power.Alternate.phi_alternate.value,...
            Aircraft.Powertrain.Secondary.Power.Climb.phi_climb.value,...
            Aircraft.Powertrain.Secondary.Power.Cruise.phi_cruise.value,...
            Aircraft.Powertrain.Secondary.Power.Descent.phi_descent.value,...
            Aircraft.Powertrain.Secondary.Power.Landing.phi_landing.value,...
            Aircraft.Powertrain.Secondary.Power.Loiter.phi_loiter.value,...
            Aircraft.Powertrain.Secondary.Power.Takeoff.phi_takeoff.value];
        if any(phi)
            
            nenginesection = 2; % Number of engine sections;
            
            Xe = ones (nenginesection,points);
            Ye = ones (nenginesection,points);
            Ze = ones (nenginesection,points);
            
            for i=1:nenginesection
                
                if i == 1
                    Xe(i,1:end) = xmin;
                end
                
                if i == nenginesection
                    Xe(i,1:end) = xmax;
                end
                secdiameter = Aircraft.Geometry.Engine.Primary.df.value;
                y = secdiameter / 2* cos(teta/57.3);
                z = secdiameter / 2* sin(teta/57.3);
                Ye(i,1:end) = y+Aircraft.Geometry.Engine.Primary.ypos.value*...
                    Aircraft.Geometry.Wing.b.value/2;
                Ze(i,1:end) = z+Aircraft.Geometry.Engine.Primary.zpos.value*df/2;
                
            end
        else
            nenginesection = 2; % Number of engine sections;
            
            Xe = ones (nenginesection,points);
            Ye = ones (nenginesection,points);
            Ze = ones (nenginesection,points);
            
            for i=1:nenginesection
                
                if i == 1
                    Xe(i,1:end) = xmin;
                end
                
                if i == nenginesection
                    Xe(i,1:end) = xmax;
                end
                secdiameter = Aircraft.Geometry.Engine.Primary.df.value;
                y = secdiameter / 2* cos(teta/57.3);
                z = secdiameter / 2* sin(teta/57.3);
                Ye(i,1:end) = y+yengine;
                Ze(i,1:end) = z+Aircraft.Geometry.Engine.Primary.zpos.value*df/2;
                
            end
        end
    else
        
        nenginesection = 2; % Number of engine sections;
        
        Xe = ones (nenginesection,points);
        Ye = ones (nenginesection,points);
        Ze = ones (nenginesection,points);
        
        for i=1:nenginesection
            
            if i == 1
                Xe(i,1:end) = xmin;
            end
            
            if i == nenginesection
                Xe(i,1:end) = xmax;
            end
            secdiameter = Aircraft.Geometry.Engine.Primary.df.value;
            y = secdiameter / 2* cos(teta/57.3);
            z = secdiameter / 2* sin(teta/57.3);
            Ye(i,1:end) = y+yengine;
            Ze(i,1:end) = z+Aircraft.Geometry.Engine.Primary.zpos.value*df/2;
            
        end
    end
    surf(Xe,Ye,Ze,'FaceColor','black','EdgeColor','none')
    surf(Xe,-Ye,Ze,'FaceColor','black','EdgeColor','none')

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
    return
    
    
end 