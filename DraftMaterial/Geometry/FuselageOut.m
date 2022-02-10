function FuselageOut(Aircraft)
%% Fuselage Geometry Output function
% This function creates figures of 3-View of the fuselage and 3D Solid PLOT
% of the Aircraft based on here imported airfoils
%
% 1) Auxiliary variables are firstly defined
% 2) 2D plot -> 3-Views of aircraft drafting
% 3) 3D plot -> 3D PLOT of aircraft structure with imported airfoils

comp = 'Fuselage';

if isfield(Aircraft.Geometry, "Fuselage") == 1
    disp ([comp, 'exists'])
    disp ('-----------------------------------------------------------------')
    disp ('-----------------------------------------------------------------')
    disp ('Fuselage Geometry utilities running')
    
    %% Auxiliary declaration
    
    df = Aircraft.Geometry.Fuselage.diameter.value;         % fuselage diameter (unit)
    dbase = 0.1*df;                                         %fuselage base diameter (unit)
    lf = Aircraft.Geometry.Fuselage.length.value;           % fuselage lenght (unit)
    Aircraft.Geometry.Fuselage.kcockpit.value = 1.0;        % nose ratio = lnose/df
    Aircraft.Geometry.Fuselage.ktail.value = 1.5;           % tail ratio = ltail/df

    %% TOP VIEW

    figure('Name','Fuselage-Top-View','NumberTitle','off');
    
    %TOP-VIEW
    % cabin
    plot([df/2, df/2],...
        [0+Aircraft.Geometry.Fuselage.kcockpit.value*df, ...
        lf-Aircraft.Geometry.Fuselage.ktail.value*df],...
        'b',"LineWidth",1.5)
    hold on
    plot([-df/2, -df/2],...
        [Aircraft.Geometry.Fuselage.kcockpit.value*df, ...
        lf-Aircraft.Geometry.Fuselage.ktail.value*df],...
        'b',"LineWidth",1.5)
    
    %cockpit - parabolic reconstruction
    x = linspace (-df/2, df/2, 50);
    a = Aircraft.Geometry.Fuselage.kcockpit.value * df / ((df/2)^2);
    plot(x,a.*x.^2,'b',"LineWidth",1.5)
    %
    %tail - linear construction
    plot([-dbase/2, dbase/2],...
        [lf, lf],'b',"LineWidth",1.5)
    plot([dbase/2, df/2],...
        [lf,...
        lf-Aircraft.Geometry.Fuselage.ktail.value*df],...
        'b',"LineWidth",1.5)
    plot([-dbase/2, -df/2],...
        [lf,...
        lf-Aircraft.Geometry.Fuselage.ktail.value*df],...
        'b',"LineWidth",1.5)
    hold off
    
    grid on
    title('Fuselage Top-View')
    xlabel({'y (m)'})
    ylabel({'x (m)'})
    axis equal
    
    saveas(gcf, 'Fuselage-Top-View.fig')
    saveas(gcf, 'Fuselage-Top-View.png')
    
    %SIDE-VIEW
    %
    figure('Name','Fuselage-Side-View','NumberTitle','off');
    
    %cockpit - parabolic construction
    plot(a.*x.^2,x,'b',"LineWidth",1.5)
    hold on
    %dbase
    plot([lf, lf],...
        [df/2, df/2-dbase],...
        'b',"LineWidth",1.5)
    %cabin
    plot([0+Aircraft.Geometry.Fuselage.kcockpit.value*df,...
        lf-Aircraft.Geometry.Fuselage.ktail.value*df],...
        [df/2, df/2],'b',"LineWidth",1.5)
    plot([0+Aircraft.Geometry.Fuselage.kcockpit.value*df, ...
        lf-Aircraft.Geometry.Fuselage.ktail.value*df],...
        [-df/2, -df/2],'b',"LineWidth",1.5)
    %tail
    plot([lf,...
        lf-Aircraft.Geometry.Fuselage.ktail.value*df],...
        [df/2, df/2],'b',"LineWidth",1.5)
    plot([lf,...
        lf-Aircraft.Geometry.Fuselage.ktail.value*df],...
        [df/2-dbase, -df/2],'b',"LineWidth",1.5)
    
    hold off
    grid on
    title('Fuselage Side-View')
    xlabel({'y (m)'})
    ylabel({'x (m)'})
    axis equal
    
    saveas(gcf, 'Fuselage-Side-View.fig')
    saveas(gcf, 'Fuselage-Side-View.png')
    
    %FRONT-VIEW
    figure('Name','Front-View','NumberTitle','off');
    teta = linspace (0,360,100);
    y = df/2 .* cos(teta/57.3);
    z = df/2 .* sin(teta/57.3);
    plot(y,z,'b',"LineWidth",1.5)
    
    grid on
    title('Fuselage Front-View')
    xlabel({'y (m)'})
    ylabel({'x (m)'})
    axis equal
    
    saveas(gcf, 'Fuselage-Front-View.fig')
    saveas(gcf, 'Fuselage-Front-View.png')
    
    %3D
        figure('Name','Fuselage-3D','NumberTitle','off');
    nsection = 30;                  % Number of fuselage sections
    nsecockpit = 20;                % Number of fuselage cockptit sections
    nsecabin = 1;                   % Number of fuselage cabin sections
    nsectail = nsection-nsecockpit-nsecabin; % Number of fuselage tail sections
    
    step_cabin = (lf-(Aircraft.Geometry.Fuselage.kcockpit.value *...
        df))/nsecabin;
    
    
    points = 100;
    teta = linspace (0,360,points);    % circumference sections
    
    
    Xf = ones (nsection,points);
    Yf = ones (nsection,points);
    Zf = ones (nsection,points);
    
    for i=1 : nsection
        
        %cockpit
        if i<=nsecockpit
            secdiameter=linspace (0,df, nsecockpit);
            y = secdiameter(i)/2* cos(teta/57.3);
            z = secdiameter(i)/2* sin(teta/57.3);
            
            %linear reconstruction
            % step_cockpit =
            % Aircraft.Geometry.Fuselage.kcockpit.value*df))/nsecockpit
            %X(i,1:end) = 0+(i-1)*step_cockpit;
            
            %parabolic reconstruction
            x_cockpit = linspace (0, df/2, nsecockpit);
            a = Aircraft.Geometry.Fuselage.kcockpit.value * df / ((df/2)^2);
            ycockpit = a*x_cockpit(i)^2;
            
            Xf(i,1:end) = ycockpit;
            Yf(i,1:end) = y;
            Zf(i,1:end) = z;
        end
        
        %cabin
        if i>nsecockpit && i<=(nsecockpit+nsecabin)
            %secdiameter = df;
            Xf(i,1:end) = (lf-...
                (Aircraft.Geometry.Fuselage.ktail.value *...
                df));
            Yf(i,1:end) = y;
            Zf(i,1:end) = z;
        end
        
        %tail
        if i>(nsecockpit+nsecabin)
            secdiameter = linspace (df, ...
                dbase,nsectail+1);
            y = secdiameter(i+1-(nsecockpit+nsecabin))/2* cos(teta/57.3);
            z = secdiameter(i+1-(nsecockpit+nsecabin))/2* sin(teta/57.3);
            
            %linear reconstruction
            step_tail = (Aircraft.Geometry.Fuselage.ktail.value*df)/nsectail;
            Xf(i,1:end) = (lf -...
                Aircraft.Geometry.Fuselage.ktail.value*df)...
                +(i-(nsecockpit+nsecabin))*step_tail;
            
            Yf(i,1:end) = y;
            tailcone_slope = ((df/2-dbase/2)/...
                (Aircraft.Geometry.Fuselage.ktail.value*df));
            Zf(i,1:end) = z +...
                ((i-(nsecockpit+nsecabin))*step_tail)*tailcone_slope;
            
            %Z(i,1:end) = z;
        end
        
    end
    
    surf(Xf,Yf,Zf,'FaceColor','blu','EdgeColor','none')
    camlight left; lighting phong
    daspect([1.0 , 1.0, 1.0])
    
    grid on
    title('Fuselage Side-View')
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
    
end
end