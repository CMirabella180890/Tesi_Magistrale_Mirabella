function AircraftOut(Aircraft)
%%Predesign Output function
% This function creates figures of 3-View of the Aircraft and 3D Solid PLOT
% of the Aircraft based on here imported airfoils
%
% 1) Auxiliary variables are firstly defined 
% 2) 2D plot -> 3-Views of aircraft drafting
% 3) 3D plot -> 3D PLOT of aircraft structure with imported airfoils

%% Auxiliary variables
%Tlar and Configurations
cfg = 'lw_tt_bm';     %aircraft configuration
Aircraft.TLAR.Type_engine.value = "TP";

% zpos = Aircraft.Geometry.Wing.zpos.value;          %wing root chord zeta position 
% hzpos = Aircraft.Geometry.Horizontal.zpos.value;   %horizontal root chord zeta position
% ezpos = Aircraft.Geometry.Engine.Primary.zpos.value;       %engine zeta position
% cfg = char(Aircraft.TLAR.Configuration.value);     %aircraft configuration
Aircraft.Geometry.Wing.zpos.value = -0.15;          %wing root chord zeta position
zpos = Aircraft.Geometry.Wing.zpos.value;          %wing root chord zeta position 
hzpos = 0.0;   %horizontal root chord zeta position

%fuselage
% df = Aircraft.Geometry.Fuselage.df.value;           % fuselage diameter (m)
% dbase = Aircraft.Geometry.Fuselage.dbase.value;     %fuselage base diameter
% lf = Aircraft.Geometry.Fuselage.lf.value;           % fuselage lenght (m)
df = 0.35;           % fuselage diameter (m)
dbase = 0.03;     %fuselage base diameter
lf = 5.0;           % fuselage lenght (m)
Aircraft.Geometry.Fuselage.kcockpit.value = 1.0;
Aircraft.Geometry.Fuselage.ktail.value = 1.5;

%wing
Aircraft.Geometry.Wing.xtip_le.value =0.5;     % wing tip leading edge %of fuselage lenght
xtip_le = Aircraft.Geometry.Wing.xtip_le.value;     % wing tip leading edge %of fuselage lenght
Aircraft.Geometry.Wing.xle.value = 0.5;         % wing root leading edge %of fuselage lenght
Aircraft.Geometry.Wing.croot.value = 0.5;       %meter
Aircraft.Geometry.Wing.ypos.value = 0.0; % % of wing semispan

%horizontal
%xtip_le = Aircraft.Geometry.Horizontal.xtip_le.value; % horizontal tip leading edge %of fuselage lenght
xtip_le_h = 0.95;
Aircraft.Geometry.Horizontal.xle.value = 0.92;
Aircraft.Geometry.Horizontal.xtip_le.value = xtip_le_h;
Aircraft.Geometry.Horizontal.b.value=1.496; 
Aircraft.Geometry.Horizontal.ypos.value = 0.0; % % of wing span

%vertical
Aircraft.Geometry.Vertical.xle.value = 0.95; %of fuselage lenght
Aircraft.Geometry.Vertical.croot.value = 0.3136; %m
Aircraft.Geometry.Vertical.ctip.value = 0.1534725; %m
Aircraft.Geometry.Vertical.xtip_le.value = 1.0; %of fuselage lenght
xtip_le_v = Aircraft.Geometry.Vertical.xtip_le.value;
Aircraft.Geometry.Vertical.b.value = 0.437502; %m
Aircraft.Geometry.Vertical.zpos.value = 1.0; % % of df

%engine 
ezpos = 1.0;       % % df engine zeta position
Aircraft.Geometry.Engine.Primary.xpos.value = 0.90;
Aircraft.Geometry.Engine.Primary.lf.value = 0.5; %m
Aircraft.Geometry.Engine.Primary.ypos.value = 0.0; % %of wing semispan
Aircraft.Geometry.Engine.Primary.df.value = 0.1; %m
Aircraft.Geometry.Engine.Primary.propdiam.value = 0.4;    %prop diameter in meters
Aircraft.Geometry.Engine.Primary.zpos.value = 1.0;  %of fus df
%landing gear
% x_main_lg = Aircraft.Geometry.Undercarriage.Main.xpos.value*lf;
% x_nose_lg = Aircraft.Geometry.Undercarriage.Nose.xpos.value*lf;
% y_main_lg = Aircraft.Geometry.Undercarriage.Main.ypos.value*df;
% y_nose_lg = Aircraft.Geometry.Undercarriage.Nose.ypos.value*df;
% z_main_lg = Aircraft.Geometry.Undercarriage.Main.zpos.value*df;
% z_nose_lg = Aircraft.Geometry.Undercarriage.Nose.zpos.value*df;
% d_wheel_main = Aircraft.Geometry.Undercarriage.Main.diameter.value;
% d_wheel_nose = Aircraft.Geometry.Undercarriage.Nose.diameter.value;
x_main_lg = 0.6*lf;
x_nose_lg = 0.1*lf;
y_main_lg = 1.1*df;
y_nose_lg = 0.0*df;
z_main_lg = -1.1*df;
z_nose_lg = -1.1*df;
d_wheel_main = 0.1;
d_wheel_nose = 0.1;
Aircraft.Geometry.Undercarriage.Main.diameter.value = 0.1; %m
Aircraft.Geometry.Undercarriage.Nose.diameter.value = 0.1; %m

wheel_width = 0.1*df; % main gear wheel width (m)

% Importing airfoils sections for wing, horizontal, vertical...
root_coord_w = importdata('_Airfoil\IRON_Root.txt'); % Root Airfoil coordiante read from file (x/c & z/c)
tip_coord_w = importdata('_Airfoil\IRON_Tip.txt'); % Root Airfoil coordiante read from file (x/c & z/c)

root_coord_h = importdata('_Airfoil\Simmetric_Root.txt'); % Root Airfoil coordiante read from file (x/c & z/c)
tip_coord_h = importdata('_Airfoil\Simmetric_Tip.txt'); % Root Airfoil coordiante read from file (x/c & z/c)

root_coord_v = importdata('_Airfoil\Simmetric_Root.txt'); % Root Airfoil coordiante read from file (x/c & z/c)
tip_coord_v = importdata('_Airfoil\Simmetric_Tip.txt'); % Root Airfoil coordiante read from file (x/c & z/c)


%% 2D Aicraft 3-views
% Planform view

figure('Name','Aircraft 3-View','NumberTitle','off');
subplot(2,2,1)

%fuselage
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
%
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


%wing 
plot([0 0],Aircraft.Geometry.Wing.xle.value*lf + ...
    [0, Aircraft.Geometry.Wing.croot.value],'r',"LineWidth",1.5)
plot([Aircraft.Geometry.Wing.b.value/2, Aircraft.Geometry.Wing.b.value/2],...
    xtip_le*lf +...
    [0, Aircraft.Geometry.Wing.ctip.value],'r',"LineWidth",1.5)
plot([-Aircraft.Geometry.Wing.b.value/2, -Aircraft.Geometry.Wing.b.value/2],...
    xtip_le*lf + ...
    [0, Aircraft.Geometry.Wing.ctip.value],'r',"LineWidth",1.5)
plot([0, Aircraft.Geometry.Wing.b.value/2],...
    [Aircraft.Geometry.Wing.xle.value*lf, xtip_le*lf],...
    'r',"LineWidth",1.5)
plot([0 -Aircraft.Geometry.Wing.b.value/2],...
    [Aircraft.Geometry.Wing.xle.value*lf, xtip_le*lf],...
    'r',"LineWidth",1.5)
plot([0, Aircraft.Geometry.Wing.b.value/2],...
    [Aircraft.Geometry.Wing.xle.value*lf+Aircraft.Geometry.Wing.croot.value, xtip_le*lf+Aircraft.Geometry.Wing.ctip.value],...
    'r',"LineWidth",1.5)
plot([0, -Aircraft.Geometry.Wing.b.value/2],...
    [Aircraft.Geometry.Wing.xle.value*lf+Aircraft.Geometry.Wing.croot.value, xtip_le*lf+Aircraft.Geometry.Wing.ctip.value],...
    'r',"LineWidth",1.5)

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

%vertical
plot([0, 0], Aircraft.Geometry.Vertical.xle.value*lf + ...
    [0, Aircraft.Geometry.Vertical.croot.value],'y',"LineWidth",2)

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
    if cfg(end-1:end)== 'bm'
        dprop = Aircraft.Geometry.Engine.Primary.propdiam.value;    %prop diameter in meters
        plot([yengine-dprop/2 yengine+dprop/2],[xmin xmin],'k',"LineWidth",1)
        plot([-yengine-dprop/2 -yengine+dprop/2],[xmin xmin],'k',"LineWidth",1)
    else
        dprop = Aircraft.Geometry.Engine.Primary.propdiam.value;    %prop diameter in meters
        plot([yengine-dprop/2 yengine+dprop/2],[xmin xmin],'k',"LineWidth",1)
        plot([-yengine-dprop/2 -yengine+dprop/2],[xmin xmin],'k',"LineWidth",1)
    end
end

grid on
title('Aircraft Top-View')
xlabel({'y (m)'})
ylabel({'x (m)'})
axis equal


% Side view
subplot(2,2,4)

%fuselage
%
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

% %wing
% plot([Aircraft.Geometry.Wing.xle.value*lf,...
%     Aircraft.Geometry.Wing.xle.value*lf+Aircraft.Geometry.Wing.croot.value]...
%     ,[zpos*df/2, zpos*df/2],'r',"LineWidth",1.5)

% %horizontal
% plot([Aircraft.Geometry.Horizontal.xle.value*lf,...
%     Aircraft.Geometry.Horizontal.xle.value*lf+Aircraft.Geometry.Horizontal.croot.value]...
%     ,[df/2+hzpos*Aircraft.Geometry.Vertical.b.value, ...
%     df/2+hzpos*Aircraft.Geometry.Vertical.b.value],'g',"LineWidth",1.5)

%vertical

plot([Aircraft.Geometry.Vertical.xle.value*lf, Aircraft.Geometry.Vertical.xle.value*lf+Aircraft.Geometry.Vertical.croot.value],...
    [df/2, df/2],'y',"LineWidth",1.5)
plot([Aircraft.Geometry.Vertical.xtip_le.value*lf, Aircraft.Geometry.Vertical.xtip_le.value*lf+Aircraft.Geometry.Vertical.ctip.value],...
    [df/2+Aircraft.Geometry.Vertical.b.value, df/2+Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)
plot([Aircraft.Geometry.Vertical.xle.value*lf, xtip_le_v*lf],...
    [df/2, df/2+Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)

plot([Aircraft.Geometry.Vertical.xle.value*lf+Aircraft.Geometry.Vertical.croot.value, xtip_le_v*lf+Aircraft.Geometry.Vertical.ctip.value],...
    [df/2, df/2+Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)

%engine
zmin = ezpos*df/2 - Aircraft.Geometry.Engine.Primary.df.value/2;
zmax = ezpos*df/2 + Aircraft.Geometry.Engine.Primary.df.value/2;


%propeller
if Aircraft.TLAR.Type_engine.value == "TP"
    if cfg(end-1:end) == "wm"
        zeta_engine = Aircraft.Geometry.Engine.Primary.zpos.value*df/2;
        plot ([xmin xmin],[zeta_engine+dprop/2 zeta_engine-dprop/2],'k',"LineWidth",1)
    else
        zeta_engine = Aircraft.Geometry.Engine.Primary.zpos.value*df/2;
        plot ([xmin xmin],[zeta_engine+dprop/2 zeta_engine-dprop/2],'k',"LineWidth",1)
%         zmin = zeta_engine-dprop/2;
%         zmax = zeta_engine+dprop/2;
    end
    
end

plot([xmin xmax],[zmin zmin],'k',"LineWidth",2)
plot([xmin xmax],[zmax zmax],'k',"LineWidth",2)
plot([xmin xmin],[zmin zmax],'k',"LineWidth",2)
plot([xmax xmax],[zmin zmax],'k',"LineWidth",2)


%undercarriage
teta = linspace (0,360,100);

x = x_main_lg  + ...
    d_wheel_main/2.* cos(teta/57.3);
z = z_main_lg + d_wheel_main/2 +...
    d_wheel_main/2.* sin(teta/57.3);
plot(x,z,'k','LineWidth',2)

x = x_nose_lg  + ...
    d_wheel_nose/2.* cos(teta/57.3);
z = z_nose_lg + d_wheel_nose/2 +...
    d_wheel_nose/2.* sin(teta/57.3);
plot(x,z,'k','LineWidth',2)


grid on
title('Aircraft Side-View')
xlabel({'x (m)'})
ylabel({'z (m)'})
axis equal


% Front view
subplot(2,2,3)

%fuselage
teta = linspace (0,360,100);
y = df/2 .* cos(teta/57.3);
z = df/2 .* sin(teta/57.3);
plot(y,z,'b',"LineWidth",1.5)

hold on
%wing
plot([0, Aircraft.Geometry.Wing.b.value/2], [zpos*df/2, ...
    zpos*df/2+...
    Aircraft.Geometry.Wing.b.value/2*tan(Aircraft.Geometry.Wing.dihedral.value/57.3)],'r',"LineWidth",1.5);
plot([0, -Aircraft.Geometry.Wing.b.value/2], [zpos*df/2, ...
    zpos*df/2+...
    Aircraft.Geometry.Wing.b.value/2*tan(Aircraft.Geometry.Wing.dihedral.value/57.3)],'r',"LineWidth",1.5);

%horizontal
if cfg(4:5)== 'ct'
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

%vertical
plot([0 0], [df/2, df/2+Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)

%engine
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

%undercarriage
%main
plot ([y_main_lg y_main_lg], [z_main_lg (z_main_lg+d_wheel_main)],'k',"LineWidth",2)
plot ([y_main_lg-wheel_width y_main_lg-wheel_width], [z_main_lg (z_main_lg+d_wheel_main)],'k',"LineWidth",2)

plot ([-y_main_lg -y_main_lg], [z_main_lg (z_main_lg+d_wheel_main)],'k',"LineWidth",2)
plot ([-y_main_lg+wheel_width -y_main_lg+wheel_width], [z_main_lg (z_main_lg+d_wheel_main)],'k',"LineWidth",2)
%nose
plot ([y_nose_lg-0.5*wheel_width y_nose_lg-0.5*wheel_width], [z_nose_lg (z_nose_lg+d_wheel_nose)],'k',"LineWidth",2)
plot ([y_nose_lg+0.5*wheel_width y_nose_lg+0.5*wheel_width], [z_nose_lg (z_nose_lg+d_wheel_nose)],'k',"LineWidth",2)


grid on
title('Aircraft Front-View')
xlabel({'y (m)'})
ylabel({'z (m)'})
axis equal

%saveas(gcf,[pwd 'Aircraft.fig']);
%saveas(gcf,[pwd 'Aircraft.png']);
saveas(gcf, 'Aircraft.fig');
saveas(gcf, 'Aircraft.png');

pwd 
hold off

%% 3D OUTPUT
% 3D model plot based on Aircraft structure and imported airfoils

%WING
%Aifoils
root_coord = root_coord_w;
tip_coord = tip_coord_w;

%root airfoil
root_coord_3D(:,1) = root_coord(:,1).*Aircraft.Geometry.Wing.croot.value +...
    Aircraft.Geometry.Wing.xle.value.*lf;

root_coord_3D(:,2) = root_coord(:,1).*Aircraft.Geometry.Wing.ypos.value*Aircraft.Geometry.Wing.b.value;

root_coord_3D(:,3) =  root_coord(:,2).*Aircraft.Geometry.Wing.croot.value +...
    Aircraft.Geometry.Wing.zpos.value.*df/2;


%tip airfoil
tip_coord_3D(:,1) = tip_coord(:,1).*Aircraft.Geometry.Wing.ctip.value +...
    Aircraft.Geometry.Wing.xtip_le.value.*lf;

tip_coord_3D(:,2) = Aircraft.Geometry.Wing.b.value/2;

tip_coord_3D(:,3) =  tip_coord(:,2).*Aircraft.Geometry.Wing.ctip.value +...
    Aircraft.Geometry.Wing.zpos.value.*df/2+...
    Aircraft.Geometry.Wing.b.value/2*tan(Aircraft.Geometry.Wing.dihedral.value/57.3);

subplot(2,2,4)
plot(root_coord_3D(:,1),root_coord_3D(:,3),'r',"LineWidth",2)
plot(tip_coord_3D(:,1),tip_coord_3D(:,3),'r',"LineWidth",2)

    XXw = [root_coord_3D(:,1), tip_coord_3D(:,1)]';
    YYw = [root_coord_3D(:,2), tip_coord_3D(:,2)]';
    ZZw = [root_coord_3D(:,3), tip_coord_3D(:,3)]';

subplot(2,2,2)
surf(XXw,YYw,ZZw,'FaceColor','red','EdgeColor','none')

hold on
surf(XXw,-YYw,ZZw,'FaceColor','red','EdgeColor','none')
daspect([1.0 , 1.0, 1.0])

%
%HORIZONTAL

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

subplot(2,2,4)
plot(root_coord_3D(:,1),root_coord_3D(:,3),'g',"LineWidth",2)
plot(tip_coord_3D(:,1),tip_coord_3D(:,3),'g',"LineWidth",2)

    XXh = [root_coord_3D(:,1), tip_coord_3D(:,1)]';
    YYh = [root_coord_3D(:,2), tip_coord_3D(:,2)]';
    ZZh = [root_coord_3D(:,3), tip_coord_3D(:,3)]';


subplot(2,2,2)
surf(XXh,YYh,ZZh,'FaceColor','green','EdgeColor','none')

surf(XXh,-YYh,ZZh,'FaceColor','gree','EdgeColor','none')


%
%VERTICAL

root_coord = root_coord_v;
tip_coord = tip_coord_v;

%root airfoil
root_coord_3D(:,1) = root_coord(:,1).*Aircraft.Geometry.Vertical.croot.value +...
    Aircraft.Geometry.Vertical.xle.value.*lf;

root_coord_3D(:,2) =  root_coord(:,2).*Aircraft.Geometry.Vertical.croot.value;

root_coord_3D(:,3) =  Aircraft.Geometry.Vertical.zpos.value*...
    df/2;

%tip airfoil
tip_coord_3D(:,1) = tip_coord(:,1).*Aircraft.Geometry.Vertical.ctip.value +...
    Aircraft.Geometry.Vertical.xtip_le.value.*lf;

tip_coord_3D(:,2) = tip_coord(:,2).*Aircraft.Geometry.Vertical.ctip.value;

tip_coord_3D(:,3) =  Aircraft.Geometry.Vertical.zpos.value*...
    df/2+Aircraft.Geometry.Vertical.b.value;

subplot(2,2,1)
plot(root_coord_3D(:,2),root_coord_3D(:,1),'y',"LineWidth",2)
plot(tip_coord_3D(:,2),tip_coord_3D(:,1),'y',"LineWidth",2)

    XXv = [root_coord_3D(:,1), tip_coord_3D(:,1)]';
    YYv = [root_coord_3D(:,2), tip_coord_3D(:,2)]';
    ZZv = [root_coord_3D(:,3), tip_coord_3D(:,3)]';

subplot(2,2,2)
surf(XXv,YYv,ZZv,'FaceColor','yellow','EdgeColor','none')



%FUSELAGE 
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

hold on

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


% % Distributed propulsion
% if isfield(Aircraft.Powertrain,'Secondary')
%     phi=[Aircraft.Powertrain.Secondary.Power.Alternate.phi_alternate.value,...
%         Aircraft.Powertrain.Secondary.Power.Climb.phi_climb.value,...
%         Aircraft.Powertrain.Secondary.Power.Cruise.phi_cruise.value,...
%         Aircraft.Powertrain.Secondary.Power.Descent.phi_descent.value,...
%         Aircraft.Powertrain.Secondary.Power.Landing.phi_landing.value,...
%         Aircraft.Powertrain.Secondary.Power.Loiter.phi_loiter.value,...
%         Aircraft.Powertrain.Secondary.Power.Takeoff.phi_takeoff.value];
%     if any(phi)
%         [Xe_dp,Ye_dp,Ze_dp]=DistributedProulsion3D_fun(Aircraft);
%         for j=1:Aircraft.Geometry.Engine.Secondary.N.value/2
%             surf(Xe_dp(:,:,j),Ye_dp(:,:,j),Ze_dp(:,:,j),'FaceColor','magenta','EdgeColor','none')
%             hold on
%             surf(Xe_dp(:,:,j),-Ye_dp(:,:,j),Ze_dp(:,:,j),'FaceColor','magenta','EdgeColor','none')
%         end
%     end
% end


%undercarriage
%main
nlgsection = 2; % Number of engine sections; 

    Xlg = ones (nlgsection,points);
    Ylg = ones (nlgsection,points);
    Zlg = ones (nlgsection,points);

    ymin = y_main_lg - wheel_width;
    ymax = y_main_lg;
    
    for i=1:nlgsection
        
        if i == 1
            Ylg(i,1:end) = ymin;
        end
        
        if i == nlgsection
            Ylg(i,1:end) = ymax;
        end
        secdiameter = Aircraft.Geometry.Undercarriage.Main.diameter.value;
        x = secdiameter / 2* cos(teta/57.3);
        z = secdiameter / 2* sin(teta/57.3);
        Xlg(i,1:end) = x+x_main_lg;
        Zlg(i,1:end) = z+z_main_lg + d_wheel_main/2;
        
    end

surf(Xlg,Ylg,Zlg,'FaceColor','black','EdgeColor','none')
surf(Xlg,-Ylg,Zlg,'FaceColor','black','EdgeColor','none')

%nose
nlgsection = 2; % Number of engine sections; 

    Xnlg = ones (nlgsection,points);
    Ynlg = ones (nlgsection,points);
    Znlg = ones (nlgsection,points);

    ymin = y_nose_lg - wheel_width/2;
    ymax = y_nose_lg + wheel_width/2;
    
    for i=1:nlgsection
        
        if i == 1
            Ynlg(i,1:end) = ymin;
        end
        
        if i == nlgsection
            Ynlg(i,1:end) = ymax;
        end
        secdiameter = Aircraft.Geometry.Undercarriage.Nose.diameter.value;
        x = secdiameter / 2* cos(teta/57.3);
        z = secdiameter / 2* sin(teta/57.3);
        Xnlg(i,1:end) = x+x_nose_lg;
        Znlg(i,1:end) = z+z_nose_lg + d_wheel_nose/2;
        
    end

surf(Xnlg,Ynlg,Znlg,'FaceColor','black','EdgeColor','none')
surf(Xnlg,-Ynlg,Znlg,'FaceColor','black','EdgeColor','none')

xlabel({'x (m)'})
ylabel({'y (m)'})
zlabel({'z (m)'})

disp('Saving Aircraft figures into /Geoemtry directory')
% saveas(gcf,[pwd 'Aircraft3D.fig']);
% saveas(gcf,[pwd 'Aircraft3D.png']);
saveas(gcf,'Aircraft3D.fig');
saveas(gcf,'Aircraft3D.png');
%copyfile ('_PreDesignOutput',pwd,'../_Utilities' )

%saving math file for CAD
cd ([pwd])
disp('Saving Aircraft structure for CAD into /Geoemtry directory')
disp('...')
%wing
if isfield (Aircraft.Geometry,'Wing') == 1
x = XXw';
y = YYw';
z = ZZw';
disp('...wing.mat saved')
save wing.mat x y z
end 

%horizontal
if isfield (Aircraft.Geometry,'Horizontal') == 1
x = XXh';
y = YYh';
z = ZZh';
disp('...horizontal.mat saved')
save horizontal.mat x y z
end 

%vertical
if isfield (Aircraft.Geometry,'Vertical') == 1
x = XXv';
y = YYv';
z = ZZv';
disp('...vertical.mat saved')
save vertical.mat x y z
end 

%fuselage
if isfield (Aircraft.Geometry,'Fuselage') == 1
x = Xf';
y = Yf';
z = Zf';
disp('...fuselage.mat saved')
save fuselage.mat x y z
end

%engine
if isfield (Aircraft.Geometry,'Engine') == 1
x = Xe';
y = Ye';
z = Ze';
disp('...engine.mat saved')
save engine.mat x y z
end

%engine
if isfield (Aircraft.Geometry,'Undercarriage') == 1
x = Xlg';
y = Ylg';
z = Zlg';
disp('...landing_gear.mat saved')
save landing_gear.mat x y z
end

%cd ('../')


end