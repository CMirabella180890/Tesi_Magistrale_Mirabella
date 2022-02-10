function AircraftOut(Aircraft)
%%Predesign Output function
% This function creates figures of 3-View of the Aircraft and 3D Solid PLOT
% of the Aircraft based on here imported airfoils
%
% 1) Auxiliary variables are firstly defined 
% 2) 2D plot -> 3-Views of aircraft drafting
% 3) 3D plot -> 3D PLOT of aircraft structure with imported airfoils

close all
disp ('-----------------------------------------------------------------')
disp ('-----------------------------------------------------------------')
disp ('Whole Aircraft Geometry utilities running')

%% Auxiliary variables
%Tlar and Configurations
cfg = 'lw_tt_bm';                           %aircraft configuration TO BE EDITED from outside
Aircraft.TLAR.Type_engine.value = "TP";     %aircraft configuration TO BE EDITED from outside
prop_config = 'pusher';                     % 'tractor' or 'pusher' to be provided from TLARs

%fuselage
df = Aircraft.Geometry.Fuselage.diameter.value;         % fuselage diameter (m)
dbase = 0.1*df;                                         %fuselage base diameter
lf = Aircraft.Geometry.Fuselage.length.value;           % fuselage lenght (m)
Aircraft.Geometry.Fuselage.kcockpit.value = 1.0;
Aircraft.Geometry.Fuselage.ktail.value = 1.5;

%wing
zpos = Aircraft.Geometry.Wing.zle.value;                            % wing root chord zeta position (unit)
xle = Aircraft.Geometry.Wing.xle.value/lf;                          % wing tip leading edge %of fuselage lenght
xtip_le = xle + ...
           (Aircraft.Geometry.Wing.b.value/2*...
           tan(Aircraft.Geometry.Wing.sweep.value/57.3))/lf;        % wing tip leading edge %of fuselage lenght
Aircraft.Geometry.Wing.ypos.value = 0.0;                            % %of wing semispan
%kink 1
ckink1 = Aircraft.Geometry.Wing.chord_kink_one.value;                                   % kink 1 chord in meters
ykink1 = Aircraft.Geometry.Wing.panel_span1.value * Aircraft.Geometry.Wing.b.value/2;   % kink 1 wing semispan in meters
%kink 2
ckink2 = Aircraft.Geometry.Wing.chord_kink_two.value;                                                                               % kink 2 chord in meters  
ykink2 = (Aircraft.Geometry.Wing.panel_span1.value + Aircraft.Geometry.Wing.panel_span2.value) * Aircraft.Geometry.Wing.b.value/2;  % kink 2 wing semispan in meters
%sweep_le
sweep_le = Aircraft.Geometry.Wing.sweep.value;

%horizontal
hzpos = (Aircraft.Geometry.Horizontal.zloc.value - df/2*Aircraft.Geometry.Vertical.zpos.value  ) / Aircraft.Geometry.Vertical.b.value;   %horizontal root chord zeta position %of vtail span (poitive upward)
Aircraft.Geometry.Horizontal.xle.value = Aircraft.Geometry.Horizontal.xloc.value/lf;      % % of lf
xtip_le_h = Aircraft.Geometry.Horizontal.xle.value + ...
    Aircraft.Geometry.Horizontal.b.value/2*tan(Aircraft.Geometry.Horizontal.sweep.value/57.3)/lf;                                   % % of lf
Aircraft.Geometry.Horizontal.xtip_le.value = xtip_le_h; % % of lf
%Aircraft.Geometry.Horizontal.b.value=1.496; 
Aircraft.Geometry.Horizontal.ypos.value = Aircraft.Geometry.Horizontal.yloc.value/Aircraft.Geometry.Wing.b.value/2; % % of wing span

%vertical
sweep_v = Aircraft.Geometry.Vertical.sweep.value; %sweep angle in degree
switch Aircraft.Geometry.Vertical.empennage_flag.value
    case 'Multiple fin'
        xle_v = xtip_le_h*lf;
        xtip_le_v = (xle_v + Aircraft.Geometry.Vertical.b.value*tan(sweep_v/57.3))/lf;
        yvloc = Aircraft.Geometry.Horizontal.b.value/2; % to be modified
        zvloc = df/2 + hzpos*Aircraft.Geometry.Vertical.b.value; %in meters!
    case 'Double fin'
        xle_v = xtip_le_h*lf;
        xtip_le_v = (xle_v + Aircraft.Geometry.Vertical.b.value*tan(sweep_v/57.3))/lf;
        yvloc = Aircraft.Geometry.Horizontal.b.value/2; % to be modified
        zvloc = df/2 + hzpos*Aircraft.Geometry.Vertical.b.value;
    case 'Single fin'
        xle_v = Aircraft.Geometry.Vertical.xle.value*lf;
        xtip_le_v = (xle_v + Aircraft.Geometry.Vertical.b.value*tan(sweep_v/57.3))/lf;
        yvloc = 0.0;
        zvloc = Aircraft.Geometry.Vertical.zpos.value*df/2;  % in meters
end

%engine 
ezpos = 1.0;       % % df engine zeta position
Aircraft.Geometry.Engine.Primary.xpos.value = 0.90;             % % lf
Aircraft.Geometry.Engine.Primary.lf.value = 0.5;                % engine lenght m
Aircraft.Geometry.Engine.Primary.ypos.value = 0.0;              % %of wing semispan
Aircraft.Geometry.Engine.Primary.df.value = 0.1;                % m
Aircraft.Geometry.Engine.Primary.propdiam.value = 0.6;          %prop diameter in meters
Aircraft.Geometry.Engine.Primary.zpos.value = 1.0;              %of fus df

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

% Importing airfoils sections for wing, horizontal, vertical...TO BE
% MODIFIED FROM INPUT data structure!!!
root_coord_w = importdata('_Airfoil\E216_104.txt'); % Root Airfoil coordiante read from file (x/c & z/c) IRON_Root.txt
tip_coord_w = importdata('_Airfoil\E216_104.txt'); % Root Airfoil coordiante read from file (x/c & z/c) IRON_Tip.txt

root_coord_h = importdata('_Airfoil\Simmetric_Root.txt'); % Root Airfoil coordiante read from file (x/c & z/c)
tip_coord_h = importdata('_Airfoil\Simmetric_Tip.txt'); % Root Airfoil coordiante read from file (x/c & z/c)

root_coord_v = importdata('_Airfoil\Simmetric_Root.txt'); % Root Airfoil coordiante read from file (x/c & z/c)
tip_coord_v = importdata('_Airfoil\Simmetric_Tip.txt'); % Root Airfoil coordiante read from file (x/c & z/c)


%% 2D Planform Aicraft 3-views
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Planform view %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Top-View','NumberTitle','off');
%subplot(2,2,1)

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
switch Aircraft.Geometry.Vertical.empennage_flag.value
    case 'Multiple fin'
        plot([yvloc, yvloc], xtip_le_h*lf + ...
            [0, Aircraft.Geometry.Vertical.croot.value],'y',"LineWidth",2)
        plot([-yvloc, -yvloc], xtip_le_h*lf + ...
            [0, Aircraft.Geometry.Vertical.croot.value],'y',"LineWidth",2)
    case 'Double fin'
        plot([yvloc, yvloc], xtip_le_h*lf + ...
            [0, Aircraft.Geometry.Vertical.croot.value],'y',"LineWidth",2)
        plot([-yvloc, -yvloc], xtip_le_h*lf + ...
            [0, Aircraft.Geometry.Vertical.croot.value],'y',"LineWidth",2)
    case 'Single fin'
        plot([yvloc, yvloc], Aircraft.Geometry.Vertical.xle.value*lf + ...
            [0, Aircraft.Geometry.Vertical.croot.value],'y',"LineWidth",2)
end


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

%MOVABLES START HERE

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
    
    xle_inner = Aircraft.Geometry.Wing.xle.value +  y_inner * tan(Aircraft.Geometry.Wing.sweep.value/57.3) +...
        c_inner - cf_inner;
    
    xle_outer = Aircraft.Geometry.Wing.xle.value +  y_outer * tan(Aircraft.Geometry.Wing.sweep.value/57.3) +...
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
    
    xle_inner = Aircraft.Geometry.Wing.xle.value +  y_inner * tan(Aircraft.Geometry.Wing.sweep.value/57.3) +...
        c_inner - cf_inner;
    
    xle_outer = Aircraft.Geometry.Wing.xle.value +  y_outer * tan(Aircraft.Geometry.Wing.sweep.value/57.3) +...
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
title('Aircraft Top-View')
xlabel({'y (m)'})
ylabel({'x (m)'})
axis equal

saveas(gcf, 'Top-View.fig')
saveas(gcf, 'Top-View.png')

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Side view %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%subplot(2,2,4)

figure('Name','Side-View','NumberTitle','off');

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

%wing
plot([xle*lf, xle*lf+Aircraft.Geometry.Wing.croot.value]...
    ,[zpos, zpos],'r',"LineWidth",1.5)

%horizontal
plot([Aircraft.Geometry.Horizontal.xle.value*lf,...
    Aircraft.Geometry.Horizontal.xle.value*lf+Aircraft.Geometry.Horizontal.croot.value]...
    ,[df/2+hzpos*Aircraft.Geometry.Vertical.b.value, ...
    df/2+hzpos*Aircraft.Geometry.Vertical.b.value],'g',"LineWidth",1.5)

%vertical
switch Aircraft.Geometry.Vertical.empennage_flag.value
    case 'Multiple fin'
        plot([xtip_le_h*lf, xtip_le_h*lf+Aircraft.Geometry.Vertical.croot.value],...
            [zvloc, zvloc],'y',"LineWidth",1.5)
        plot([xtip_le_v*lf, xtip_le_v*lf+Aircraft.Geometry.Vertical.ctip.value],...
            [zvloc + Aircraft.Geometry.Vertical.b.value, zvloc + Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)
        plot([xtip_le_h*lf, xtip_le_v*lf],...
            [zvloc, zvloc + Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)        
        plot([xtip_le_h*lf+Aircraft.Geometry.Vertical.croot.value, xtip_le_v*lf+Aircraft.Geometry.Vertical.ctip.value],...
            [zvloc, zvloc + Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)
    case 'Double fin'
        plot([xtip_le_h*lf, xtip_le_h*lf+Aircraft.Geometry.Vertical.croot.value],...
            [zvloc, zvloc],'y',"LineWidth",1.5)
         plot([xtip_le_v*lf, xtip_le_v*lf+Aircraft.Geometry.Vertical.ctip.value],...
            [zvloc + Aircraft.Geometry.Vertical.b.value, zvloc + Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)
          plot([xtip_le_h*lf, xtip_le_v*lf],...
            [zvloc, zvloc + Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)        
        plot([xtip_le_h*lf+Aircraft.Geometry.Vertical.croot.value, xtip_le_v*lf+Aircraft.Geometry.Vertical.ctip.value],...
            [zvloc, zvloc + Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)
    case 'Single fin'
        plot([xtip_le_h*lf, xtip_le_h*lf+Aircraft.Geometry.Vertical.croot.value],...
            [zvloc, zvloc],'y',"LineWidth",1.5)
        plot([xtip_le_v*lf, xtip_le_v*lf+Aircraft.Geometry.Vertical.ctip.value],...
            [zvloc + Aircraft.Geometry.Vertical.b.value, zvloc + Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)
           plot([xtip_le_h*lf, xtip_le_v*lf],...
            [zvloc, zvloc + Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)        
        plot([xtip_le_h*lf+Aircraft.Geometry.Vertical.croot.value, xtip_le_v*lf+Aircraft.Geometry.Vertical.ctip.value],...
            [zvloc, zvloc + Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)
end


%RUDDER
%linear chord law
c0 = Aircraft.Geometry.Vertical.croot.value;          % root chord (unit)
c1 = Aircraft.Geometry.Vertical.ctip.value;           % tip chord (unit)
span0 = 0.0;
span1 = Aircraft.Geometry.Vertical.b.value;
[c_0,c_y] = chord_linear(c0 ,c1, span0, span1);
%
y_inner = Aircraft.Geometry.Rudder.eta_inner.value * span1;
y_outer = Aircraft.Geometry.Rudder.eta_outer.value * span1;
%
c_inner = c_0 - c_y * y_inner;                      %vtail chord at inner station (unit)
c_outer = c_0 - c_y * y_outer;                      %vtail chord at inner station (unit)
%
cf_inner = c_inner * Aircraft.Geometry.Rudder.cr_c_root.value;   %rudder chord at inner station (unit)
cf_outer = c_outer* Aircraft.Geometry.Rudder.cr_c_tip.value;     %rudder chord at inner station (unit)
%
xle_inner = xtip_le_h*lf +  y_inner * tan(sweep_v/57.3) + (c_inner - cf_inner);
xle_outer = xtip_le_h*lf +  y_outer * tan(sweep_v/57.3) + (c_outer - cf_outer);
%

plot([xle_inner  xle_inner+cf_inner],...
     zvloc + [y_inner y_inner],'y',"LineWidth",1.5)
plot([xle_outer xle_outer+cf_outer],...
     zvloc + [y_outer y_outer],'y',"LineWidth",1.5)
 plot([xle_inner xle_outer],...
     zvloc + [y_inner y_outer],'y',"LineWidth",1.5)

fill ([xle_inner xle_outer xle_outer+cf_outer xle_inner+cf_inner],...
      [zvloc + y_inner  zvloc + y_outer...
       zvloc + y_outer  zvloc + y_inner],'y')

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

%LANDING GEAR
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

saveas(gcf, 'Side-View.fig')
saveas(gcf, 'Side-View.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Front view
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subplot(2,2,3)
figure('Name','Front-View','NumberTitle','off');

%fuselage
teta = linspace (0,360,100);
y = df/2 .* cos(teta/57.3);
z = df/2 .* sin(teta/57.3);
plot(y,z,'b',"LineWidth",1.5)

hold on
%wing
plot([0, Aircraft.Geometry.Wing.b.value/2], [zpos, ...
    zpos+...
    Aircraft.Geometry.Wing.b.value/2*tan(Aircraft.Geometry.Wing.dihedral.value/57.3)],'r',"LineWidth",1.5);
plot([0, -Aircraft.Geometry.Wing.b.value/2], [zpos, ...
    zpos+...
    Aircraft.Geometry.Wing.b.value/2*tan(Aircraft.Geometry.Wing.dihedral.value/57.3)],'r',"LineWidth",1.5);

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

%vertical
switch Aircraft.Geometry.Vertical.empennage_flag.value
    case 'Multiple fin'
        plot([yvloc yvloc], [zvloc, zvloc + Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)
        plot([-yvloc -yvloc], [zvloc, zvloc + Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)
    case 'Double fin'
        plot([yvloc yvloc], [zvloc, zvloc + Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)
        plot([-yvloc -yvloc], [zvloc, zvloc + Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)
    case 'Single fin'
        plot([yvloc yvloc], [zvloc, zvloc + Aircraft.Geometry.Vertical.b.value],'y',"LineWidth",1.5)
end

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
saveas(gcf, 'Front-View.fig')
saveas(gcf, 'Front-View.png')

pwd 
hold off

%% 3D OUTPUT
% 3D model plot based on Aircraft structure and imported airfoils

figure('Name','3D-View','NumberTitle','off');

%WING
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

%subplot(2,2,4)
% plot(root_coord_3D(:,1),root_coord_3D(:,3),'r',"LineWidth",2)
% plot(tip_coord_3D(:,1),tip_coord_3D(:,3),'r',"LineWidth",2)

    XXw = [root_coord_3D(:,1), tip_coord_3D(:,1)]';
    YYw = [root_coord_3D(:,2), tip_coord_3D(:,2)]';
    ZZw = [root_coord_3D(:,3), tip_coord_3D(:,3)]';

%subplot(2,2,2)
surf(XXw,YYw,ZZw,'FaceColor','red','EdgeColor','none')

hold on
surf(XXw,-YYw,ZZw,'FaceColor','red','EdgeColor','none')
daspect([1.0 , 1.0, 1.0])

%
%HORIZONTAL
clear root_coord_3D;
clear tip_coord_3D;

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

% %subplot(2,2,4)
% plot(root_coord_3D(:,1),root_coord_3D(:,3),'g',"LineWidth",2)
% plot(tip_coord_3D(:,1),tip_coord_3D(:,3),'g',"LineWidth",2)

    XXh = [root_coord_3D(:,1), tip_coord_3D(:,1)]';
    YYh = [root_coord_3D(:,2), tip_coord_3D(:,2)]';
    ZZh = [root_coord_3D(:,3), tip_coord_3D(:,3)]';


%subplot(2,2,2)
surf(XXh,YYh,ZZh,'FaceColor','green','EdgeColor','none')
surf(XXh,-YYh,ZZh,'FaceColor','gree','EdgeColor','none')


%
%VERTICAL
clear root_coord_3D;
clear tip_coord_3D;

root_coord = root_coord_v;
tip_coord = tip_coord_v;


switch Aircraft.Geometry.Vertical.empennage_flag.value
    case 'Multiple fin'
        %x-root
        root_coord_3D(:,1) = root_coord(:,1).*Aircraft.Geometry.Vertical.croot.value +...
            xle_v;
        %y-root
        root_coord_3D(:,2) = yvloc + root_coord(:,2).*Aircraft.Geometry.Vertical.croot.value;
        %z-root
        root_coord_3D(:,3) =  zvloc;
        %x-tip
        tip_coord_3D(:,1) = tip_coord(:,1).*Aircraft.Geometry.Vertical.ctip.value +...
            xtip_le_v;
        %y-tip
        tip_coord_3D(:,2) = yvloc + tip_coord(:,2).*Aircraft.Geometry.Vertical.ctip.value;
        %z-tip
        tip_coord_3D(:,3) =  zvloc + Aircraft.Geometry.Vertical.b.value;
        
        %subplot(2,2,1)
%         plot(root_coord_3D(:,2),root_coord_3D(:,1),'y',"LineWidth",2)
%         plot(tip_coord_3D(:,2),tip_coord_3D(:,1),'y',"LineWidth",2)
        
        XXv = [root_coord_3D(:,1), tip_coord_3D(:,1)]';
        YYv = [root_coord_3D(:,2), tip_coord_3D(:,2)]';
        ZZv = [root_coord_3D(:,3), tip_coord_3D(:,3)]';
        
        %subplot(2,2,2)
        surf(XXv,YYv,ZZv,'FaceColor','yellow','EdgeColor','none')
        surf(XXv,-YYv,ZZv,'FaceColor','yellow','EdgeColor','none')
        
    case 'Double fin'
        %x-root
        root_coord_3D(:,1) = root_coord(:,1).*Aircraft.Geometry.Vertical.croot.value +...
            xtip_le_h.*lf; %same of htail tip (assumptions)
        %y-root
        root_coord_3D(:,2) = yvloc + root_coord(:,2).*Aircraft.Geometry.Vertical.croot.value;
        %z-root
        root_coord_3D(:,3) =  zvloc;
        %x-tip
        tip_coord_3D(:,1) = tip_coord(:,1).*Aircraft.Geometry.Vertical.ctip.value +...
            xtip_le_v.*lf;
        %y-tip
        tip_coord_3D(:,2) = yvloc + tip_coord(:,2).*Aircraft.Geometry.Vertical.ctip.value;
        %z-tip
        tip_coord_3D(:,3) =  zvloc + Aircraft.Geometry.Vertical.b.value;
        
        %subplot(2,2,1)
%         plot(root_coord_3D(:,2),root_coord_3D(:,1),'y',"LineWidth",2)
%         plot(tip_coord_3D(:,2),tip_coord_3D(:,1),'y',"LineWidth",2)
        
        XXv = [root_coord_3D(:,1), tip_coord_3D(:,1)]';
        YYv = [root_coord_3D(:,2), tip_coord_3D(:,2)]';
        ZZv = [root_coord_3D(:,3), tip_coord_3D(:,3)]';
        
        %subplot(2,2,2)
        surf(XXv,YYv,ZZv,'FaceColor','yellow','EdgeColor','none')
        surf(XXv,-YYv,ZZv,'FaceColor','yellow','EdgeColor','none')
        
    case 'Single fin'
        %x-root
        root_coord_3D(:,1) = root_coord(:,1).*Aircraft.Geometry.Vertical.croot.value +...
            xle_v;
        %y-root
        root_coord_3D(:,2) = yvloc + root_coord(:,2).*Aircraft.Geometry.Vertical.croot.value;
        %z-root
        root_coord_3D(:,3) =  zvloc;
        %x-tip
        tip_coord_3D(:,1) = tip_coord(:,1).*Aircraft.Geometry.Vertical.ctip.value +...
            xtip_le_v.*lf;
        %y-tip
        tip_coord_3D(:,2) = yvloc + tip_coord(:,2).*Aircraft.Geometry.Vertical.ctip.value;
        %z-tip
        tip_coord_3D(:,3) =  zvloc + Aircraft.Geometry.Vertical.b.value;
        
        %subplot(2,2,1)
%         plot(root_coord_3D(:,2),root_coord_3D(:,1),'y',"LineWidth",2)
%         plot(tip_coord_3D(:,2),tip_coord_3D(:,1),'y',"LineWidth",2)
        
        XXv = [root_coord_3D(:,1), tip_coord_3D(:,1)]';
        YYv = [root_coord_3D(:,2), tip_coord_3D(:,2)]';
        ZZv = [root_coord_3D(:,3), tip_coord_3D(:,3)]';
        
        %subplot(2,2,2)
        surf(XXv,YYv,ZZv,'FaceColor','yellow','EdgeColor','none')
        surf(XXv,-YYv,ZZv,'FaceColor','yellow','EdgeColor','none')        
end




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

saveas(gcf,'Aircraft3D.fig');
saveas(gcf,'Aircraft3D.png');
%copyfile ('_PreDesignOutput',pwd,'../_Utilities' )
hold off

%% 3D SOLID 4 - VIEWS

%SIDE VIEW
open ('Aircraft3D.fig')
%figure('Name','Aircraft Top-View','NumberTitle','off');
title('Aircraft 3D Side-View')
view(0,0.5)
saveas(gcf,'SideView3D.fig');
saveas(gcf,'SideView3D.png');


%FRONT-VIEW
open ('Aircraft3D.fig')
%figure('Name','Aircraft Front-View','NumberTitle','off');
title('Aircraft 3D Front-View')
view(-90,0.0)
saveas(gcf,'FrontView3D.fig');
saveas(gcf,'FrontView3D.png');

%FRONT-VIEW
open ('Aircraft3D.fig')
%figure('Name','Aircraft Top-View','NumberTitle','off');
title('Aircraft 3D Top-View')
view(0,90.0)
saveas(gcf,'TopView3D.fig');
saveas(gcf,'TopView3D.png');



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

  %% moving results
    movefile('Top-View.png', Aircraft.res_dir);
    movefile('Top-View.fig', Aircraft.res_dir);
    movefile('Side-View.fig', Aircraft.res_dir);
    movefile('Side-View.png', Aircraft.res_dir);
    movefile('Front-View.fig', Aircraft.res_dir);
    movefile('Front-View.png', Aircraft.res_dir);
    movefile('Aircraft3D.fig', Aircraft.res_dir);
    movefile('Aircraft3D.png', Aircraft.res_dir);
    movefile('FrontView3D.fig', Aircraft.res_dir);
    movefile('FrontView3D.png', Aircraft.res_dir);
    movefile('SideView3D.fig', Aircraft.res_dir);
    movefile('SideView3D.png', Aircraft.res_dir);
    movefile('TopView3D.fig', Aircraft.res_dir);
    movefile('TopView3D.png', Aircraft.res_dir);
    
    if exist('engine.mat') == 2
        movefile('engine.mat', Aircraft.res_dir);
    end
    if exist('fuselage.mat') == 2
        movefile('fuselage.mat', Aircraft.res_dir);
    end
    if exist('wing.mat') == 2
        movefile('wing.mat', Aircraft.res_dir);
    end
    if exist('vertical.mat') == 2
        movefile('vertical.mat', Aircraft.res_dir);
    end
    if exist('horizontal.mat') == 2
        movefile('horizontal.mat', Aircraft.res_dir);
    end
    if exist('landing_gear.mat') == 2
        movefile('landing_gear.mat', Aircraft.res_dir);
    end
    
disp ('-----------------------------------------------------------------')
disp ('-----------------------------------------------------------------')
disp ('End Geometry utilities run')

end

% % % %WING
% % % surf(XXw,YYw,ZZw,'FaceColor','red','EdgeColor','none')
% % % daspect([1.0 1.0 1.0])
% % % hold on
% % % surf(XXw,-YYw,ZZw,'FaceColor','red','EdgeColor','none')
% % % camlight left; lighting phong
% % % 
% % % %HORIZONTAL
% % % surf(XXh,YYh,ZZh,'FaceColor','green','EdgeColor','none')
% % % surf(XXh,-YYh,ZZh,'FaceColor','gree','EdgeColor','none')
% % % 
% % % %VERTICAL
% % % surf(XXv,YYv,ZZv,'FaceColor','yellow','EdgeColor','none')
% % % surf(XXv,-YYv,ZZv,'FaceColor','yellow','EdgeColor','none')
% % % 
% % % %FUSELAGE
% % % surf(Xf,Yf,Zf,'FaceColor','blu','EdgeColor','none')
% % % camlight left; lighting phong
% % % 
% % % %ENGINE
% % % surf(Xe,Ye,Ze,'FaceColor','black','EdgeColor','none')
% % % surf(Xe,-Ye,Ze,'FaceColor','black','EdgeColor','none')
% % % 
% % % %LANDING GEAR
% % % surf(Xlg,Ylg,Zlg,'FaceColor','black','EdgeColor','none')
% % % surf(Xlg,-Ylg,Zlg,'FaceColor','black','EdgeColor','none')
% % % surf(Xnlg,Ynlg,Znlg,'FaceColor','black','EdgeColor','none')
% % % surf(Xnlg,-Ynlg,Znlg,'FaceColor','black','EdgeColor','none')
% % % camlight left; lighting phong
% % % 
% % % title('Aircraft 3D Side-View')
% % % xlabel({'x (m)'})
% % % ylabel({'y (m)'})
% % % zlabel({'z (m)'})
% % % axis equal
% % % view(0,0.5)
% % % 
% % % saveas(gcf,'SideView3D.fig');
% % % saveas(gcf,'SideView3D.png');
% % % hold off