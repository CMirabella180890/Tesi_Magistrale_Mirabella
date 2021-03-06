%% CHAPTER 15 - Power plant
%Ref EASA/Tecnam p2006 report

import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

% REGULATIONS
requirement = Aircraft.Certification.Regulation.value;
csvla_333d  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.Attributes.cs;
csvla_361a1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Takeoff.Mean_torque.Attributes.cs;
csvla_361a2 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Max_Continous.Mean_torque.Attributes.cs;
csvla_361b  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Takeoff.Limit_torque.Attributes.cs;
csvla_363   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Side_loads.Total_side_load.Attributes.cs;

% chapter_number = 12;
% ch = strcat('ch' , num2str(chapter_number)); 
ch = Chapter();
ch.Title = 'Power plant';
disp(['Chapter 15', (' "'), ch.Title,('" ') ,'writing...' ])
% -------------------------------------------------------------------------
% moving to another path for figure
cd ..
cd ..
%  regulation = Aircraft.Certification.Regulation.value;
 results_path = [pwd '\utilities\Geometry\DroneVLA_results\'];
cd(RepDir);

fig = FormalImage([results_path,'Engine-Front-View.png']);
         fig.Caption = 'Engine, front view.';
         fig.Height = '5in';
         fig.LinkTarget='engine_front_view';
         add(ch,fig);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% moving to another path for figure
cd ..
cd ..
%  regulation = Aircraft.Certification.Regulation.value;
 results_path = [pwd '\utilities\Geometry\DroneVLA_results\'];
cd(RepDir);

fig = FormalImage([results_path,'Engine-Side-View.png']);
         fig.Caption = 'Engine, side view.';
         fig.Height = '5in';
         fig.LinkTarget='engine_side_view';
         add(ch,fig);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% moving to another path for figure
cd ..
cd ..
%  regulation = Aircraft.Certification.Regulation.value;
 results_path = [pwd '\utilities\Geometry\DroneVLA_results\'];
cd(RepDir);

fig = FormalImage([results_path,'Engine-Top-View.png']);
         fig.Caption = 'Engine, top view.';
         fig.Height = '5in';
         fig.LinkTarget='engine_top_view';
         add(ch,fig);
% -------------------------------------------------------------------------
str = ['The engine mount and its supporting structure must be designed' ...
       ' for the effects of:'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(ch,para);
% -------------------------------------------------------------------------
        % -----------------------------------------------------------------
        %ordered list
        ref1 = ['a limit engine torque corresponding to takeoff' ...
            ' power and propeller speed acting simultaneously' ...
            ' with 75% of the limit loads from flight condition A' ...
            ' of' ...
            (' ') ...
            char(requirement) ...
            (' ') ...
            char(csvla_333d) ...
            ', according to ' ...
            (' ') ...
            char(requirement) ...
            (' ') ...
            char(csvla_361a1) ...
            ';'];
        ref2 = ['the limit engine torque as specified in' ...
            (' ') ...
            char(requirement) ...
            (' ') ...
            char(csvla_361b) ...
            ' acting simultaneously with the limit loads' ...
            ' from flight condition A of' ...
            (' ') ...
            char(requirement) ...
            (' ') ...
            char(csvla_333d) ...
            ', according to' ...
            (' ') ...
            char(requirement) ...
            (' ') ...
            char(csvla_361a2) ...
            '.'];
        ol = OrderedList({ref1, ref2});
        %list
        append(ch,ol);
        % -----------------------------------------------------------------
str = ['The limit engine torque to be considered under' ...
    (' ') ...
    char(requirement) ...
    (' ') ...
    char(csvla_361a2) ...
    ' must be obtained by multiplying the mean torque for maximum' ...
    ' continuous power by a factor determined as follows:'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(ch,para);
        % -----------------------------------------------------------------
        %ordered list
        ref1 = ['for a four-stroke engines, (i) 1.33 for engines' ...
            ' with five or more cylinders, (ii) 2, 3, 4 or 8 for' ...
            ' engines with four, three, two or one cylinders, respectively;'];
        ref2 = ['for a two-stroke engines, (i) 2 for engines' ...
            ' with three or more cylinders, (ii) 3 or 6 for' ...
            ' engines with two or one cylinders, respectively.'];
        ol = OrderedList({ref1, ref2});
        %list
        append(ch,ol);
        % -----------------------------------------------------------------
str = ['The engine mount and its supporting structure must be' ...
    ' designed for a limit load factor in a lateral direction,' ...
    ' for the side load on the engine mount not less than 1.33 and' ...
    ' this side load may be assumed to be independent of other' ...
    ' flight conditions, according to' ...
    (' ') ...
    char(requirement) ...
    (' ') ...
    char(csvla_363) ...
    '.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(ch,para);
        % -------------------------------------------------------------------------

%sec
sec = Section();
sec.Title = 'Engine torque';

takeoff_power                      = Aircraft.Engine.Takeoff.Power.value;
max_continuous_power               = Aircraft.Engine.Max_Continous.Power.value;
takeoff_engine_rpm                 = Aircraft.Engine.Takeoff.RPM.value;
reduction_ratio                    = Aircraft.Engine.Reduction_ratio.value;
propeller_rotational_speed         = Aircraft.Engine.Takeoff.Propeller_rotational_speed.value;
mean_engine_torque_max_power       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Takeoff.Mean_torque.value;
mean_engine_torque_max_continuous  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Max_Continous.Mean_torque.value;
limit_engine_torque_max_power      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Takeoff.Limit_torque.value;
limit_engine_torque_max_continuous = Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Max_Continous.Limit_torque.value;

str = ['The engine takeoff power is ' ...
    strcat(num2str(takeoff_power,4)) ... 
    (' ') ...
    char(Aircraft.Engine.Takeoff.Power.Attributes.unit) ...
    ' at ' ...
    strcat(num2str(takeoff_engine_rpm,4)) ...
    (' ') ...
    char(Aircraft.Engine.Takeoff.RPM.Attributes.unit) ...
    '. The rotational speed of the propeller is ' ...
    strcat(num2str(takeoff_engine_rpm,4)) ...
    '/' ...
    strcat(num2str(reduction_ratio,4)) ...
    ' = ' ...
    strcat(num2str(propeller_rotational_speed,4)) ...
    (' ') ...
    char(Aircraft.Engine.Takeoff.Propeller_rotational_speed.Attributes.unit) ...
    '. The maximum continuous power is ' ...
    (' ') ...
    strcat(num2str(max_continuous_power,4)) ...
    (' ') ...
    char(Aircraft.Engine.Max_Continous.Power.Attributes.unit) ...
    '. The mean engine torque is ' ...
    strcat(num2str(mean_engine_torque_max_power,4)) ...
    (' ') ...
    char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Takeoff.Mean_torque.Attributes.unit) ...
    '. Using a factor of ' ...
    strcat(num2str(Aircraft.Engine.Correction_factor.value,4)) ... 
    ' for a four cylinder engine, the limit torque will be ' ...
    strcat(num2str(limit_engine_torque_max_power,4)) ...
    (' ') ...
    char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Takeoff.Limit_torque.Attributes.unit) ...
    '. This limit torque acts simultaneously with the 75 % of the inertia limit load. ' ...
    ' The mean engine torque at max continuous power is ' ...
    (' ') ...
    strcat(num2str(mean_engine_torque_max_continuous,4))  ...
    (' ') ...
    char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Max_Continous.Mean_torque.Attributes.unit) ...
    '. Using a factor of ' ...
    (' ') ...
    strcat(num2str(Aircraft.Engine.Correction_factor.value,4)) ...
    ' for a four cylinder engine, the limit torque will be ' ...
    (' ') ...
    strcat(num2str(limit_engine_torque_max_continuous,4)) ...
    (' ') ...
    char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Max_Continous.Limit_torque.Attributes.unit) ...
    (' ') ...
    ' which acts simultaneously with the 100 % of the inertia limit load.' ];
          para = Paragraph(str);
          para.Style = {HAlign('justify')};
          add(sec,para);
          add(ch,sec);
          
        % MAX CONTINUOUS  
        myNumEq = strcat( ' = ', ...
                                num2str(max_continuous_power,4), ' * \frac{1000}', ...
                                '{ \frac{2 * 3.14 * ', num2str(propeller_rotational_speed,4), ...
                                '}{60}} =', num2str(mean_engine_torque_max_continuous,4), ...
                                '\,\,', Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Max_Continous.Mean_torque.Attributes.unit);
        % latex interprete with $ simbol
        myEq = "$ MT_{continuous} = P_{continuous} * \frac{1000}{\frac{2\pi*RPM_{prop}}{60}}";
        eq = Equation(strcat(myEq, myNumEq));
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg);
        % -----------------------------------------------------------------
        myNumEq = strcat( ' = ', num2str(reduction_ratio,4), ' * ', ...
                          num2str(mean_engine_torque_max_continuous,4), ...
                          ' = ', num2str(limit_engine_torque_max_continuous,4), ...
                          Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Max_Continous.Limit_torque.Attributes.unit);
        % latex interprete with $ simbol
        myEq = "$ LT_{continuous} = RR_{prop} * MT_{continuous}";
        eq = Equation(strcat(myEq, myNumEq));
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg);
        
        % TAKEOFF
        myNumEq = strcat( ' = ', ...
                                num2str(takeoff_power,4), ' * \frac{1000}', ...
                                '{ \frac{2 * 3.14 * ', num2str(propeller_rotational_speed,4), ...
                                '}{60}} =', num2str(mean_engine_torque_max_power,4), ...
                                '\,\,', Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Takeoff.Mean_torque.Attributes.unit);
        % latex interprete with $ simbol
        myEq = "$ MT_{takeoff} = P_{takeoff} * \frac{1000}{\frac{2\pi*RPM_{prop}}{60}}";
        eq = Equation(strcat(myEq, myNumEq));
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg);
        % -----------------------------------------------------------------
        myNumEq = strcat( ' = ', num2str(reduction_ratio,4), ' * ', ...
                          num2str(mean_engine_torque_max_power,4), ...
                          ' = ', num2str(limit_engine_torque_max_power,4), ...
                          Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Takeoff.Limit_torque.Attributes.unit);
        % latex interprete with $ simbol
        myEq = "$ LT_{takeoff} = RR_{prop} * MT_{takeoff}";
        eq = Equation(strcat(myEq, myNumEq));
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg);
        
% str = ['ADD HERE details '];

% para = Paragraph(str);
% add(sec,para);
% add(ch,sec);

        % MEAN TORQUE - MAX POWER
        myEq = "$MT_{takeoff} = \mathrm{mean torque at takeoff power (N * m); and}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg1 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref1 = eqImg1;         

        % LIMIT TORQUE - MAX POWER
        myEq = "$LT_{takeoff} = \mathrm{limit torque at takeoff power (N * m); and}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg2 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref2 = eqImg2; 
        
        % MEAN TORQUE - MAX CONTINUOUS
        myEq = "$MT_{continuous} = \mathrm{mean torque at max continuous power (N * m); and}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg3 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref3 = eqImg3;         

        % LIMIT TORQUE - MAX CONTINUOUS
        myEq = "$LT_{continuous} = \mathrm{limit torque at max continuous power (N * m).}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg4 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref4 = eqImg4;      
        
        ol = UnorderedList({ref1, ref2, ref3, ref4});
%         ol = UnorderedList({ref1,ref2,ref3,...
%             ref4,ref5,ref6, ref7,ref8});
%         ol = UnorderedList({ref1, ref2, ref3,...
%             ref4,ref5,ref6, ref7, ref8, ref9});
        append(sec,ol);


Engine_mount_mass       = Aircraft.Engine.Engine_mount_mass.value;
Engine_accessories_mass = Aircraft.Engine.Engine_accessories_mass.value;
Propeller_spinner_mass  = Aircraft.Propeller.Propeller_spinner_mass.value;
Engine_block_mass       = Engine_mount_mass + Engine_accessories_mass + Propeller_spinner_mass;
g                       = Aircraft.Constants.g.value;
Total_side_load         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Side_loads.Total_side_load.value;

%sec
sec = Section();
sec.Title = 'Side load on engine mount';
% str = strcat(' The limit load factor in a lateral direction is 1.33. The mass ', ...
%        ' of the engine group is ', num2str(Engine_block_mass,4), ...
%        Aircraft.Engine.Engine_block_mass.Attributes.unit, '.', ...
%        ' The side load results is 1.33 * ', num2str(Engine_block_mass,4), ...
%        ' * ', num2str(g,4), ' * (1/10) = ', num2str(Total_side_load,4), ...
%        Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Side_loads.Total_side_load.Attributes.unit);
str = ['The limit load factor in a lateral direction is 1.33. The mass ' ...
    'of the engine group is ' ...
    strcat(num2str(Engine_block_mass,4)) ...
    (' ') ...
    char(Aircraft.Engine.Engine_block_mass.Attributes.unit) ...
    '. The side load results is 1.33*' ...
    strcat(num2str(Engine_block_mass,4)) ...
    '*' ...
    strcat(num2str(g,4)) ...
    '*(1/10) = ' ...
    strcat(num2str(Total_side_load,4)) ...
    (' ') ...
    char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Side_loads.Total_side_load.Attributes.unit)];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);
%add(ch,para);

add(ch,sec);


Gust_limit_load              = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.nA1.value;
Inertia_load_on_engine_mount = Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Inertia_load_on_engine_mount.value;

%sec
sec = Section();
sec.Title = 'Inertia load on engine mount';

% str = strcat('The inertia load is equal to the maximum limit load factor ', ...
%     'times the engine group weight: ', num2str(Gust_limit_load,4), ...
%     '*', num2str(Engine_block_mass,4), '*', num2str(g,4), '*(1/10) = ', ...
%     num2str(Inertia_load_on_engine_mount,4), ...
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Inertia_load_on_engine_mount.Attributes.unit);
str = ['The inertia load is equal to the maximum limit load factor ' ...
    'times the engine group weight: ' ...
    (' ') ...
    strcat(num2str(Gust_limit_load,4)) ...
    '*' ...
    strcat(num2str(Engine_block_mass,4)) ...
    '*' ...
    strcat(num2str(g,4)) ...
    '*(1/10) = ' ...
    (' ') ...
    strcat(num2str(Inertia_load_on_engine_mount,4)) ...
    (' ') ...
    char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Inertia_load_on_engine_mount.Attributes.unit)];

para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);

add(ch,sec);

%sec
sec = Section();
sec.Title = 'Gyroscopic loads';

gyro_airworth_reg = char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Gyroscopic_loads.Attributes.cs); 

str = ['According to ' ...
    (' ') ...
    gyro_airworth_reg ...
    ', for a two blade propeller, the maximum gyroscopic couple is given by: '];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);

        % GYSCOPIC COUPLE EQUATION AMC 23.371(a)
        % latex interprete with $ simbol
        myEq = "$ 2 * I_{p} * \omega_1 * \omega_2 ";
        eq = Equation(strcat(myEq));
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg);
        
str = ['Where: '];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);

        % MEAN TORQUE - MAX POWER
        myEq = "$ I_p = \mathrm{polar moment of inertia of the propeller (kg * m^2); and}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg1 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref1 = eqImg1;         

        % LIMIT TORQUE - MAX POWER
        myEq = "$ \omega_{1} = \mathrm{propeller rotation speed (rad / sec); and}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg2 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref2 = eqImg2; 
        
        % MEAN TORQUE - MAX CONTINUOUS
        myEq = "$ \omega_{2} = \mathrm{rate of pitch or yaw (rad / sec).}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg3 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref3 = eqImg3; 
        ol = UnorderedList({ref1, ref2, ref3});
%         ol = UnorderedList({ref1,ref2,ref3,...
%             ref4,ref5,ref6, ref7,ref8});
%         ol = UnorderedList({ref1, ref2, ref3,...
%             ref4,ref5,ref6, ref7, ref8, ref9});
        append(sec,ol);
     
Engine_normal_load_factor = Aircraft.Engine.Engine_normal_load_factor.value;
Engine_gyro_inertia_load  = Engine_normal_load_factor * Engine_block_mass;
        
str = ['The asymmetric flow through the propeller disc is discounted ' ...
    'because the propeller diameter is less than 2.74 m as established by' ...
    (' ') ...
    gyro_airworth_reg ...
    '. The polar moment of inertia of the propeller is' ...
    (' ') ...
    strcat(num2str(Aircraft.Propeller.Propeller_polar_moment.value,4)) ...
    (' ') ...
    char(Aircraft.Propeller.Propeller_polar_moment.Attributes.unit) ...
    '. The rate of pitch or yaw is established as 1.0 rad/sec and 2.5 rad/sec' ...
    ' respectively , the load factor is ' ...
    (' ') ...
    strcat(num2str(Engine_normal_load_factor,4)) ...
    ' and the power condition is max continuous power as prescribed in' ...
    (' ') ...
    gyro_airworth_reg ...
    '. Therefore, the inertial load is equal to' ...
    (' ') ...
    strcat(num2str(Engine_normal_load_factor,4)) ...
    '*' ...
    strcat(num2str(Engine_block_mass,4)) ...
    ' = ' ...
    (' ') ...
    strcat(num2str(Engine_gyro_inertia_load,4)) ...
    (' ') ...
    char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Inertia_load_on_engine_mount.Attributes.unit) ...
    '. The rotation of the propeller at maximum continuous power is ' ...
    (' ') ...
    strcat(num2str(Aircraft.Engine.Max_Continous.RPM.value,4)) ...
    '/' ...
    strcat(num2str(reduction_ratio,4)) ...
    ' = ' ...
    strcat(num2str(Aircraft.Engine.Max_Continous.Propeller_rotational_speed.value,4)) ...
    '. The gyroscopic couple is: '];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);    

        %ordered list
        ref1 = ['Yaw case: ' ...
            (' ') ...
            '2*' ...
            strcat(num2str(Aircraft.Propeller.Propeller_polar_moment.value,4)) ...
            '*2.5*(2*3.14/60)*' ...
            strcat(num2str(Aircraft.Engine.Max_Continous.Propeller_rotational_speed.value,4)) ...
            ' = ' ...
            strcat(num2str(Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Gyroscopic_loads.yaw_case.value,4)) ...
            (' ') ...
            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Gyroscopic_loads.yaw_case.Attributes.unit)];
        ref2 = ['Pitch case: ' ...
            (' ') ...
            '2*' ...
            strcat(num2str(Aircraft.Propeller.Propeller_polar_moment.value,4)) ...
            '*1.0*(2*3.14/60)*' ...
            strcat(num2str(Aircraft.Engine.Max_Continous.Propeller_rotational_speed.value,4)) ...
            ' = ' ...
            strcat(num2str(Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Gyroscopic_loads.pitch_case.value,4)) ...
            (' ') ...
            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Gyroscopic_loads.pitch_case.Attributes.unit)];
        ol = OrderedList({ref1, ref2});
        
%         %1
%         para = Paragraph(str);
%         para.Style = {HAlign('justify')};
%         add(sec,para)
%         %2
%         para2 = Paragraph(str2);
%         para2.Style = {HAlign('justify')};
%         add(sec,para2)
        %list
        append(sec,ol);
        
add(ch,sec);

%moving to another path for figure
% cd ..
% cd ..
%  regulation = Aircraft.Certification.Regulation.value;
%  results_path = [pwd '\' regulation '\Output\'];
% 
%  cd (RepDir);
% 
% fig = FormalImage([results_path,'Finalenvelope.pdf']);
%          fig.Caption = 'Maneuver and Gust load factors and diagram';
%          fig.Height = '5in';
%          fig.LinkTarget='maneuver_ref';
%          add(ch,fig);

%add(ch7,sec1);


%% END chapter
%Adding chapters
add(rpt,ch);