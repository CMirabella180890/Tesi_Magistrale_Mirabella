%% CHAPTER 15 - Power plant
%Ref EASA/Tecnam p2006 report

import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

% chapter_number = 12;
% ch = strcat('ch' , num2str(chapter_number)); 
ch = Chapter();
ch.Title = 'Power plant';
disp(['Chapter 15', (' "'), ch.Title,('" ') ,'writing...' ])

str = ['ADD HERE details '];
para = Paragraph(str);
para.Style = {HAlign('justify')};

add(ch,para);

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
          
str = strcat('The engine takeoff power is ', ...
              strcat(' - ', num2str(takeoff_power,4)), ...
              Aircraft.Engine.Takeoff.Power.Attributes.unit,...
              ' at ', strcat(' - ', num2str(takeoff_engine_rpm,4)), Aircraft.Engine.Takeoff.RPM.Attributes.unit, '.', ...
              ' The rotational speed of the propeller is ',  strcat(' - ', num2str(takeoff_engine_rpm,4)), ...
              '/', num2str(reduction_ratio,4), ' = ', num2str(propeller_rotational_speed,4), ... 
              ' ', Aircraft.Engine.Takeoff.Propeller_rotational_speed.Attributes.unit, '.', ...
              ' The maximum continuous power is ', strcat(' - ', num2str(max_continuous_power,4)), ...
              Aircraft.Engine.Max_Continous.Power.Attributes.unit, ...
              '.', ' The mean engine torque is ',  ... 
              strcat(' - ', num2str(mean_engine_torque_max_power,4)), ...
              Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Takeoff.Mean_torque.Attributes.unit, '.', ...
              ' Using a factor of ', strcat(' - ', num2str(Aircraft.Engine.Correction_factor.value,4)), ' for a four cylinder engine, ', ...
              ' the limit torque will be ', ...
              strcat(' - ', num2str(limit_engine_torque_max_power,4)), ...
              Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Takeoff.Limit_torque.Attributes.unit, '.', ...
              ' This limit torque acts simultaneously with the 75 % of the inertia limit load. ', ...
              ' The mean engine torque at max continuous power is ', ...
              strcat(' - ', num2str(mean_engine_torque_max_continuous,4)), ' ', ...
              Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Max_Continous.Mean_torque.Attributes.unit, '.', ...
              ' Using a factor of ', strcat(' - ', num2str(Aircraft.Engine.Correction_factor.value,4)), ' for a four cylinder engine, ', ...
              ' the limit torque will be ', ...
              strcat(' - ', num2str(limit_engine_torque_max_continuous,4)), ...
              ' ', Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Max_Continous.Limit_torque.Attributes.unit, ...
              ' which acts simultaneously with the 100 % of the inertia limit load.');
        % MAX CONTINUOUS  
        myNumEq = strcat( ' = ', ...
                                num2str(max_continuous_power,4), ' * \frac{1000}', ...
                                '{ \frac{2 * 3.14 * ', num2str(propeller_rotational_speed,4), ...
                                '}{60}} =', num2str(mean_engine_torque_max_continuous), ...
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
                          ' = ', num2str(limit_engine_torque_max_continuous), ...
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
                                '}{60}} =', num2str(mean_engine_torque_max_power), ...
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
                          ' = ', num2str(limit_engine_torque_max_power), ...
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

para = Paragraph(str);
add(sec,para);
add(ch,sec);

Engine_mount_mass       = Aircraft.Engine.Engine_mount_mass.value;
Engine_accessories_mass = Aircraft.Engine.Engine_accessories_mass.value;
Propeller_spinner_mass  = Aircraft.Propeller.Propeller_spinner_mass.value;
Engine_block_mass       = Engine_mount_mass + Engine_accessories_mass + Propeller_spinner_mass;
g                       = Aircraft.Constants.g.value;
Total_side_load         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Side_loads.Total_side_load.value;

%sec
sec = Section();
sec.Title = 'Side load on engine mount';
str = strcat(' The limit load factor in a lateral direction is 1.33. The mass ', ...
       ' of the engine group is ', num2str(Engine_block_mass,4), ...
       Aircraft.Engine.Engine_block_mass.Attributes.unit, '.', ...
       ' The side load results is 1.33 * ', num2str(Engine_block_mass,4), ...
       ' * ', num2str(g,4), ' * (1/10) = ', num2str(Total_side_load,4), ...
       Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Side_loads.Total_side_load.Attributes.unit);
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);
%add(ch,para);

add(ch,sec);


Gust_limit_load              = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.nA1.value;
Inertia_load_on_engine_mount = Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Inertia_load_on_engine_mount.value;

%sec
sec = Section();
sec.Title = 'Intertia load on engine mount';

str = strcat('The inertia load is equal to the maximum limit load factor ', ...
    'times the engine group weight: ', num2str(Gust_limit_load,4), ...
    '*', num2str(Engine_block_mass,4), '*', num2str(g,4), '*(1/10) = ', ...
    num2str(Inertia_load_on_engine_mount,4), ...
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Inertia_load_on_engine_mount.Attributes.unit);

para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);

add(ch,sec);

%sec
sec = Section();
sec.Title = 'Gyroscopic loads';
str = ['ADD HERE details '];
para = Paragraph(str);
add(ch,para);

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