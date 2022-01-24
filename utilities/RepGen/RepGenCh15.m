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
add(ch,para);

%sec
sec = Section();
sec.Title = 'Engine torque';

str = strcat('The engine takeoff power is ', ...
              strcat(' - ', num2str(Aircraft.Engine.Takeoff.Power.value)), ...
              Aircraft.Engine.Takeoff.Power.Attributes.unit,...
              ' at ', strcat(' - ', num2str(Aircraft.Engine.Takeoff.RPM.value)), Aircraft.Engine.Takeoff.RPM.Attributes.unit, '.', ...
              ' The rotational speed of the propeller is ',  strcat(' - ', num2str(Aircraft.Engine.Takeoff.RPM.value)), ...
              '/', num2str(Aircraft.Engine.Reduction_ratio.value), ' = ', num2str(Aircraft.Engine.Takeoff.Propeller_rotational_speed.value), ... 
              ' ', Aircraft.Engine.Takeoff.Propeller_rotational_speed.Attributes.unit, '.', ...
              ' The maximum continuous power is ', strcat(' - ', num2str(Aircraft.Engine.Max_Continous.Power.value)), ...
              Aircraft.Engine.Max_Continous.Power.Attributes.unit, ...
              '.', ' The mean engine torque is ',  ... 
              strcat(' - ', num2str(Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Takeoff.Mean_torque.value)), ...
              Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Takeoff.Mean_torque.Attributes.unit, '.', ...
              ' Using a factor of ', strcat(' - ', num2str(Aircraft.Engine.Correction_factor.value)), ' for a four cylinder engine, ', ...
              ' the limit torque will be ', ...
              strcat(' - ', num2str(Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Takeoff.Limit_torque.value)), ...
              Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Takeoff.Limit_torque.Attributes.unit, '.', ...
              ' This limit torque acts simultaneously with the 75 % of the inertia limit load. ', ...
              ' The mean engine torque at max continuous power is ', ...
              strcat(' - ', num2str(Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Max_Continous.Mean_torque.value)), ' ', ...
              Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Max_Continous.Mean_torque.Attributes.unit, '.', ...
              ' Using a factor of ', strcat(' - ', num2str(Aircraft.Engine.Correction_factor.value)), ' for a four cylinder engine, ', ...
              ' the limit torque will be ', ...
              strcat(' - ', num2str(Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Max_Continous.Limit_torque.value)), ...
              ' ', Aircraft.Certification.Regulation.SubpartC.Flightloads.Engine_loads.Max_Continous.Limit_torque.Attributes.unit, ...
              ' which acts simultaneously with the 100 % of the inertia limit load.');

% str = ['ADD HERE details '];

para = Paragraph(str);
add(sec,para);
add(ch,sec);

%sec
sec = Section();
sec.Title = 'Side load on engine mount';
str = ['ADD HERE details '];
para = Paragraph(str);
add(ch,para);

add(ch,sec);


%sec
sec = Section();
sec.Title = 'Intertia load on engine mount';
str = ['ADD HERE details '];
para = Paragraph(str);
add(ch,para);

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