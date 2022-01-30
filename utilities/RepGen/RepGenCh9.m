%% CHAPTER 9 - Loads on the aeroplane
%Ref EASA: WING LOAD CALCULATION
%(Example document for LSA applicants – v1 of 08.03.16)
%Date of issue: DD/MM/YYYY
%Document reference: ABCD-FL-57-00

import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

% (a) Strength requirements are specified in terms of limit loads (the maximum loads to be expected
% in service) and ultimate loads (limit loads multiplied by prescribed factors of safety). Unless
% otherwise provided, prescribed loads are limit loads.

% (b) Unless otherwise provided, the air, ground, and water loads must be placed in equilibrium with
% inertia forces, considering each item of mass in the aeroplane. These loads must be distributed
% to conservatively approximate or closely represent actual conditions.

% (c) If deflections under load would significantly change the distribution ofexternal or internal loads,
% this redistribution must be taken into account.

% (d) Simplified structural design criteria given in this Subpart C and its appendices may be used only
% for aeroplanes with conventional configurations. IfAppendixA is used, the entire appendix must
% be substituted for the corresponding paragraphs of this subpart, i.e. CS-VLA 321 to 459

% chapter_number = 9;
% ch = strcat('ch' , num2str(chapter_number)); 
ch = Chapter();
ch.Title = 'Loads on the aeroplane';
disp(['Chapter 9', (' "'), ch.Title,('" ') ,'writing...' ])

str = ['Strength requirements are specified in terms of limit loads' ...
    ' (the maximum loads to be expected in service) and ultimate loads' ...
    ' (limit loads multiplied by prescribed factors of safety). Unless' ...
    ' otherwise provided, prescribed loads are limit loads.' ...
    ' Unless otherwise provided, the air, ground, and water loads' ...
    ' must be placed in equilibrium with inertia forces, considering' ...
    ' each item of mass in the aeroplane. These loads must be' ...
    ' distributed to conservatively approximate or closely represent actual conditions.'];
para = Paragraph(str);
add(ch,para);

%sec
sec = Section();
sec.Title = 'Reference axes and sign convention';
str = ['ADD HERE details for balancing Equation'];
para = Paragraph(str);
add(ch,para);
subsec = Section();
subsec.Title = 'aaaaa';

add(sec,subsec);
add(ch,sec);

%sec
sec = Section();
sec.Title = 'Symmetrical flight conditions';
str = ['ADD HERE details for balancing Equation'];
para = Paragraph(str);
add(ch,para);
add(ch,sec);

%sec
sec = Section();
sec.Title = 'Aerodynamic centre';
str = ['ADD HERE details for balancing Equation'];
para = Paragraph(str);
add(ch,para);
add(ch,sec);

%sec
sec = Section();
sec.Title = 'Pitching moment of the wing';
str = ['ADD HERE details for balancing Equation'];
para = Paragraph(str);
add(ch,para);
add(ch,sec);


%moving to another path for figure
cd ..
cd ..
 regulation = Aircraft.Certification.Regulation.value;
 results_path = [pwd '\' regulation '\Output\'];

 cd (RepDir);
 
 fig = FormalImage([results_path,'Wingairloads.png']);
 fig.Caption = 'Wing airloads';
 fig.Height = '5in';
 fig.LinkTarget='wing_loads';
 add(sec,fig);

  fig = FormalImage([results_path,'Balancingloads.png']);
 fig.Caption = 'Balancing loads';
 fig.Height = '5in';
 fig.LinkTarget='bala_loads';
 add(sec,fig);





%% END chapter
%Adding chapters
add(rpt,ch);