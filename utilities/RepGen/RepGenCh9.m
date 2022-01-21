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

% chapter_number = 9;
% ch = strcat('ch' , num2str(chapter_number)); 
ch = Chapter();
ch.Title = 'Loads on the aeroplane';
disp(['Chapter 9', (' "'), ch.Title,('" ') ,'writing...' ])

str = ['ADD HERE details for balancing Equation'];
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