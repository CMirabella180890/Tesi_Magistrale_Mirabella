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
ch.Title = 'Loads on the wing';

str = ['ADD HERE details for balancing Equation'];
para = Paragraph(str);
add(ch,para);


%sec
sec = Section();
sec.Title = 'Influence of the fuselage';
str = ['ADD HERE details for fuselage effect how are they accounted?'];
para = Paragraph(str);
add(ch,para);
add(ch,sec);

%sec
sec = Section();
sec.Title = 'Forces and moments acting on the wings';
str = ['ADD HERE details for balancing Equation'];
para = Paragraph(str);
add(ch,para);
%sub
subsec = Section();
subsec.Title = 'SpanWise Airloads Distribution';

add(sec,subsec);

%sub
subsec = Section();
subsec.Title = 'Normal and parallel component';

add(sec,subsec);

%sub
subsec = Section();
subsec.Title = 'Shear, Bending and Torsion';

add(sec,subsec);

%sub
subsec = Section();
subsec.Title = 'Critical load condition';

add(sec,subsec);
add(ch,sec);

%sec
sec = Section();
sec.Title = 'Unsymmetrical loads';
str = ['ADD HERE details for uns loads'];
para = Paragraph(str);
add(ch,para);
%sub
subsec = Section();
subsec.Title = 'Rolling condition';
add(sec,subsec);

%sub
subsec = Section();
subsec.Title = 'Effect of aileron displacement on the wing torsion';
add(sec,subsec);

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