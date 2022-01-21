%% CHAPTER 11 - Loads on the htail
%Ref EASA: WING LOAD CALCULATION
%(Example document for LSA applicants – v1 of 08.03.16)
%Date of issue: DD/MM/YYYY
%Document reference: ABCD-FL-57-00

import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

% chapter_number = 11;
% ch = strcat('ch' , num2str(chapter_number)); 
ch = Chapter();
ch.Title = 'Loads on the horizontal tail';
disp(['Chapter 11', (' "'), ch.Title,('" ') ,'writing...' ])

str = ['ADD HERE details '];
para = Paragraph(str);
add(ch,para);


%sec
sec = Section();
sec.Title = 'Balancing loads';
str = ['ADD HERE details '];
para = Paragraph(str);
add(ch,para);

fig = FormalImage([results_path,'Balancingloads.png']);
fig.Caption = 'Balancing loads';
fig.Height = '5in';
fig.LinkTarget='bala_loads';
add(sec,fig);

add(ch,sec);

%sec
sec = Section();
sec.Title = 'Manouevring loads';
str = ['ADD HERE details '];
para = Paragraph(str);
add(ch,para);
%sub
subsec = Section();
subsec.Title = 'Unchecked manoeuvre';

add(sec,subsec);

%sub
subsec = Section();
subsec.Title = 'Checked manoeuvre';

add(sec,subsec);

%sub
subsec = Section();
subsec.Title = 'Gust loads';

add(sec,subsec);


add(ch,sec);

%sec
sec = Section();
sec.Title = 'Horizontal tail loads summary';
str = ['ADD HERE details '];
para = Paragraph(str);
add(ch,para);

add(ch,sec);


%sec
sec = Section();
sec.Title = 'Unsysmmetrical loads';
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