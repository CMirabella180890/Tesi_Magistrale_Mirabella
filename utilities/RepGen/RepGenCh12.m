<<<<<<< HEAD
%% CHAPTER 12 - Loads on the vertical tail
%Ref EASA/Tecnam p2006 report

import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

% chapter_number = 12;
% ch = strcat('ch' , num2str(chapter_number)); 
ch = Chapter();
ch.Title = 'Loads on the vertical tail';
disp(['Chapter 12', (' "'), ch.Title,('" ') ,'writing...' ])


str = ['ADD HERE details '];
para = Paragraph(str);
add(ch,para);


% %sec
% sec = Section();
% sec.Title = 'Balancing loads';
% str = ['ADD HERE details '];
% para = Paragraph(str);
% add(ch,para);
% add(ch,sec);

%sec
sec = Section();
sec.Title = 'Manouevring loads';
str = ['ADD HERE details '];
para = Paragraph(str);
add(ch,para);
%sub
subsec = Section();
subsec.Title = 'a(1)';

add(sec,subsec);

%sub
subsec = Section();
subsec.Title = 'a(2)';

add(sec,subsec);

%sub
subsec = Section();
subsec.Title = 'a(3)';

add(sec,subsec);

%sub
subsec = Section();
subsec.Title = 'Gust loads';

add(sec,subsec);

add(ch,sec);

%sec
sec = Section();
sec.Title = 'Vertical tail loads summary';
str = ['ADD HERE details '];
para = Paragraph(str);
add(ch,para);

add(ch,sec);

%sec
sec = Section();
sec.Title = 'Combined loads';
str = ['ADD HERE details on h-v combined loads'];
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
=======
%% CHAPTER 12 - Loads on the vertical tail
%Ref EASA/Tecnam p2006 report

import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

% chapter_number = 12;
% ch = strcat('ch' , num2str(chapter_number)); 
ch = Chapter();
ch.Title = 'Loads on the vertical tail';
disp(['Chapter 12', (' "'), ch.Title,('" ') ,'writing...' ])


str = ['ADD HERE details '];
para = Paragraph(str);
add(ch,para);


% %sec
% sec = Section();
% sec.Title = 'Balancing loads';
% str = ['ADD HERE details '];
% para = Paragraph(str);
% add(ch,para);
% add(ch,sec);

%sec
sec = Section();
sec.Title = 'Manouevring loads';
str = ['ADD HERE details '];
para = Paragraph(str);
add(ch,para);
%sub
subsec = Section();
subsec.Title = 'a(1)';

add(sec,subsec);

%sub
subsec = Section();
subsec.Title = 'a(2)';

add(sec,subsec);

%sub
subsec = Section();
subsec.Title = 'a(3)';

add(sec,subsec);

%sub
subsec = Section();
subsec.Title = 'Gust loads';

add(sec,subsec);

add(ch,sec);

%sec
sec = Section();
sec.Title = 'Vertical tail loads summary';
str = ['ADD HERE details '];
para = Paragraph(str);
add(ch,para);

add(ch,sec);

%sec
sec = Section();
sec.Title = 'Combined loads';
str = ['ADD HERE details on h-v combined loads'];
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
>>>>>>> eb25f8fea66ec6552ebf7b0bde47154c0a1a06a0
add(rpt,ch);