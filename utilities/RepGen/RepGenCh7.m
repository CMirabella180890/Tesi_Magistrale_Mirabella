%% CHAPTER 7 - Manoeuvring and Gust load factors n
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

ch7 = Chapter();
ch7.Title = 'Manoeuvring and Gust load factors n';

str = ['ADD HERE Manoeuvring and Gust load factors n, figures, tables....ecc. ecc. '];
para = Paragraph(str);

%moving to another path for figure
cd ..
cd ..
 regulation = Aircraft.Certification.Regulation.value;
 results_path = [pwd '\' regulation '\Output\'];

 cd (RepDir);

fig = FormalImage([results_path,'Gustenvelope.pdf']);
         fig.Caption = 'Maneuver and Gust load factors and diagram';
         fig.Height = '5in';
         fig.LinkTarget='maneuver_ref';
         add(ch7,fig);

%add(ch7,sec1);



add(ch7,para)
%% END chapter
%Adding chapters
add(rpt,ch7);