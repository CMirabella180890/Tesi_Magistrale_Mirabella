%% CHAPTER 3 - List of Abbreviations
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

ch1 = Chapter();
ch1.Title = 'Introduction';

ch3 = Chapter();
ch3.Title = 'List of Abbreviations';

str = ['ADD HERE list of abbreviations as a formatted table....to be created'];
para = Paragraph(str);
% append(para,InternalLink('tlarTableRef','refTabella'));
add(ch3,para)
%% END chapter
%Adding chapters
add(rpt,ch3);