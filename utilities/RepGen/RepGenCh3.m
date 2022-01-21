%% CHAPTER 3 - List of Abbreviations
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

ch = Chapter();
ch.Title = 'List of Abbreviations';
disp(['Chapter 3', (' "'), ch.Title,('" ') ,'writing...' ])

str = ['ADD HERE list of abbreviations as a formatted table....to be created'];

% %List un-ordered
 item1 = 'CL = lift coefficient';
 item2 = 'CD....';
 item3 = '...';
 item4 = '...';
 item5 = '...';
 item6 = '...';
 item7 = '...';

ol = UnorderedList({item1, item2, item3,...
    item4, item5,item6,item7});

append(ch,ol);

para = Paragraph(str);
% append(para,InternalLink('tlarTableRef','refTabella'));
add(ch,para)
%% END chapter
%Adding chapters
add(rpt,ch);