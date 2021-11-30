%% CHAPTER 6 - Altitude
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

ch = Chapter();
ch.Title = 'Altitude';

str = strcat('The maximum permissible operational altitude is ADD H13000ft.',...
            'Despite the CS-LSA requirements do not require to accounts for the effects of altitude, ',... 
    'such effects have been considered up to 10000 ft. In fact the gust load factor have been ',... 
    'calculated at such altitude. This is considered acceptable since it covers the operational ',...
    'range within which the aeroplane will fly most of the time.');
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(ch,para)

str = strcat('(Note: the CS-LSA requirement does not require to account for the effects of altitude. ',...
            'Calculating the loads at sea level would be acceptable.',... 
            'In this case, the choice to consider such effect up to 10000 ft is a decision of a designer, which would be accepted by the team.)');
para = Paragraph(str);
para.Style = {HAlign('justify')};

add(ch,para)

%% END chapter
%Adding chapters
add(rpt,ch);