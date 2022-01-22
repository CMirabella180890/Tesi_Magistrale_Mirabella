<<<<<<< HEAD
%% CHAPTER 2 - References
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

ch = Chapter();
ch.Title = 'References';
disp(['Chapter 2', (' "'), ch.Title,('" ') ,'writing...' ])

str = ['HERE BELOW AN EXAMPLE OF REFERENCES TO BE EDITED'];
para = Paragraph(str);

%references: CAN BE MOVED IN OTHER position
ref1 = 'ASTM F2245-12d,” ASTM.“ASTM F2245-12d, ASTM.';
ref2 = 'ABCD-FL-57-00 Wing Load Calculation, EASA.';
ref3 = 'ISO 2533:1975, International Standardization Organization, 1975.';
ref4 = 'CS-LSA Certification Specifications and Acceptable Means of Compliance, Amnd.1 29.Jul.2013, EASA, 2013.';
ref5 = '“ABCD-FTR-01-00 Flight Test Report,” EASA.';
ref6 = 'L. Smith, “NACA technical note 1945, ‘Aerodynamic characteristics of 15 NACA airfoil sections at seven Reynolds numbers from 0.7x10E6 to 9x10E6,” 1949.';
ref7 = 'ABCD-WB-08-00 Weight and Balance Report, EASA.';

ol = OrderedList({ref1, ref2, ref3,...
    ref4, ref5,ref6,ref7});

append(ch,ol);

% append(para,InternalLink('tlarTableRef','refTabella'));
add(ch,para)

%% END chapter
%Adding chapters
=======
%% CHAPTER 2 - References
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

ch = Chapter();
ch.Title = 'References';
disp(['Chapter 2', (' "'), ch.Title,('" ') ,'writing...' ])

str = ['HERE BELOW AN EXAMPLE OF REFERENCES TO BE EDITED'];
para = Paragraph(str);

%references: CAN BE MOVED IN OTHER position
ref1 = 'ASTM F2245-12d,” ASTM.“ASTM F2245-12d, ASTM.';
ref2 = 'ABCD-FL-57-00 Wing Load Calculation, EASA.';
ref3 = 'ISO 2533:1975, International Standardization Organization, 1975.';
ref4 = 'CS-LSA Certification Specifications and Acceptable Means of Compliance, Amnd.1 29.Jul.2013, EASA, 2013.';
ref5 = '“ABCD-FTR-01-00 Flight Test Report,” EASA.';
ref6 = 'L. Smith, “NACA technical note 1945, ‘Aerodynamic characteristics of 15 NACA airfoil sections at seven Reynolds numbers from 0.7x10E6 to 9x10E6,” 1949.';
ref7 = 'ABCD-WB-08-00 Weight and Balance Report, EASA.';

ol = OrderedList({ref1, ref2, ref3,...
    ref4, ref5,ref6,ref7});

append(ch,ol);

% append(para,InternalLink('tlarTableRef','refTabella'));
add(ch,para)

%% END chapter
%Adding chapters
>>>>>>> eb25f8fea66ec6552ebf7b0bde47154c0a1a06a0
add(rpt,ch);