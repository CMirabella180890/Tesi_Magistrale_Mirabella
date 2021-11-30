function RepGenMain(Aircraft)
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

%%EXAMPLES (figure, table, equation)
% figure
% cross reference in para
% append(para,InternalLink('wingTableRef','Ref:wing'));
% add(sec,para)
% fig = FormalImage([PATH,'Polars.pdf']);
%          fig.Caption = 'Wing-Body reference Aerodynamics';
%          fig.Height = '5in';
%          fig.LinkTarget='polars_wb';
% % ADD TO SECTION
%          add(sec,fig);
% %ADD TO CHAPTER
%          add(ch,sec);
% 

% % table
% 3 COULUMNS TABLE OF WING PARAMETERS
%         wing = [wing, fieldValue, fieldUnit];
%         header = {'Wing parameters', 'Value', 'Measure unit'};
%
%          tbl = FormalTable(header,wing);
%         tbl.Style = {...
%             RowSep('solid','black','2px'),...
%             ColSep('solid','black','1px'),...
%             Border('ridge','black','5px')...
%             };
%         tbl.Header.Style = {...
%             RowSep('solid','black','2px'),...
%             ColSep('solid','black','2px'),...
%             Border('ridge','black','5px')...
%             };
%         tbl.TableEntriesStyle = {HAlign('left')};
%         % In order to put a table with a caption, the API Report denomination should
%         % be used, the other options are from API DOM. In order to solve the problem,
%         % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).      
%         tbl = BaseTable(tbl);
%         tbl.Title = 'Wing Geometrical Parameters';
%         tbl.LinkTarget = 'wingTableRef';
%ADD TO SECTION         
%         add(sec,tbl);
% ADD TO CHAPTER
% add(ch,sec);

% equation
% latex interprete with $ simbol
% eq = Equation("$ V_{S} = \sqrt{\frac{2 W_{MTOM}}{\rho_{0}C_{L_{MAX_{Clean}}}S}}");
% eq.DisplayInline = true;
% eq.FontSize = 12;
% eqImg = getImpl(eq,rpt);
% if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
%     eqImg.Style = {VerticalAlign("-30%")};
% elseif(rpt.Type == "docx") 
%     eqImg.Style = {VerticalAlign("-5pt")};
% end
% append(sec,eqImg);

% %List ordered/unordered
% item1 = 'ASTM F2245-12d,” ASTM.“ASTM F2245-12d, ASTM.';
% item2 = 'ABCD-FL-57-00 Wing Load Calculation, EASA.';
% item3 = 'ISO 2533:1975, International Standardization Organization, 1975.';
% item4 = 'CS-LSA Certification Specifications and Acceptable Means of Compliance, Amnd.1 29.Jul.2013, EASA, 2013.';
% item5 = '“ABCD-FTR-01-00 Flight Test Report,” EASA.';
% item6 = 'L. Smith, “NACA technical note 1945, ‘Aerodynamic characteristics of 15 NACA airfoil sections at seven Reynolds numbers from 0.7x10E6 to 9x10E6,” 1949.';
% item7 = 'ABCD-WB-08-00 Weight and Balance Report, EASA.';
% 
% ol = OrderedList({item1, item2, item3,...
%     item4, item5,item6,item7});
% 
% % ol = UnorderedList({item1, item2, item3,...
%     item4, item5,item6,item7});
% append(ch,ol);



%% Report set-up
% Flight loads user's settings
Aircraft.Report.author = 'Pierluigi Della Vecchia';
str = strcat("The author is " , " ", Aircraft.Report.author);
disp(str)
disp(' ')
Aircraft.Report.society = 'Design of Aircraft and Flight Technologies, DAF';
str = strcat("The author affiliation is " , " ", Aircraft.Report.society);
disp(str)
%Template available
% 'DAF_template'
% 'DAF_SMUP_template'
% 'SMUP_template'
Aircraft.Report.template = 'No Template';
%Mytemplate.template = 'DAF_template'; % Other availabe
str = strcat("The report template is " , " ",Aircraft.Report.template);
disp(str)

%report format
Aircraft.Report.Type = 'pdf';
str = strcat("The report format is " , " ",Aircraft.Report.Type);
disp(str)
%pdf
%docx
%html

%%%

rpt = Report('Flight Loads',Aircraft.Report.Type);

disp('Writing report...')
rpt.Locale = 'en';

%% Flight Loads Report Generator
% INPUT:  Aircraft data structure
% OUTPUT: FLIGHT LOADS.xxx report

RepDir = pwd;

%% TITLE PAGE
%Title
tp = TitlePage();
% if any(Phi)
%     tp.Title = 'Conceptual Design of Hybrid Electric Aircraft';
% elseif any(phi)
%     tp.Title = 'Conceptual Design of Turbo-Electric Aircraft';
% else
Title = ['Flight Loads: ', ' ',Aircraft.Certification.Aircraft_Name.value,' aircraft'];
tp.Title = Title;
%end
%% Cover image
%tp.Image = '_figures/cover.jpg';
cd ..
geometry_path = [pwd '/Geometry/'];
tp.Image = [geometry_path 'Aircraft3D.png'];
tp.Publisher = Aircraft.Report.society;
tp.Author = Aircraft.Report.author;
tp.PubDate = date();
cd (RepDir);

add(rpt,tp);


%% TABLE OF CONTENT
toc = TableOfContents();
add(rpt,toc);
% List of figures
lof = ListOfFigures();
lof.Title = "List of Figures";
append(rpt,lof);
% List of tables
lot = ListOfTables();
lot.Title = "List of Tables";
append(rpt,lot);

%% CHAPTERS
RepGenCh1       %Introduction
RepGenCh2       %References
RepGenCh3       %List of Abbreviations
RepGenCh4       %Aircraft data
RepGenCh5       %Design Airspeeds
RepGenCh6       %Altitude
RepGenCh7       %Manoeuvring and Gust load factors n
RepGenCh8       %V-n Envelope
RepGenCh9       %Loads on the aeroplane
RepGenCh10      %Loads on the wing
RepGenCh11      %Loads on the horizontal tail
RepGenCh12      %Loads on the vertical tail
RepGenCh13      %Loads on the wing flaps
RepGenCh14      %Loads on the control surfaces
RepGenCh15      %Power plant



%% END OF REPORT
close(rpt)
rptview(rpt)
disp('Report successfully written. Execution terminated.')
end
