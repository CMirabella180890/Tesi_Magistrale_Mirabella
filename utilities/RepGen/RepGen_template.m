function RepGen(Aircraft)
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

%% Report details
% Propeller Tool user's settings (to be moved)
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
Aircraft.Report.template = 'DAF_template';
%Mytemplate.template = 'DAF_template'; % Other availabe
str = strcat("The report template is " , " ",Aircraft.Report.template);
disp(str)

%%%

rpt = Report('Flight Loads','pdf',Aircraft.Report.template);
%pdf --> MYTEMPLATE.pdftx
disp('Writing report...')
rpt.Locale = 'en';

%% Propeller Report Generator
% INPUT:  Aircraft data structure
% OUTPUT: FLIGHT LOADS .pdf report

RepDir = pwd;

%% TITLE PAGE
%Title
tp = TitlePage('TemplateSrc',Aircraft.Report.template,'TemplateName','TitlePage');
% if any(Phi)
%     tp.Title = 'Conceptual Design of Hybrid Electric Aircraft';
% elseif any(phi)
%     tp.Title = 'Conceptual Design of Turbo-Electric Aircraft';
% else
Title = ['Flight Loads: ', ' ',Aircraft.Certification.Aircraft_Name.value,' aircraft'];
tp.Title = Title;
%end
%tp.Image = '_figures/cover.jpg';
cd ..
results_path = [pwd '/Geometry/DroneVLA_results/'];
tp.Image = [results_path 'Aircraft3D.png'];
tp.Publisher = Aircraft.Report.society;
tp.Author = Aircraft.Report.author;
tp.PubDate = date();
cd (RepDir);

add(rpt,tp);


%% TABLE OF CONTENT
toc = TableOfContents('TemplateSrc',Aircraft.Report.template,'TemplateName','TableOfContents');
add(rpt,toc);

lof = ListOfFigures();
lof.Title = "List of Figures";
append(rpt,lof);

lot = ListOfTables();
lot.Title = "List of Tables";
append(rpt,lot);

% switch prop.case_to_do
%     case {1,4}
%% CHAPTER 1 - INTRODUCTION
ch1 = Chapter('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
ch1.Title = 'Introduction';

para = Paragraph(strcat('This document defines the SUBPART C - Structure - Flight Loads of the:',' ', char(Aircraft.Certification.Aircraft_Name.value),...
    'The boundaries of the flight envelope will be defined within this document. All speeds are calibrated airspeeds (CAS) (requirement 4.4 [1])',...
    'and given in knots if not stated otherwise.',...
    'All other units used are metric (SI units).',...
    'The weights are given in mass units (kg) but the formulas require force units as input,',...
    'therefore these are calculated in place wherever they are used.', ...
    'Note: The speeds defined within this document should be used for the placards,',...
    'speed markings, aeroplane flight manual (limitations), load calculations and need to be verified by flight test.'));
% append(para,InternalLink('tlarTableRef','refTabella'));
add(ch1,para)
%% END chapter
%Adding chapters
add(rpt,ch1);

%% CHAPTER 2 - References
ch2 = Chapter('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
ch2.Title = 'References';

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


append(ch2,ol);

% append(para,InternalLink('tlarTableRef','refTabella'));
add(ch2,para)
%% END chapter
%Adding chapters
add(rpt,ch2);


%% CHAPTER 3 - List of Abbreviations
ch3 = Chapter('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
ch3.Title = 'List of Abbreviations';

str = ['ADD HERE list of abbreviations as a formatted table....to be created'];
para = Paragraph(str);
% append(para,InternalLink('tlarTableRef','refTabella'));
add(ch3,para)
%% END chapter
%Adding chapters
add(rpt,ch3);


%% CHAPTER 4 - Aircraft data
ch4 = Chapter('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
ch4.Title = 'Aircraft data';

str = ['Add here all the aircraft geometrical, aero and inertial and masses data useful for following paragraph'];
para = Paragraph(str);
% append(para,InternalLink('tlarTableRef','refTabella'));
add(ch4,para)

%         %% CHAPTER 4 - SECTION 1
%         % geo
sec1 = Section('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
sec1.Title = 'Geometry';

para = Paragraph('The aircraft reference geometry is summarized in table:');
 append(para,InternalLink('wingTableRef','refTabella'));
add(sec1,para)
         wing = fieldnames(Aircraft.Geometry.Wing);
         fieldValue = cell(length(wing),1);
         fieldUnit = cell(length(wing),1);
         for i = 1:length(wing)
             fieldValue{i} = Aircraft.Geometry.Wing.(wing{i}).value;
             % significant digits
             if isnumeric(fieldValue{i})
                 fieldValue{i} = num2str(fieldValue{i},5);
             end
             % Not every field has Attributes. If not, they have only one field
             if length(fieldnames(Aircraft.Geometry.Wing.(wing{i}))) > 1
                 fieldUnit{i} = Aircraft.Geometry.Wing.(wing{i}).Attributes.unit;
             else % Void cells cannot be converted in strings
                 fieldUnit{i} = '-';
             end
         end
         wing = [wing, fieldValue, fieldUnit];
         header = {'Wing parameters', 'Value', 'Measure unit'};
         tbl = FormalTable(header,wing);
        tbl.Style = {...
            RowSep('solid','black','2px'),...
            ColSep('solid','black','1px'),...
            Border('ridge','black','5px')...
            };
        tbl.Header.Style = {...
            RowSep('solid','black','2px'),...
            ColSep('solid','black','2px'),...
            Border('ridge','black','5px')...
            };
        tbl.TableEntriesStyle = {HAlign('left')};
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
        tbl = BaseTable(tbl);
        tbl.Title = 'Wing Geometrical Parameters';
        tbl.LinkTarget = 'wingTableRef';

        add(sec1,tbl);
%


add(ch4,sec1);

%         %% CHAPTER 4 - SECTION 2
%         % aero
sec2 = Section('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
sec2.Title = 'Aerodynamic';

para = Paragraph('The aircraft reference aerodynamic is in table:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec2,para)



add(ch4,sec2);



%% END chapter
%Adding chapters
add(rpt,ch4);

%% CHAPTER 5 - Design Airspeeds
ch5 = Chapter('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
ch5.Title = 'Design Airspeeds';

str = ['This chapter defines the operating and design airspeeds as required for certification REFFFF'];
para = Paragraph(str);
% append(para,InternalLink('tlarTableRef','refTabella'));
add(ch5,para)

%         %% CHAPTER 5 - SECTION 1
%         % VH
sec1 = Section('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
sec1.Title = 'Maximum speed in level flight VH';

para = Paragraph(strcat('According to flight tests [5] at maximum weight and maximum continuous',...
    'power at sea level conditions, the maximum speed in level flight has been determined:',...
    'V_H=130\ kts'));
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec1,para)

add(ch5,sec1);

%         %% CHAPTER 5 - SECTION 2
%         % VS, VS0, VS1
sec2 = Section('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
sec2.Title = 'Stall speeds VS, VS0, VS1';

para = Paragraph('ADD TEXTS:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec2,para)

%Equation Example
eq = Equation("\int_{0}^{2} x^2\sin(x) dx");
eq.DisplayInline = true;
eq.FontSize = 14;
eqImg = getImpl(eq,rpt);
if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
    eqImg.Style = {VerticalAlign("-30%")};
elseif(rpt.Type == "docx") 
    eqImg.Style = {VerticalAlign("-5pt")};
end
append(sec2,eqImg);

add(ch5,sec2);

%         %% CHAPTER 5 - SECTION 3
%         % VA
sec3 = Section('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
sec3.Title = 'Design manoeuvring speed VA         ';

para = Paragraph('ADD TEXTS:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec3,para)



add(ch5,sec3);


%         %% CHAPTER 5 - SECTION 4
%         %
sec4 = Section('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
sec4.Title = 'Flaps maximum operating speed VF';

para = Paragraph('ADD TEXTS:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec4,para)


add(ch5,sec4);

%         %% CHAPTER 5 - SECTION 5
%         %
sec5 = Section('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
sec5.Title = 'Flaps maximum operating speed VFE';

para = Paragraph('ADD TEXTS:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec5,para)



add(ch5,sec5);

%         %% CHAPTER 5 - SECTION 6
%         %
sec6 = Section('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
sec6.Title = 'Design cruising speed VC';

para = Paragraph('ADD TEXTS:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec6,para)


add(ch5,sec6);


%         %% CHAPTER 5 - SECTION 7
%         %
sec7 = Section('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
sec7.Title = 'Design dive speed VD';

para = Paragraph('ADD TEXTS:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec7,para)


add(ch5,sec7);


%         %% CHAPTER 5 - SECTION 8
%         %
sec8 = Section('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
sec8.Title = 'Demonstrated dive speed VDF';

para = Paragraph('ADD TEXTS:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec8,para)


add(ch5,sec8);

%         %% CHAPTER 5 - SECTION 9
%         %
sec9 = Section('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
sec9.Title = 'Never exceed speed VNE';

para = Paragraph('ADD TEXTS:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec9,para)


add(ch5,sec9);



%% END chapter
%Adding chapters
add(rpt,ch5);

%% CHAPTER 6 - Altitude
ch6 = Chapter('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
ch6.Title = 'Altitude';

str = ['ADD HERE ALTITUDE DETAILS'];
para = Paragraph(str);
% append(para,InternalLink('tlarTableRef','refTabella'));
add(ch6,para)
%% END chapter
%Adding chapters
add(rpt,ch6);


%% CHAPTER 7 - Manoeuvring and Gust load factors n
ch7 = Chapter('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
ch7.Title = 'Manoeuvring and Gust load factors n';

str = ['ADD HERE Manoeuvring and Gust load factors n, figures, tables....ecc. ecc. '];
para = Paragraph(str);
% append(para,InternalLink('tlarTableRef','refTabella'));
add(ch7,para)
%% END chapter
%Adding chapters
add(rpt,ch7);


%% CHAPTER 8 - Manoeuvring and Gust load factors n
ch8 = Chapter('TemplateSrc',Aircraft.Report.template,'TemplateName','Section');
ch8.Title = 'V-n Envelope';

str = ['ADD HERE V-n Envelope'];
para = Paragraph(str);
% append(para,InternalLink('tlarTableRef','refTabella'));
add(ch8,para)
%% END chapter
%Adding chapters
add(rpt,ch8);


%
%         tlpr = fieldnames(prop.TLPR);
%         fieldValue = cell(length(tlpr),1);
%         fieldUnit = cell(length(tlpr),1);
%         for i = 1:length(tlpr)
%             fieldValue{i} = prop.TLPR.(tlpr{i}).value;
%             % significant digits
%             if isnumeric(fieldValue{i})
%                 fieldValue{i} = num2str(fieldValue{i},5);
%             end
%             % Not every field has Attributes. If not, they have only one field
%             if length(fieldnames(prop.TLPR.(tlpr{i}))) > 1
%                 fieldUnit{i} = prop.TLPR.(tlpr{i}).Attributes.unit;
%             else % Void cells cannot be converted in strings
%                 fieldUnit{i} = '-';
%             end
%         end
%         tlpr = [tlpr, fieldValue, fieldUnit];
%         header = {'TLPR', 'Value', 'Measure unit'};
%         tbl = FormalTable(header,tlpr);
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
%         tbl.Title = 'Top Level Propeller Requirements';
%         tbl.LinkTarget = 'tlarTableRef';
%
%         add(sec1,tbl);
%
%         % blank line
%         para = Paragraph();
%         append(para, Text(sprintf(' \n ')));
%         add(sec1,para)
%         % end blank line
%
%         if prop.Thrust == 0
%             para2 = Paragraph('The propeller design is by Power');
%             add(sec1,para2)
%             para3 = Paragraph(['The Power value is' ' ' num2str(prop.TLPR.Power.value) ' ' prop.TLPR.Power.Attributes.unit]);
%             add(sec1,para3)
%         else
%             para2 = Paragraph('The propeller design is by Thrust');
%             add(sec1,para2)
%             para3 = Paragraph(['The Thrust value is '  ' ' num2str(prop.TLPR.Thrust.value) ' ' prop.TLPR.Thrust.Attributes.unit]);
%             add(sec1,para3)
%         end
%         % % append(para,InternalLink('tlarTableRef','refTabella'));
%
%         % blank line
%         para = Paragraph();
%         append(para, Text(sprintf(' \n ')));
%         add(sec1,para)
%         % end blank line
%
%         % J=V/nD
%         para4 = Paragraph(['The propeller design point will be at J=V/nD equal to ', num2str(prop.J.value)]);
%         add(sec1,para4)
%
%         %propeller MODE of design
%         %prop.mode       = 'POT';propeller design and analysis mode:
%         % GRAD -> Graded Momentum Formulation
%         % POT  -> Potential formulation
%         % VRTX -> Vortex formulation
%
%         if prop.mode == "GRAD"
%             para5 = Paragraph('The propeller design method is "Graded Momentum Formulation"' );
%             add(sec1,para5)
%         else
%             if prop.mode == "POT"
%                 para5 = Paragraph('The propeller design method is "Potential formulation"' );
%                 add(sec1,para5)
%             else
%                 para5 = Paragraph('The propeller design method is "Vortex formulation"' );
%                 add(sec1,para5)
%             end
%         end
%

%
%         %% CHAPTER 1 - SECTION 2
%         sec2 = Section('TemplateSrc',prop.template,'TemplateName','Section');
%         sec2.Title = 'Blade geometry';
%         para = Paragraph(...
%             ['The main designed blade geometrical characteristics are shown in the following two figures \n',...
%             'The geometry of the propeller is defined according to mininum induced load method. ']);
%         add(sec2,para)
%
%         % Propeller geometrical characteristics c/R and pitch angle
%         fig = FormalImage([results_path,'PropChords.png']);
%         fig.Caption = 'Propeller chords law';
%         fig.Height = '5in';
%         fig.LinkTarget='PropChords';
%         add(sec2,fig);
%         %
%         fig = FormalImage([results_path,'PropPitch.png']);
%         fig.Caption = 'Propeller pitch angle law';
%         fig.Height = '5in';
%         fig.LinkTarget='PropPitch';
%         add(sec2,fig)
%         %
%
%         add(ch1,sec2)
%         % % %add(rpt,ch1) % moved to the end of chapter 2 to include a table of chap2
%         % % %end chapter 1


%
%                %% CHAPTER  - Prop variable pitch analysis
%         ch4 = Chapter('TemplateSrc',prop.template,'TemplateName','Section');
%         ch4.Title = 'Propeller CAD';
%
%         para = Paragraph(strcat('This chapter contains the CAD rendering of designed propeller'));
%         % append(para,InternalLink('tlarTableRef','refTabella'));
%         add(ch4,para)
%
%         sec1 = Section('TemplateSrc',prop.template,'TemplateName','Section');
%         sec1.Title = 'CAD shape';
%
%         % Propeller BLADE shape
%         fig = FormalImage([results_path,'PropBLADE.png']);
%         fig.Caption = 'Propeller BLADE';
%         fig.Height = '5in';
%         fig.LinkTarget='PropBlade';
%         add(sec1,fig);
%
%         % Propeller CAD shape
%         fig = FormalImage([results_path,'PropCAD.png']);
%         fig.Caption = 'Propeller CAD';
%         fig.Height = '5in';
%         fig.LinkTarget='PropCad';
%         add(sec1,fig);
%
%
%         add(ch4,sec1);


%         add(rpt,ch4);





%     %case {2,5}
%         %% CHAPTER 1 - TLPR
%         ch1 = Chapter('TemplateSrc',prop.template,'TemplateName','Section');
%         ch1.Title = 'Propeller Design';
%
%         para = Paragraph(strcat('This capter contains the design of propeller named:',' ', char(prop.name)));
%         % append(para,InternalLink('tlarTableRef','refTabella'));
%         add(ch1,para)
%
%         sec1 = Section('TemplateSrc',prop.template,'TemplateName','Section');
%         sec1.Title = 'Top Level Propeller Requirements';
%
%         %updating Propeller structure for reporting
%         prop.TLPR.nblades.value = prop.nblades;
%         prop.TLPR.nblades.Attributes.unit = '-';
%         prop.TLPR.radius.value = prop.radius;
%         prop.TLPR.radius.Attributes.unit = 'meters';
%         prop.TLPR.hub_radius.value = prop.hub_radius;
%         prop.TLPR.hub_radius.Attributes.unit = 'meters';
%         prop.TLPR.whub_radius.value = prop.whub_radius;
%         prop.TLPR.whub_radius.Attributes.unit = 'meters';
%         prop.TLPR.des_speed.value = prop.des_speed;
%         prop.TLPR.des_speed.Attributes.unit = 'm/s';
%         prop.TLPR.adva.value = prop.adva;
%         prop.TLPR.adva.Attributes.unit = '-';
%         prop.TLPR.rpm.value = prop.rpm;
%         prop.TLPR.rpm.Attributes.unit = 'RPM';
%         prop.TLPR.Thrust.value = prop.Thrust;
%         prop.TLPR.Thrust.Attributes.unit = 'Newton';
%         prop.TLPR.Power.value = prop.Power;
%         prop.TLPR.Power.Attributes.unit = 'Watt';
%         prop.TLPR.altitude.value = prop.altitude;
%         prop.TLPR.altitude.Attributes.unit = 'Km';
%         prop.TLPR.cl.value = prop.cl;
%         prop.TLPR.cl.Attributes.unit = '-';
%
%         % prop.J.value = prop.TLPR.des_speed.value/(prop.TLPR.rpm.value/60*2*prop.TLPR.radius.value);
%
%         %% CHAPTER 1 - SECTION 1
%         % TLPR
%         para = Paragraph('The propeller assigned requirements are listed in the following table:');
%         % append(para,InternalLink('tlarTableRef','refTabella'));
%         add(sec1,para)
%
%         tlpr = fieldnames(prop.TLPR);
%         fieldValue = cell(length(tlpr),1);
%         fieldUnit = cell(length(tlpr),1);
%         for i = 1:length(tlpr)
%             fieldValue{i} = prop.TLPR.(tlpr{i}).value;
%             % significant digits
%             if isnumeric(fieldValue{i})
%                 fieldValue{i} = num2str(fieldValue{i},5);
%             end
%             % Not every field has Attributes. If not, they have only one field
%             if length(fieldnames(prop.TLPR.(tlpr{i}))) > 1
%                 fieldUnit{i} = prop.TLPR.(tlpr{i}).Attributes.unit;
%             else % Void cells cannot be converted in strings
%                 fieldUnit{i} = '-';
%             end
%         end
%         tlpr = [tlpr, fieldValue, fieldUnit];
%         header = {'TLPR', 'Value', 'Measure unit'};
%         tbl = FormalTable(header,tlpr);
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
%         tbl.Title = 'Top Level Propeller Requirements';
%         tbl.LinkTarget = 'tlarTableRef';
%
%         add(sec1,tbl);
%
%         % blank line
%         para = Paragraph();
%         append(para, Text(sprintf(' \n ')));
%         add(sec1,para)
%         % end blank line
%
%         if prop.Thrust == 0
%             para2 = Paragraph('The propeller is designed assigned the Power');
%             add(sec1,para2)
%             para3 = Paragraph(['The Power value is' ' ' num2str(prop.TLPR.Power.value) ' ' prop.TLPR.Power.Attributes.unit]);
%             add(sec1,para3)
%         else
%             para2 = Paragraph('The propeller is designed assigned the Thrust');
%             add(sec1,para2)
%             para3 = Paragraph(['The Thrust value is '  ' ' num2str(prop.TLPR.Thrust.value) ' ' prop.TLPR.Thrust.Attributes.unit]);
%             add(sec1,para3)
%         end
%         % % append(para,InternalLink('tlarTableRef','refTabella'));
%
%         % blank line
%         para = Paragraph();
%         append(para, Text(sprintf(' \n ')));
%         add(sec1,para)
%         % end blank line
%
%         % J=V/nD
%         para4 = Paragraph(['The propeller design point will be at J=V/nD equal to ', num2str(prop.J.value)]);
%         add(sec1,para4)
%
%         %propeller MODE of design
%         %prop.mode       = 'POT';propeller design and analysis mode:
%         % GRAD -> Graded Momentum Formulation
%         % POT  -> Potential formulation
%         % VRTX -> Vortex formulation
%
%         if prop.mode == "GRAD"
%             para5 = Paragraph('The propeller design method is "Graded Momentum Formulation"' );
%             add(sec1,para5)
%         else
%             if prop.mode == "POT"
%                 para5 = Paragraph('The propeller design method is "Potential formulation"' );
%                 add(sec1,para5)
%             else
%                 para5 = Paragraph('The propeller design method is "Vortex formulation"' );
%                 add(sec1,para5)
%             end
%         end
%
%         add(ch1,sec1);
%
%         %% CHAPTER 1 - SECTION 2
%         sec2 = Section('TemplateSrc',prop.template,'TemplateName','Section');
%         sec2.Title = 'Blade geometry';
%         para = Paragraph(...
%             ['The main designed blade geometrical characteristics are shown in the following two figures \n',...
%             'The geometry of the propeller is defined according to mininum induced load method. ']);
%         add(sec2,para)
%
%         % Propeller geometrical characteristics c/R and pitch angle
%         fig = FormalImage([results_path,'PropChords.png']);
%         fig.Caption = 'Propeller chords law';
%         fig.Height = '5in';
%         fig.LinkTarget='PropChords';
%         add(sec2,fig);
%         %
%         fig = FormalImage([results_path,'PropPitch.png']);
%         fig.Caption = 'Propeller pitch angle law';
%         fig.Height = '5in';
%         fig.LinkTarget='PropPitch';
%         add(sec2,fig)
%         %
%
%         add(ch1,sec2)
%         % % %add(rpt,ch1) % moved to the end of chapter 2 to include a table of chap2
%         % % %end chapter 1
%
%         %% CHAPTER 2 - PropAnalysis
%         ch2 = Chapter('TemplateSrc',prop.template,'TemplateName','Section');
%         ch2.Title = 'Propeller Analysis';
%         para = Paragraph(strcat('This chapter contains the analysis of propeller named:',' ', char(prop.name)));
%         % append(para,InternalLink('tlarTableRef','refTabella'));
%         add(ch2,para)
%
%         %sec1
%         sec1 = Section('TemplateSrc',prop.template,'TemplateName','Section');
%         sec1.Title = 'Analysis results';
%
%         para = Paragraph(...
%             'The propeller coefficients at design point are shown in the following figures \n');
%         add(sec1,para)
%
%         % Propeller results: CT, CP, eta
%         fig = FormalImage([results_path,'PropCT.png']);
%         fig.Caption = 'Propeller thrust coefficient';
%         fig.Height = '5in';
%         fig.LinkTarget='PropThrust';
%         add(sec1,fig);
%         %
%         fig = FormalImage([results_path,'PropCP.png']);
%         fig.Caption = 'Propeller power coefficient';
%         fig.Height = '5in';
%         fig.LinkTarget='PropPower';
%         add(sec1,fig);
%         %
%         fig = FormalImage([results_path,'PropEta.png']);
%         fig.Caption = 'Propeller efficiency';
%         fig.Height = '5in';
%         fig.LinkTarget='PropEta';
%         add(sec1,fig);
%
%         %table with coefficients
%         results = compose('%.3f',load([results_path, 'Results_', prop.name, '_star.csv']));
%         header = {'J', 'CT', 'CP', 'Eta'};
%         tbl = FormalTable(header,results);
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
%         tbl.Title = 'Propeller analysis results';
%         tbl.LinkTarget = 'PropRes';
%         add(sec1,tbl);
%
%         add(ch2,sec1);
%         %end chapter 2
%         %% CHAPTER - Prop variable pitch analysis
%         ch4 = Chapter('TemplateSrc',prop.template,'TemplateName','Section');
%         ch4.Title = 'Propeller CAD';
%
%         para = Paragraph(strcat('This chapter contains the CAD rendering of designed propeller'));
%         % append(para,InternalLink('tlarTableRef','refTabella'));
%         add(ch4,para)
%
%         sec1 = Section('TemplateSrc',prop.template,'TemplateName','Section');
%         sec1.Title = 'CAD shape';
%
%         % Propeller BLADE shape
%         fig = FormalImage([results_path,'PropBLADE.png']);
%         fig.Caption = 'Propeller BLADE';
%         fig.Height = '5in';
%         fig.LinkTarget='PropBlade';
%         add(sec1,fig);
%
%         % Propeller CAD shape
%         fig = FormalImage([results_path,'PropCAD.png']);
%         fig.Caption = 'Propeller CAD';
%         fig.Height = '5in';
%         fig.LinkTarget='PropCad';
%         add(sec1,fig);
%
%
%         add(ch4,sec1);
%
%         %end chapter 4
%
%
%
%                 %% END
%         %Adding chapters
%         add(rpt,ch1);
%         add(rpt,ch2);
%         add(rpt,ch4);
%
%
%         close(rpt)
%         rptview(rpt)
%         disp('Report successfully written. Execution terminated.')
%
% %     %case {3}
% %         %% CHAPTER 1 - TLPR
% %         ch1 = Chapter('TemplateSrc',prop.template,'TemplateName','Section');
% %         ch1.Title = 'Propeller Design';
% %
% %         para = Paragraph(strcat('This capter contains the design of propeller named:',' ', char(prop.name)));
% %         % append(para,InternalLink('tlarTableRef','refTabella'));
% %         add(ch1,para)
% %
% %         sec1 = Section('TemplateSrc',prop.template,'TemplateName','Section');
% %         sec1.Title = 'Top Level Propeller Requirements';
% %
% %         %updating Propeller structure for reporting
% %         prop.TLPR.nblades.value = prop.nblades;
% %         prop.TLPR.nblades.Attributes.unit = '-';
% %         prop.TLPR.radius.value = prop.radius;
% %         prop.TLPR.radius.Attributes.unit = 'meters';
% %         prop.TLPR.hub_radius.value = prop.hub_radius;
% %         prop.TLPR.hub_radius.Attributes.unit = 'meters';
% %         prop.TLPR.whub_radius.value = prop.whub_radius;
% %         prop.TLPR.whub_radius.Attributes.unit = 'meters';
% %         prop.TLPR.des_speed.value = prop.des_speed;
% %         prop.TLPR.des_speed.Attributes.unit = 'm/s';
% %         prop.TLPR.adva.value = prop.adva;
% %         prop.TLPR.adva.Attributes.unit = '-';
% %         prop.TLPR.rpm.value = prop.rpm;
% %         prop.TLPR.rpm.Attributes.unit = 'RPM';
% %         prop.TLPR.Thrust.value = prop.Thrust;
% %         prop.TLPR.Thrust.Attributes.unit = 'Newton';
% %         prop.TLPR.Power.value = prop.Power;
% %         prop.TLPR.Power.Attributes.unit = 'Watt';
% %         prop.TLPR.altitude.value = prop.altitude;
% %         prop.TLPR.altitude.Attributes.unit = 'Km';
% %         prop.TLPR.cl.value = prop.cl;
% %         prop.TLPR.cl.Attributes.unit = '-';
% %
% %         % prop.J.value = prop.TLPR.des_speed.value/(prop.TLPR.rpm.value/60*2*prop.TLPR.radius.value);
% %
% %         %% CHAPTER 1 - SECTION 1
% %         % TLPR
% %         para = Paragraph('The propeller assigned requirements are listed in the following table:');
% %         % append(para,InternalLink('tlarTableRef','refTabella'));
% %         add(sec1,para)
% %
% %         tlpr = fieldnames(prop.TLPR);
% %         fieldValue = cell(length(tlpr),1);
% %         fieldUnit = cell(length(tlpr),1);
% %         for i = 1:length(tlpr)
% %             fieldValue{i} = prop.TLPR.(tlpr{i}).value;
% %             % significant digits
% %             if isnumeric(fieldValue{i})
% %                 fieldValue{i} = num2str(fieldValue{i},5);
% %             end
% %             % Not every field has Attributes. If not, they have only one field
% %             if length(fieldnames(prop.TLPR.(tlpr{i}))) > 1
% %                 fieldUnit{i} = prop.TLPR.(tlpr{i}).Attributes.unit;
% %             else % Void cells cannot be converted in strings
% %                 fieldUnit{i} = '-';
% %             end
% %         end
% %         tlpr = [tlpr, fieldValue, fieldUnit];
% %         header = {'TLPR', 'Value', 'Measure unit'};
% %         tbl = FormalTable(header,tlpr);
% %         tbl.Style = {...
% %             RowSep('solid','black','2px'),...
% %             ColSep('solid','black','1px'),...
% %             Border('ridge','black','5px')...
% %             };
% %         tbl.Header.Style = {...
% %             RowSep('solid','black','2px'),...
% %             ColSep('solid','black','2px'),...
% %             Border('ridge','black','5px')...
% %             };
% %         tbl.TableEntriesStyle = {HAlign('left')};
% %         % In order to put a table with a caption, the API Report denomination should
% %         % be used, the other options are from API DOM. In order to solve the problem,
% %         % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
% %         tbl = BaseTable(tbl);
% %         tbl.Title = 'Top Level Propeller Requirements';
% %         tbl.LinkTarget = 'tlarTableRef';
% %
% %         add(sec1,tbl);
% %
% %         % blank line
% %         para = Paragraph();
% %         append(para, Text(sprintf(' \n ')));
% %         add(sec1,para)
% %         % end blank line
% %
% %         if prop.Thrust == 0
% %             para2 = Paragraph('The propeller is designed assigned the Power');
% %             add(sec1,para2)
% %             para3 = Paragraph(['The Power value is' ' ' num2str(prop.TLPR.Power.value) ' ' prop.TLPR.Power.Attributes.unit]);
% %             add(sec1,para3)
% %         else
% %             para2 = Paragraph('The propeller is designed assigned the Thrust');
% %             add(sec1,para2)
% %             para3 = Paragraph(['The Thrust value is '  ' ' num2str(prop.TLPR.Thrust.value) ' ' prop.TLPR.Thrust.Attributes.unit]);
% %             add(sec1,para3)
% %         end
% %         % % append(para,InternalLink('tlarTableRef','refTabella'));
% %
% %         % blank line
% %         para = Paragraph();
% %         append(para, Text(sprintf(' \n ')));
% %         add(sec1,para)
% %         % end blank line
% %
% %         % J=V/nD
% %         para4 = Paragraph(['The propeller design point will be at J=V/nD equal to ', num2str(prop.J.value)]);
% %         add(sec1,para4)
% %
% %         %propeller MODE of design
% %         %prop.mode       = 'POT';propeller design and analysis mode:
% %         % GRAD -> Graded Momentum Formulation
% %         % POT  -> Potential formulation
% %         % VRTX -> Vortex formulation
% %
% %         if prop.mode == "GRAD"
% %             para5 = Paragraph('The propeller design method is "Graded Momentum Formulation"' );
% %             add(sec1,para5)
% %         else
% %             if prop.mode == "POT"
% %                 para5 = Paragraph('The propeller design method is "Potential formulation"' );
% %                 add(sec1,para5)
% %             else
% %                 para5 = Paragraph('The propeller design method is "Vortex formulation"' );
% %                 add(sec1,para5)
% %             end
% %         end
% %
% %         add(ch1,sec1);
% %
% %         %% CHAPTER 1 - SECTION 2
% %         sec2 = Section('TemplateSrc',prop.template,'TemplateName','Section');
% %         sec2.Title = 'Blade geometry';
% %         para = Paragraph(...
% %             ['The main designed blade geometrical characteristics are shown in the following two figures \n',...
% %             'The geometry of the propeller is defined according to mininum induced load method. ']);
% %         add(sec2,para)
% %
% %         % Propeller geometrical characteristics c/R and pitch angle
% %         fig = FormalImage([results_path,'PropChords.png']);
% %         fig.Caption = 'Propeller chords law';
% %         fig.Height = '5in';
% %         fig.LinkTarget='PropChords';
% %         add(sec2,fig);
% %         %
% %         fig = FormalImage([results_path,'PropPitch.png']);
% %         fig.Caption = 'Propeller pitch angle law';
% %         fig.Height = '5in';
% %         fig.LinkTarget='PropPitch';
% %         add(sec2,fig)
% %         %
% %
% %         add(ch1,sec2)
% %         % % %add(rpt,ch1) % moved to the end of chapter 2 to include a table of chap2
% %         % % %end chapter 1
% %
% %         %% CHAPTER 2 - PropAnalysis
% %         ch2 = Chapter('TemplateSrc',prop.template,'TemplateName','Section');
% %         ch2.Title = 'Propeller Analysis';
% %         para = Paragraph(strcat('This chapter contains the analysis of propeller named:',' ', char(prop.name)));
% %         % append(para,InternalLink('tlarTableRef','refTabella'));
% %         add(ch2,para)
% %
% %         %sec1
% %         sec1 = Section('TemplateSrc',prop.template,'TemplateName','Section');
% %         sec1.Title = 'Analysis results';
% %
% %         para = Paragraph(...
% %             'The propeller coefficients at design point are shown in the following figures \n');
% %         add(sec1,para)
% %
% %         % Propeller results: CT, CP, eta
% %         fig = FormalImage([results_path,'PropCT.png']);
% %         fig.Caption = 'Propeller thrust coefficient';
% %         fig.Height = '5in';
% %         fig.LinkTarget='PropThrust';
% %         add(sec1,fig);
% %         %
% %         fig = FormalImage([results_path,'PropCP.png']);
% %         fig.Caption = 'Propeller power coefficient';
% %         fig.Height = '5in';
% %         fig.LinkTarget='PropPower';
% %         add(sec1,fig);
% %         %
% %         fig = FormalImage([results_path,'PropEta.png']);
% %         fig.Caption = 'Propeller efficiency';
% %         fig.Height = '5in';
% %         fig.LinkTarget='PropEta';
% %         add(sec1,fig);
% %
% %         %table with coefficients
% %         results = compose('%.3f',load([results_path, 'Results_', prop.name, '_star.csv']));
% %         header = {'J', 'CT', 'CP', 'Eta'};
% %         tbl = FormalTable(header,results);
% %         tbl.Style = {...
% %             RowSep('solid','black','2px'),...
% %             ColSep('solid','black','1px'),...
% %             Border('ridge','black','5px')...
% %             };
% %         tbl.Header.Style = {...
% %             RowSep('solid','black','2px'),...
% %             ColSep('solid','black','2px'),...
% %             Border('ridge','black','5px')...
% %             };
% %         tbl.TableEntriesStyle = {HAlign('left')};
% %         % In order to put a table with a caption, the API Report denomination should
% %         % be used, the other options are from API DOM. In order to solve the problem,
% %         % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
% %         tbl = BaseTable(tbl);
% %         tbl.Title = 'Propeller analysis results';
% %         tbl.LinkTarget = 'PropRes';
% %         add(sec1,tbl);
% %
% %         add(ch2,sec1);
% %         %end chapter 2
% %
% %         %% CHAPTER 3 - Prop variable pitch analysis
% %         ch3 = Chapter('TemplateSrc',prop.template,'TemplateName','Section');
% %         ch3.Title = 'Propeller Variable Pitch';
% %
% %         para = Paragraph('This chapter contains the analysis of designed propeller with pitch angle variation');
% %         % append(para,InternalLink('tlarTableRef','refTabella'));
% %         add(ch3,para)
% %
% %         %sec1
% %         sec1 = Section('TemplateSrc',prop.template,'TemplateName','Section');
% %         sec1.Title = 'Analysis variable pitch results';
% %
% %         para = Paragraph(...
% %             'The propeller coefficients at variable pitch angles are shown in the following figures \n');
% %         add(sec1,para)
% %
% %         % Propeller results: CT, CP, eta
% %         fig = FormalImage([results_path,'PropCT_vp.png']);
% %         fig.Caption = 'Propeller thrust coefficient variable pitch';
% %         fig.Height = '5in';
% %         fig.LinkTarget='PropThrustvp';
% %         add(sec1,fig);
% %         %
% %         fig = FormalImage([results_path,'PropCP_vp.png']);
% %         fig.Caption = 'Propeller power coefficient variable pitch';
% %         fig.Height = '5in';
% %         fig.LinkTarget='PropPowervp';
% %         add(sec1,fig);
% %         %
% %         fig = FormalImage([results_path,'PropEta_vp.png']);
% %         fig.Caption = 'Propeller efficiency variable pitch';
% %         fig.Height = '5in';
% %         fig.LinkTarget='PropEtavp';
% %         add(sec1,fig);
% %
% %         %table with coefficients
% %         % results = compose('%.3f',load([results_path, 'Results_', prop.name, '_star.csv']));
% %         % header = {'J', 'CT', 'CP', 'Eta'};
% %         % tbl = FormalTable(header,results);
% %         % tbl.Style = {...
% %         % RowSep('solid','black','2px'),...
% %         % ColSep('solid','black','1px'),...
% %         % Border('ridge','black','5px')...
% %         % };
% %         % tbl.Header.Style = {...
% %         % RowSep('solid','black','2px'),...
% %         % ColSep('solid','black','2px'),...
% %         % Border('ridge','black','5px')...
% %         % };
% %         % tbl.TableEntriesStyle = {HAlign('left')};
% %         % % In order to put a table with a caption, the API Report denomination should
% %         % % be used, the other options are from API DOM. In order to solve the problem,
% %         % % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
% %         % tbl = BaseTable(tbl);
% %         % tbl.Title = 'Propeller analysis results';
% %         % tbl.LinkTarget = 'PropRes';
% %         % add(sec1,tbl);
% %
% %         add(ch3,sec1);
% %
% %
% %
% %
% %         %end chapter 3
% %
% %         %% CHAPTER 4 - Prop variable pitch analysis
% %         ch4 = Chapter('TemplateSrc',prop.template,'TemplateName','Section');
% %         ch4.Title = 'Propeller CAD';
% %
% %         para = Paragraph(strcat('This chapter contains the CAD rendering of designed propeller'));
% %         % append(para,InternalLink('tlarTableRef','refTabella'));
% %         add(ch4,para)
% %
% %         sec1 = Section('TemplateSrc',prop.template,'TemplateName','Section');
% %         sec1.Title = 'CAD shape';
% %
% %         % Propeller BLADE shape
% %         fig = FormalImage([results_path,'PropBLADE.png']);
% %         fig.Caption = 'Propeller BLADE';
% %         fig.Height = '5in';
% %         fig.LinkTarget='PropBlade';
% %         add(sec1,fig);
% %
% %         % Propeller CAD shape
% %         fig = FormalImage([results_path,'PropCAD.png']);
% %         fig.Caption = 'Propeller CAD';
% %         fig.Height = '5in';
% %         fig.LinkTarget='PropCad';
% %         add(sec1,fig);
% %
% %
% %         add(ch4,sec1);
% %
% %         %end chapter 4
% %
% %
% %         %% END
% %         %Adding chapters
% %         add(rpt,ch1);
% %         add(rpt,ch2);
% %         add(rpt,ch3);
% %         add(rpt,ch4);
% %
% %
% %         close(rpt)
% %         rptview(rpt)
% %         disp('Report successfully written. Execution terminated.')
% %
% end
%
close(rpt)
rptview(rpt)
disp('Report successfully written. Execution terminated.')
end
