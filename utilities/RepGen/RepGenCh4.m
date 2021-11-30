%% CHAPTER 4 - Aircraft data
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

ch = Chapter();
ch.Title = 'Aircraft data';

str = ['Add here all the aircraft geometrical, aero and inertial and masses data useful for following paragraph'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
% append(para,InternalLink('tlarTableRef','refTabella'));
add(ch,para)

%         %% CHAPTER 4 - SECTION 1
%         % geometry
sec = Section();
sec.Title = 'Geometry';

para = Paragraph('The aircraft reference geometry is summarized in table:');

% wing
if isfield(Aircraft.Geometry, 'Wing')==1
append(para,InternalLink('wingTableRef','Ref:wing'));
add(sec,para)
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
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).      
        tbl = BaseTable(tbl);
        tbl.Title = 'Wing Geometrical Parameters';
        tbl.LinkTarget = 'wingTableRef';     
        add(sec,tbl);
end
%
%horizontal
if isfield(Aircraft.Geometry, 'Horizontal')==1
append(para,InternalLink('horiTableRef','Ref:hori'));
add(sec,para)
         horizontal = fieldnames(Aircraft.Geometry.Horizontal);
         fieldValue = cell(length(horizontal),1);
         fieldUnit = cell(length(horizontal),1);
         for i = 1:length(horizontal)
             fieldValue{i} = Aircraft.Geometry.Horizontal.(horizontal{i}).value;
             % significant digits
             if isnumeric(fieldValue{i})
                 fieldValue{i} = num2str(fieldValue{i},5);
             end
             % Not every field has Attributes. If not, they have only one field
             if length(fieldnames(Aircraft.Geometry.Horizontal.(horizontal{i}))) > 1
                 fieldUnit{i} = Aircraft.Geometry.Horizontal.(horizontal{i}).Attributes.unit;
             else % Void cells cannot be converted in strings
                 fieldUnit{i} = '-';
             end
         end
         horizontal = [horizontal, fieldValue, fieldUnit];
         header = {'horizontal parameters', 'Value', 'Measure unit'};
         tbl = FormalTable(header,horizontal);
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
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).      
        tbl = BaseTable(tbl);
        tbl.Title = 'horizontal Geometrical Parameters';
        tbl.LinkTarget = 'horizontalTableRef';     
        add(sec,tbl);
end

%vertical
if isfield(Aircraft.Geometry, 'Vertical')==1
append(para,InternalLink('verTableRef','Ref:vert'));
add(sec,para)
         vertical = fieldnames(Aircraft.Geometry.Vertical);
         fieldValue = cell(length(vertical),1);
         fieldUnit = cell(length(vertical),1);
         for i = 1:length(vertical)
             fieldValue{i} = Aircraft.Geometry.Vertical.(vertical{i}).value;
             % significant digits
             if isnumeric(fieldValue{i})
                 fieldValue{i} = num2str(fieldValue{i},5);
             end
             % Not every field has Attributes. If not, they have only one field
             if length(fieldnames(Aircraft.Geometry.Vertical.(vertical{i}))) > 1
                 fieldUnit{i} = Aircraft.Geometry.Vertical.(vertical{i}).Attributes.unit;
             else % Void cells cannot be converted in strings
                 fieldUnit{i} = '-';
             end
         end
         vertical = [vertical, fieldValue, fieldUnit];
         header = {'vertical parameters', 'Value', 'Measure unit'};
         tbl = FormalTable(header,vertical);
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
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).      
        tbl = BaseTable(tbl);
        tbl.Title = 'Vertical Geometrical Parameters';
        tbl.LinkTarget = 'verticalTableRef';     
        add(sec,tbl);
end

%fuselage
if isfield(Aircraft.Geometry, 'Fuselage')==1
    %remove fields
     Aircraft.Geometry.Fuselage = rmfield (Aircraft.Geometry.Fuselage, 'id');
     Aircraft.Geometry.Fuselage = rmfield (Aircraft.Geometry.Fuselage, 'type');
         
append(para,InternalLink('verTableRef','Ref:fus'));
add(sec,para)
         fuselage = fieldnames(Aircraft.Geometry.Fuselage);
         fieldValue = cell(length(fuselage),1);
         fieldUnit = cell(length(fuselage),1);
         for i = 1:length(fuselage)

             fieldValue{i} = Aircraft.Geometry.Fuselage.(fuselage{i}).value;
             % significant digits
             if isnumeric(fieldValue{i})
                 fieldValue{i} = num2str(fieldValue{i},5);
             end
             % Not every field has Attributes. If not, they have only one field
             if length(fieldnames(Aircraft.Geometry.Fuselage.(fuselage{i}))) > 1
                 fieldUnit{i} = Aircraft.Geometry.Fuselage.(fuselage{i}).Attributes.unit;
             else % Void cells cannot be converted in strings
                 fieldUnit{i} = '-';
             end
         end
         fuselage = [fuselage, fieldValue, fieldUnit];
         header = {'Fuselage parameters', 'Value', 'Measure unit'};
         tbl = FormalTable(header,fuselage);
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
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).      
        tbl = BaseTable(tbl);
        tbl.Title = 'Fuselage Geometrical Parameters';
        tbl.LinkTarget = 'fuselageTableRef';     
        add(sec,tbl);
end

%landing gear
if isfield(Aircraft.Geometry, 'LandGear')==1
append(para,InternalLink('lgTableRef','Ref:lg'));
add(sec,para)
        landGear = fieldnames(Aircraft.Geometry.LandGear);
         fieldValue = cell(length(landGear),1);
         fieldUnit = cell(length(fuselage),1);
         for i = 1:length(landGear)
             fieldValue{i} = Aircraft.Geometry.LandGear.(landGear{i}).value;
             % significant digits
             if isnumeric(fieldValue{i})
                 fieldValue{i} = num2str(fieldValue{i},5);
             end
             % Not every field has Attributes. If not, they have only one field
             if length(fieldnames(Aircraft.Geometry.LandGear.(landGear{i}))) > 1
                 fieldUnit{i} = Aircraft.Geometry.landGear.(landGear{i}).Attributes.unit;
             else % Void cells cannot be converted in strings
                 fieldUnit{i} = '-';
             end
         end
         header = {'Fuselage parameters', 'Value', 'Measure unit'};
         tbl = FormalTable(header,landGear);
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
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).      
        tbl = BaseTable(tbl);
        tbl.Title = 'Landing Gear Geometrical Parameters';
        tbl.LinkTarget = 'landGearTableRef';     
        add(sec,tbl);
end

add(ch,sec);

%         %% CHAPTER 4 - SECTION 2
%         % aero
sec = Section();
sec.Title = 'Aerodynamic';

para = Paragraph('The aircraft reference aerodynamic is in figure:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec,para)

append(para,InternalLink('lgTableRef','polars_wb'));

%moving to another path for figure
cd ..
cd ..
 regulation = Aircraft.Certification.Regulation.value;
 results_path = [pwd '\' regulation '\Output\'];

cd (RepDir);

fig = FormalImage([results_path,'Polars.png']);
         fig.Caption = 'Wing-Body reference Aerodynamics';
         fig.Height = '5in';
         fig.LinkTarget='polars_wb';
         add(sec,fig);

add(ch,sec);


%% END chapter
%Adding chapters
add(rpt,ch);