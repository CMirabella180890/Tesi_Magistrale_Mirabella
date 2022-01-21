%% CHAPTER 4 - Aircraft data
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

ch = Chapter();
ch.Title = 'Aircraft data';

str = ['The aircraft geometrical, masses, inertial and aerodynamic data, useful for flight loads estimation are summarized in this chapter.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
% append(para,InternalLink('tlarTableRef','refTabella'));
add(ch,para)

%         %% CHAPTER 4 - SECTION 1
%         % geometry
sec = Section();
sec.Title = 'Geometry';
disp(['Chapter 4', (' "'), ch.Title,('" ') ,'writing...' ])

para = Paragraph('The aircraft reference geometrical characteristics are summarized in the following tables.');

% wing
if isfield(Aircraft.Geometry, 'Wing')==1
append(para,InternalLink('wingTableRef','Wing parameters'));
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
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).      
        tbl = BaseTable(tbl);
        tbl.Title = 'Wing parameters';
        tbl.LinkTarget = 'wingTableRef';     
        add(sec,tbl);
end
%
%horizontal
if isfield(Aircraft.Geometry, 'Horizontal')==1
%append(para,InternalLink('horiTableRef','Horizontal Tail parameters'));
%add(sec,para)
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
         header = {'Horizontal parameters', 'Value', 'Measure unit'};
         tbl = FormalTable(header,horizontal);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).      
        tbl = BaseTable(tbl);
        tbl.Title = 'Horizontal Tail parameters';
        tbl.LinkTarget = 'horizontalTableRef';     
        add(sec,tbl);
end

%vertical
if isfield(Aircraft.Geometry, 'Vertical')==1
% append(para,InternalLink('verTableRef','Vertical Tail parameters'));
% add(sec,para)
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
         header = {'Vertical parameters', 'Value', 'Measure unit'};
         tbl = FormalTable(header,vertical);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).      
        tbl = BaseTable(tbl);
        tbl.Title = 'Vertical Tail parameters';
        tbl.LinkTarget = 'verticalTableRef';     
        add(sec,tbl);
end

%fuselage
if isfield(Aircraft.Geometry, 'Fuselage')==1
    %remove fields
     Aircraft.Geometry.Fuselage = rmfield (Aircraft.Geometry.Fuselage, 'id');
     Aircraft.Geometry.Fuselage = rmfield (Aircraft.Geometry.Fuselage, 'type');
     Aircraft.Geometry.Fuselage = rmfield (Aircraft.Geometry.Fuselage, 'empennage');
     
% append(para,InternalLink('verTableRef','Fuselage parameters'));
% add(sec,para)
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
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).      
        tbl = BaseTable(tbl);
        tbl.Title = 'Fuselage parameters';
        tbl.LinkTarget = 'fuselageTableRef';     
        add(sec,tbl);
end

%landing gear
if isfield(Aircraft.Geometry, 'LandGear')==1
% append(para,InternalLink('lgTableRef','Landing Gear parameters'));
% add(sec,para)
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
         header = {'Landing gear parameters', 'Value', 'Measure unit'};
         tbl = FormalTable(header,landGear);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).      
        tbl = BaseTable(tbl);
        tbl.Title = 'Landing Gear parameters';
        tbl.LinkTarget = 'landGearTableRef';     
        add(sec,tbl);
end

%Elevator
if isfield(Aircraft.Geometry, 'Elevator')==1
% append(para,InternalLink('elTableRef','Elevator parameters'));
% add(sec,para)
        elevator = fieldnames(Aircraft.Geometry.Elevator);
         fieldValue = cell(length(elevator),1);
         fieldUnit = cell(length(elevator),1);
         for i = 1:length(elevator)
             fieldValue{i} = Aircraft.Geometry.Elevator.(elevator{i}).value;
             % significant digits
             if isnumeric(fieldValue{i})
                 fieldValue{i} = num2str(fieldValue{i},5);
             end
             % Not every field has Attributes. If not, they have only one field
             if length(fieldnames(Aircraft.Geometry.Elevator.(elevator{i}))) > 1
                 fieldUnit{i} = Aircraft.Geometry.Elevator.(elevator{i}).Attributes.unit;
             else % Void cells cannot be converted in strings
                 fieldUnit{i} = '-';
             end
         end
         elevator = [elevator, fieldValue, fieldUnit];
         header = {'Elevator parameters', 'Value', 'Measure unit'};
         tbl = FormalTable(header,elevator);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).      
        tbl = BaseTable(tbl);
        tbl.Title = 'Elevator parameters';
        tbl.LinkTarget = 'elTableRef';     
        add(sec,tbl);
end

%Rudder
if isfield(Aircraft.Geometry, 'Rudder')==1
% append(para,InternalLink('rudTableRef','Rudder parameters'));
% add(sec,para)
        rudder = fieldnames(Aircraft.Geometry.Rudder);
         fieldValue = cell(length(rudder),1);
         fieldUnit = cell(length(rudder),1);
         for i = 1:length(rudder)
             fieldValue{i} = Aircraft.Geometry.Rudder.(rudder{i}).value;
             % significant digits
             if isnumeric(fieldValue{i})
                 fieldValue{i} = num2str(fieldValue{i},5);
             end
             % Not every field has Attributes. If not, they have only one field
             if length(fieldnames(Aircraft.Geometry.Rudder.(rudder{i}))) > 1
                 fieldUnit{i} = Aircraft.Geometry.Rudder.(rudder{i}).Attributes.unit;
             else % Void cells cannot be converted in strings
                 fieldUnit{i} = '-';
             end
         end
         rudder = [rudder, fieldValue, fieldUnit];
         header = {'Rudder parameters', 'Value', 'Measure unit'};
         tbl = FormalTable(header,rudder);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).      
        tbl = BaseTable(tbl);
        tbl.Title = 'Rudder parameters';
        tbl.LinkTarget = 'rudTableRef';     
        add(sec,tbl);
end


%Aileron
if isfield(Aircraft.Geometry, 'Aileron')==1
% append(para,InternalLink('ailTableRef','Aileron parameters'));
% add(sec,para)
        aileron = fieldnames(Aircraft.Geometry.Aileron);
         fieldValue = cell(length(aileron),1);
         fieldUnit = cell(length(aileron),1);
         for i = 1:length(aileron)
             fieldValue{i} = Aircraft.Geometry.Aileron.(aileron{i}).value;
             % significant digits
             if isnumeric(fieldValue{i})
                 fieldValue{i} = num2str(fieldValue{i},5);
             end
             % Not every field has Attributes. If not, they have only one field
             if length(fieldnames(Aircraft.Geometry.Aileron.(aileron{i}))) > 1
                 fieldUnit{i} = Aircraft.Geometry.Aileron.(aileron{i}).Attributes.unit;
             else % Void cells cannot be converted in strings
                 fieldUnit{i} = '-';
             end
         end
         aileron = [aileron, fieldValue, fieldUnit];
         header = {'Aileron parameters', 'Value', 'Measure unit'};
         tbl = FormalTable(header,aileron);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).      
        tbl = BaseTable(tbl);
        tbl.Title = 'Aileron parameters';
        tbl.LinkTarget = 'ailTableRef';     
        add(sec,tbl);
end

%sec
%masses and inertia
add(ch,sec);

sec = Section();
sec.Title = 'Masses and inertia';

para = Paragraph('The aircraft reference masses and inertia are summarized in this subsection');
add(sec,para)

%Aileron
if isfield(Aircraft.Weight, 'I_Level')==1
 
 para = Paragraph ('The Aircraft masses and inertia are summarized in Table: ');
 append(para,InternalLink('weiTableRef','Weight parameters'));
 add(sec,para)
        weight = fieldnames(Aircraft.Weight.I_Level);
         fieldValue = cell(length(weight),1);
         fieldUnit = cell(length(weight),1);
         for i = 1:length(weight)
             fieldValue{i} = Aircraft.Weight.I_Level.(weight{i}).value;
             % significant digits
             if isnumeric(fieldValue{i})
                 fieldValue{i} = num2str(fieldValue{i},5);
             end
             % Not every field has Attributes. If not, they have only one field
             if length(fieldnames(Aircraft.Weight.I_Level.(weight{i}))) > 1
                 fieldUnit{i} = Aircraft.Weight.I_Level.(weight{i}).Attributes.unit;
             else % Void cells cannot be converted in strings
                 fieldUnit{i} = '-';
             end
         end
         weight = [weight, fieldValue, fieldUnit];
         header = {'Weight', 'Value', 'Measure unit'};
         tbl = FormalTable(header,weight);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).      
        tbl = BaseTable(tbl);
        tbl.Title = 'Weight parameters';
        tbl.LinkTarget = 'weiTableRef';     
        add(sec,tbl);
end



add(ch,sec);



%         %% CHAPTER 4 - SECTION 3
%         % aero
sec = Section();
sec.Title = 'Aerodynamic';

para = Paragraph('The aircraft reference aerodynamic is in figure: ');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec,para)
append(para,InternalLink('polars_wb',' Wing-Body reference Aerodynamics'));

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