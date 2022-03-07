%% CHAPTER 4 - Aircraft data
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

% CAP 4:
% 		
% aggiungere le figure di aerodinamica CL, CD e CMWB
% 		
% aggiungere i principali coefficienti 3D velivolo completo in tabella (CLMAX, CLMAXTO, CLMAX LAN, CLalha)

%% FIGURES TO ADD 

% RepGenFigCh4 = figure(47);

alpha = Aircraft.Certification.Aerodynamic_data.alpha.value;
CL    = Aircraft.Certification.Aerodynamic_data.CL.value;
CD    = Aircraft.Certification.Aerodynamic_data.CD.value;
CM    = Aircraft.Certification.Aerodynamic_data.CM.value;
% -------------------------------------------------------------------------
RepGenFigCh4 = figure(47);
subplot(2,2,1);
plot(str2num(alpha), str2num(CL));
grid on; grid minor;
xlabel("Angle of attack - $\alpha$ (deg)", "Interpreter", "latex")
ylabel("Lift coefficient - $C_L$ ", "Interpreter", "latex")
title('$C_L$ vs $\alpha$', "Interpreter", "latex")

subplot(2,2,2);
plot(str2num(CL), str2num(CD));
grid on; grid minor;
xlabel("Lift coefficient - $C_L$ ", "Interpreter", "latex")
ylabel("Drag coefficient - $C_D$ ", "Interpreter", "latex")
title('$C_L$ vs $C_D$', "Interpreter", "latex")

subplot(2,2,[3,4]);
plot(str2num(alpha), str2num(CM));
grid on; grid minor;
xlabel("Angle of attack - $\alpha$ (deg)", "Interpreter", "latex")
ylabel("Pitching moment coefficient - $C_M$ ", "Interpreter", "latex")
title('$C_M$ vs $\alpha$', "Interpreter", "latex")

% -------------------------------------------------------------------------
% Saving figures inside correct folder
cd .. 
cd ..
cd csvla
dir = pwd;
fprintf("--------------------------------------");
fprintf('\n');
fprintf('### Saving outpus inside correct Folder ###');
fprintf('\n');
SaveFolder = strcat(dir,'\Output');
fprintf('Saving RepGenFigCh4.png in: ');
fprintf('\n');      
fprintf('%s\n', SaveFolder);
% EXPORT FIGURE
% exportgraphics(final_envelope, 'Finalenvelope.pdf', 'ContentType', 'vector')
exportgraphics(RepGenFigCh4, 'RepGenFigCh4.png', 'ContentType', 'vector')

% Moving file inside correct folder
% movefile Finalenvelope.pdf Output
movefile RepGenFigCh4.png Output 

cd ..
cd utilities/RepGen

close(RepGenFigCh4);

%% REPORT GEN CHAPTER 4

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
%      Aircraft.Geometry.Fuselage = rmfield (Aircraft.Geometry.Fuselage, 'empennage');
     
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


%Flaps
if isfield(Aircraft.Geometry, 'Flaps')==1
% append(para,InternalLink('elTableRef','Elevator parameters'));
% add(sec,para)
        flaps = fieldnames(Aircraft.Geometry.Flaps);
         fieldValue = cell(length(flaps),1);
         fieldUnit = cell(length(flaps),1);
         for i = 1:length(flaps)
             fieldValue{i} = Aircraft.Geometry.Flaps.(flaps{i}).value;
             % significant digits
             if isnumeric(fieldValue{i})
                 fieldValue{i} = num2str(fieldValue{i},5);
             end
             % Not every field has Attributes. If not, they have only one field
             if length(fieldnames(Aircraft.Geometry.Flaps.(flaps{i}))) > 1
                 fieldUnit{i} = Aircraft.Geometry.Flaps.(flaps{i}).Attributes.unit;
             else % Void cells cannot be converted in strings
                 fieldUnit{i} = '-';
             end
         end
         flaps = [flaps, fieldValue, fieldUnit];
         header = {'Flaps parameters', 'Value', 'Measure unit'};
         tbl = FormalTable(header,flaps);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).      
        tbl = BaseTable(tbl);
        tbl.Title = 'Flaps parameters';
        tbl.LinkTarget = 'flTableRef';     
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
% -------------------------------------------------------------------------
        % -----------------------------------------------------------------
        % moving to another path for figure
        cd ..
        cd ..
        regulation = Aircraft.Certification.Regulation.value;
        results_path = [pwd '\' regulation '\Output\'];
        
        cd (RepDir);
        
        fig = FormalImage([results_path,'FullLiftModelInterpolation.png']);
        fig.Caption = 'Lift coefficient of the 3D wing-body configuration.';
        fig.Height = '4in';
        fig.Width = '4in';
        fig.LinkTarget='3DwingbodyaerodynamicLift';
        add(sec,fig);
        fig = FormalImage([results_path,'FullDragModelInterpolation.png']);
        fig.Caption = 'Drag coefficient of the 3D wing-body configuration.';
        fig.Height = '4in';
        fig.Width = '4in';
        fig.LinkTarget='3DwingbodyaerodynamicDrag';
        add(sec,fig);
        fig = FormalImage([results_path,'FullPitchMomModelInterpolation.png']);
        fig.Caption = 'Pitch moment coefficient of the 3D wing-body configuration.';
        fig.Height = '4in';
        fig.Width = '4in';
        fig.LinkTarget='3DwingbodyaerodynamicPitchMom';
        add(sec,fig);
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
% -------------------------------------------------------------------------
para = Paragraph('The aircraft reference aerodynamic is shown in figure: ');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec,para)
append(para,InternalLink('polars_wb',' Wing-Body reference Aerodynamics'));
% -------------------------------------------------------------------------
        str = ['The previous figures refer to wing-body aerodynamics,' ...
            ' obtained via computational fluidynamics simulations.' ...
            ' Full-vehicle aerodynamics is not available for this aircraft.' ...
            ' The calculations illustrated in this document are based on an' ...
            ' aerodynamic reference model, obtained from wing-body aerodynamics' ...
            ' through interpolation (red curves in figure). The interpolation' ...
            ' procedure is in good agreement with present data. Further' ...
            ' examination of this data through wind-tunnel' ...
            ' or in-flight tests must be performed.'];
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(sec,para);  
        % -----------------------------------------------------------------
%moving to another path for figure
cd ..
cd ..
 regulation = Aircraft.Certification.Regulation.value;
 results_path = [pwd '\' regulation '\Output\'];

cd (RepDir);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
Max_Lift_Coefficient          = Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value;
Max_Lift_Coefficient_unit     = Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.Attributes.unit;
Max_Lift_Coefficient_inverted = Aircraft.Certification.Aerodynamic_data.Min_Lift_Coefficient.value;
CLMAX_takeoff                 = Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_takeoff.value;
CLMAX_takeoff_unit            = Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_takeoff.Attributes.unit;
CLMAX_landing                 = Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_landing.value;
CLMAX_landing_unit            = Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_landing.Attributes.unit;
CL0                           = Aircraft.Certification.Aerodynamic_data.CL0.value;
CL0_unit                      = Aircraft.Certification.Aerodynamic_data.CL0.Attributes.unit;
CLALFA_rad                    = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
CLALFA_rad_unit               = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.Attributes.unit;
CLALFA_deg                    = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;
CLALFA_deg_unit               = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.Attributes.unit;
CD0                           = Aircraft.Certification.Aerodynamic_data.CD0.value;
CD0_unit                      = Aircraft.Certification.Aerodynamic_data.CD0.Attributes.unit;
CDGEAR                        = Aircraft.Certification.Aerodynamic_data.CD_landing_gear.value;
CDGEAR_UNIT                   = Aircraft.Certification.Aerodynamic_data.CD_landing_gear.Attributes.unit;
CM0                           = Aircraft.Certification.Aerodynamic_data.CM0.value;
CM0_unit                      = Aircraft.Certification.Aerodynamic_data.CM0.Attributes.unit;
CMCL                          = Aircraft.Certification.Aerodynamic_data.CMCL.value;
CMCL_unit                     = Aircraft.Certification.Aerodynamic_data.CMCL.Attributes.unit;
CM_landing_gear               = Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value;
CM_landing_gear_unit          = Aircraft.Certification.Aerodynamic_data.CM_landing_gear.Attributes.unit;
Alpha_zero_lift               = Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
Alpha_zero_lift_unit          = Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.Attributes.unit;
Oswald_efficiency_factor      = Aircraft.Certification.Aerodynamic_data.e.value;
Oswald_efficiency_factor_un   = Aircraft.Certification.Aerodynamic_data.e.Attributes.unit;
% -------------------------------------------------------------------------
        %table gust calculation        
        str = ['TABLE TO BE CHECKED!!!'];
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        para.BackgroundColor = "red";
        add(ch,para);  
        
        header = {'Quantity', 'Value', 'Units'};
        %each table row needs of a fieldValue
        %gust calculation table
        %1 
        cl_wb_max_clean = {'CLMAX',              num2str(Max_Lift_Coefficient,4),          char(Max_Lift_Coefficient_unit)};    
        cl_wb_inv_clean = {'CLMAX inverted',     num2str(Max_Lift_Coefficient_inverted,4), char(Max_Lift_Coefficient_unit)};
        cl_wb_takeoff   = {'CLMAX takeoff',      num2str(CLMAX_takeoff,4),                 char(CLMAX_takeoff_unit)};
        cl_wb_landing   = {'CLMAX landing',      num2str(CLMAX_landing,4),                 char(CLMAX_landing_unit)};
        cd0             = {'CD0',                num2str(CD0,4),                           char(CD0_unit)};
        osw_eff_factor  = {'e',                  num2str(Oswald_efficiency_factor,4),      char(Oswald_efficiency_factor_un)};
        cl_0            = {'CL0 wb',             num2str(CL0,4),                           char(CL0_unit)};
        cl_alfa_rad     = {'CLALFA wb (rad)',    num2str(CLALFA_rad,4),                    char(CLALFA_rad_unit)};
        cl_alfa_deg     = {'CLALFA wb (deg)',    num2str(CLALFA_deg,4),                    char(CLALFA_deg_unit)};
        cd_gear         = {'CD wb landing gear', num2str(CDGEAR,4),                        char(CDGEAR_UNIT)};
        cm0             = {'CM0 wb',             num2str(CM0,4),                           char(CM0_unit)};
        cm_gear         = {'CM wb landing gear', num2str(CM_landing_gear,4),               char(CM_landing_gear_unit)};
        cmcl            = {'CMCL wb',            num2str(CMCL,4),                          char(CMCL_unit)};
        alfa0l          = {'Alfa wb zero lift',  num2str(Alpha_zero_lift,4),               char(Alpha_zero_lift_unit)};

        fieldValue = [cl_wb_max_clean; cl_wb_inv_clean; cl_wb_takeoff; ...
            cl_wb_landing; cl_0; cl_alfa_rad; cl_alfa_deg; cd0; cd_gear; ...
            cm0; cm_gear; cmcl; alfa0l; osw_eff_factor];
    
          
        tbl = FormalTable(header,fieldValue);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
        tbl = BaseTable(tbl);
        tbl.Title = strcat('Aerodynamic reference values.');
        tbl.LinkTarget = 'aerodynamictableRef';
        add(sec,tbl);
% -------------------------------------------------------------------------


add(ch,sec);


%% END chapter
%Adding chapters
add(rpt,ch);