%% CHAPTER 6 - Altitude
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

%% DECLARATION, DATA AND ASSUMPTIONTS CH5:

%reference AMC
requirement = strcat('-',Aircraft.Certification.Regulation.value);

op_altitude = Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.value;
op_altitude_un = Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.Attributes.unit;

%% END DECLARATION

ch = Chapter();
ch.Title = 'Altitude';
disp(['Chapter 6', (' "'), ch.Title,('" ') ,'writing...' ])

switch requirement
    
    % CASE 1: Very Light Aircraft
    case '-CSVLA'
        str = strcat('The maximum permissible operational altitude for the aircrat is ',...
            {' '},...
            num2str(op_altitude),...
            op_altitude_un,...
            '. Despite the',...
            requirement,...
            {' '},...
            'requirements do not require to accounts for the effects of altitude, ',...
            'such effects have been considered up to ',...
            {' '},...
            num2str(op_altitude),...
            op_altitude_un,...
            '. In fact the gust load factor have been ',...
            {' '},...
            ' calculated at such altitude. This is considered acceptable since it covers the operational ',...
            {' '},...
            ' range within which the aeroplane will fly most of the time.');
        %2
        str2 = strcat('(Note: the',...
            requirement,...
            ' requirement does not require to account for the effects of altitude. ',...
            ' Calculating the loads at sea level would be acceptable.',...
            {' '},...
            ' In this case, the choice to consider such effect up to ',...
            {' '},...
            num2str(op_altitude),...
            op_altitude_un,...
            {' '},...
            ' is a decision of a designer, which would be accepted by the team.)');
        
        %1
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(ch,para)
        %2
        para2 = Paragraph(str2);
        para2.Style = {HAlign('justify')};
        para2.BackgroundColor = "green";% = {Underline('yellow')};
        add(ch,para2)
% -------------------------------------------------------------------------
        %table gust calculation        
        str = ['TABLE TO BE CHECKED!!!'];
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        para.BackgroundColor = "red";
        add(ch,para);  
        
        % -----------------------------------------------------------------
        rho0     = Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value;
        rho_unit = Aircraft.Certification.ISA_Condition.Sea_Level.rho0.Attributes.unit;
        h0       = Aircraft.Certification.ISA_Condition.Sea_Level.Altitude.value;
        h_unit   = Aircraft.Certification.ISA_Condition.Sea_Level.Altitude.Attributes.unit;
        T0       = Aircraft.Certification.ISA_Condition.Sea_Level.T0.value;
        T_unit   = Aircraft.Certification.ISA_Condition.Sea_Level.T0.Attributes.unit;
        p0       = Aircraft.Certification.ISA_Condition.Sea_Level.p0.value;
        p_unit   = Aircraft.Certification.ISA_Condition.Sea_Level.p0.Attributes.unit;
        a0       = Aircraft.Certification.ISA_Condition.Sea_Level.Speed_of_sound0.value;
        a_unit   = Aircraft.Certification.ISA_Condition.Sea_Level.Speed_of_sound0.Attributes.unit;
        % -----------------------------------------------------------------
        header = {strcat('rho(',rho_unit,')'), strcat('h(',h_unit,')'), ...
           strcat('T(',T_unit,')'), strcat('p(',p_unit,')'), ... 
           strcat('a(', a_unit,')')};
        %each table row needs of a fieldValue
        %gust calculation table
        %1 
        sea_level = {num2str(rho0,4), num2str(h0,4), ...
                     num2str(T0,4), num2str(p0, '% 6.1f'), ...
                     num2str(a0,4)};    
%         PointA     = {'A', num2str(va,4), num2str(na,4)};     
%         PointC     = {'C', num2str(vc,4), num2str(nc,4)};        
%         PointD     = {'D', num2str(vd,4), num2str(nd,4)};
%         %2 
%         PointS_inv = {'S_inv', num2str(vs_inv,4), num2str(ns_inv,4)};    
%         PointG     = {'G', num2str(vg,4), num2str(ng,4)};     
%         PointF     = {'F', num2str(vf,4), num2str(nf,4)};        
%         PointE     = {'E', num2str(ve,4), num2str(ne,4)};
        fieldValue = [sea_level];
    
          
        tbl = FormalTable(header,fieldValue);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
        tbl = BaseTable(tbl);
        tbl.Title = strcat('Sea level reference values');
        tbl.LinkTarget = 'sealeveltable';
        add(ch,tbl);
% -------------------------------------------------------------------------
 
        % -----------------------------------------------------------------
        rho_op   = Aircraft.Certification.ISA_Condition.Operative_ceiling.rho0.value;
        h_op     = Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.value;
        T_op     = Aircraft.Certification.ISA_Condition.Operative_ceiling.T0.value;
        p_op     = Aircraft.Certification.ISA_Condition.Operative_ceiling.p0.value;
        a_op     = Aircraft.Certification.ISA_Condition.Operative_ceiling.Speed_of_sound0.value;
        % -----------------------------------------------------------------
        header = {strcat('rho(',rho_unit,')'), strcat('h(',h_unit,')'), ...
           strcat('T(',T_unit,')'), strcat('p(',p_unit,')'), ... 
           strcat('a(', a_unit,')')};
        %each table row needs of a fieldValue
        %gust calculation table
        %1 
        operative = {num2str(rho_op,4), num2str(h_op,4), ...
                     num2str(T_op,4), num2str(p_op, '% 6.1f'), ...
                     num2str(a_op,4)};    
%         PointA     = {'A', num2str(va,4), num2str(na,4)};     
%         PointC     = {'C', num2str(vc,4), num2str(nc,4)};        
%         PointD     = {'D', num2str(vd,4), num2str(nd,4)};
%         %2 
%         PointS_inv = {'S_inv', num2str(vs_inv,4), num2str(ns_inv,4)};    
%         PointG     = {'G', num2str(vg,4), num2str(ng,4)};     
%         PointF     = {'F', num2str(vf,4), num2str(nf,4)};        
%         PointE     = {'E', num2str(ve,4), num2str(ne,4)};
        fieldValue = [operative];
    
          
        tbl = FormalTable(header,fieldValue);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
        tbl = BaseTable(tbl);
        tbl.Title = strcat('Operative reference values');
        tbl.LinkTarget = 'operativeleveltable';
        add(ch,tbl);
% -------------------------------------------------------------------------

        % -----------------------------------------------------------------
        rho_t   = Aircraft.Certification.ISA_Condition.Theoretical_ceiling.rho0.value;
        h_t     = Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.value;
        T_t     = Aircraft.Certification.ISA_Condition.Theoretical_ceiling.T0.value;
        p_t     = Aircraft.Certification.ISA_Condition.Theoretical_ceiling.p0.value;
        a_t     = Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Speed_of_sound0.value;
        % -----------------------------------------------------------------
        header = {strcat('rho(',rho_unit,')'), strcat('h(',h_unit,')'), ...
           strcat('T(',T_unit,')'), strcat('p(',p_unit,')'), ... 
           strcat('a(', a_unit,')')};
        %each table row needs of a fieldValue
        %gust calculation table
        %1 
        theoretical = {num2str(rho_t,4), num2str(h_t,4), ...
                     num2str(T_t,4), num2str(p_t, '% 6.1f'), ...
                     num2str(a_t,4)};    
%         PointA     = {'A', num2str(va,4), num2str(na,4)};     
%         PointC     = {'C', num2str(vc,4), num2str(nc,4)};        
%         PointD     = {'D', num2str(vd,4), num2str(nd,4)};
%         %2 
%         PointS_inv = {'S_inv', num2str(vs_inv,4), num2str(ns_inv,4)};    
%         PointG     = {'G', num2str(vg,4), num2str(ng,4)};     
%         PointF     = {'F', num2str(vf,4), num2str(nf,4)};        
%         PointE     = {'E', num2str(ve,4), num2str(ne,4)};
        fieldValue = [theoretical];
    
          
        tbl = FormalTable(header,fieldValue);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
        tbl = BaseTable(tbl);
        tbl.Title = strcat('Theoretical ceiling reference values');
        tbl.LinkTarget = 'theoreticalleveltable';
        add(ch,tbl);
% -------------------------------------------------------------------------

        % CASE CS23
    case '-CS23'
        str = strcat('Check regulation to update ');
        str2 = strcat('Check regulation to update ');
        %1
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(ch,para)
        %2
        para2 = Paragraph(str2);
        para2.Style = {HAlign('justify')};
        add(ch,para2)
        
        % CASE 3: CS25
    case '-CS25'
        str = strcat('Check regulation to update ');
        str2 = strcat('Check regulation to update ');
        %1
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(ch,para)
        %2
        para2 = Paragraph(str2);
        para2.Style = {HAlign('justify')};
        add(ch,para2)

    % CASE CS22
    case '-CS22'
        str = strcat('Check regulation to update ');
        str2 = strcat('Check regulation to update ');
        %1
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(ch,para)
        %2
        para2 = Paragraph(str2);
        para2.Style = {HAlign('justify')};
        add(ch,para2)        

    % CASE CSLSA
    case '-CSLSA'
        str = strcat('Check regulation to update ');
        str2 = strcat('Check regulation to update ');
        %1
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(ch,para)
        %2
        para2 = Paragraph(str2);
        para2.Style = {HAlign('justify')};
        add(ch,para2)  
        
end
        


%% END chapter
%Adding chapters
add(rpt,ch);