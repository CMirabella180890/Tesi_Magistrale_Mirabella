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
op_altitude_un = Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.Attribute.unit;

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
            'Calculating the loads at sea level would be acceptable.',...
            {' '},...
            'In this case, the choice to consider such effect up to ',...
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