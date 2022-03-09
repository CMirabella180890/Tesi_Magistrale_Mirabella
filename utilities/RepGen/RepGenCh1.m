%% CHAPTER 1 - INTRODUCTION
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

requirement = strcat('-', Aircraft.Certification.Regulation.value);

ch = Chapter();
ch.Title = 'Introduction';
disp(['Chapter 1', (' "'), ch.Title,('" ') ,'writing...' ])
para = Paragraph();

append(para, Text(strcat('This document defines the SUBPART C - Structure - Flight Loads of the:',' ', char(Aircraft.Certification.Aircraft_Name.value),'.',...
    ' The boundaries of the flight envelope will be defined within this document. All speeds are calibrated airspeeds (CAS) (requirement 4.4 [1])',...
    ' and given in knots if not stated otherwise.',...
    ' All other units used are metric (SI units).',...
    ' The weights are given in mass units (kg) but the formulas require force units as input,',...
    ' therefore these are calculated in place wherever they are used.', ...
    ' Note: The speeds defined within this document should be used for the placards,',...
    ' speed markings, aeroplane flight manual (limitations), load calculations and need to be verified by flight test.')));
para.WhiteSpace = 'preserve';
para.Style = {HAlign('justify')};
add(ch,para)
%% STRUCTURE OF THIS DOCUMENT 
switch requirement
    % CASE 1: Very Light Aircraft
    case '-csvla'
% -------------------------------------------------------------------------
%sec
sec = Section();
sec.Title = 'Structure of this document';
str = ['This document is mainly focused on flight loads determination ' ...
    ' and structural sizing loads acting on the main wing and horizontal' ...
    ' empennage, according to' ...
       (' ') ...
       char(requirement(2:end)) ...
       (' ') ...
    ' airworthiness regulations. The structure of this document is derived' ...
    ' from the actual Subpart C subdivision of the regulations. These' ...
    ' regulations are applicable to aeroplanes with a single engine, with' ...
    ' no more than two seats, with a Maximum Certificated Takeoff Weight' ...
    ' of no more than 750 kg and a stalling speed in the landing configuration' ...
    ' of no more than 82 km/h (45 knots) calibrated airspeed (CAS); this aircraft' ...
    ' category is approved for daylight flight operations only, under Visual' ...
    ' Flight Rules (VFR). The subdivision of this document is summarized' ...
    ' as follows: '];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);
% -------------------------------------------------------------------------
%ordered list
ref1 = ['chapters 2 and 3 give to the reader chapter 4 reports a' ...
    ' summary of bibliographical information, list of abbreviation and symbols;' ...
    ' general data, geometrical arrengement, weight, mass and inertia, and' ...
    ' aerodynamic properties of the examined aircraft' ...
        ';'];
ref2 = ['chapter 5 lists the pertinent design airspeeds of the examined aircraft' ...
    ' according to the applied regulations; chapters 6 collects information about' ...
    ' standard atmosphere properties at sea level, operational and theoretical' ...
    ' ceiling for the examined aircraft;'];
ref3 = ['chapter 7 and 8 present all the information related to the flight, gust' ...
    ' and combined flight envelope diagrams of the examined aircraft, according to' ...
    ' the regulations applied;'];
ref4 = ['chapter 7 and 8 present all the information related to the flight, gust' ...
    ' and combined flight envelope diagrams of the examined aircraft, according to' ...
    ' the regulations applied;'];
ref5 = ['chapter 9 and 10 illustrate data and diagrams pertaining to actual loads' ...
    ' applied on the aircraft structural elements during flight;'];
ref6 = ['chapter 11, 12, 13 and 14 show manoeuvring and gust loads applied on aircraft' ...
    ' empennage, flaps, hinge and their structural support elements;'];
ref7 = ['chapter 15 illustrates loads applied on the engine and its support' ...
    ' structures.'];
ol = OrderedList({ref1, ref2, ref3, ref4, ref5, ref6, ref7});
append(sec,ol);
% -------------------------------------------------------------------------
%sec
sec = Section();
sec.Title = 'Exceptions and limitations';
str = ['The' ...
(' ') ...
char(requirement(2:end)) ...
(' ') ...
' airworthiness requirements are considered applicable to conventional aircraft.' ...
' The following non-conventional designs are not included:'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);
% -------------------------------------------------------------------------
%ordered list
ref1 = ['canards;'];
ref2 = ['tandem wings;'];
ref3 = ['winglets;'];
ref4 = ['other non-standard aerodynamic surfaces.'];
ol = OrderedList({ref1, ref2, ref3, ref4});
append(sec,ol);
% -------------------------------------------------------------------------
str = ['These design characteristics may need additional requirements and' ...
    ' analysis. The' ...
    (' ') ...
    char(requirement(2:end)) ...
    (' ') ...
    ' applies to aircraft intended for non-aerobatic operation only. This includes:'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);
% -------------------------------------------------------------------------
%ordered list
ref1 = ['any manoeuvre incident to normal flight;'];
ref2 = ['stalls, with the exception of whip stalls;'];
ref3 = ['Lazy eight, chandelles, steep turns in which the bank angle is no' ...
        ' more than sixty degrees.'];
ol = OrderedList({ref1, ref2, ref3});
append(sec,ol);
% -------------------------------------------------------------------------
add(ch,sec);
% -------------------------------------------------------------------------
    % CASE: CS23
    case '-cs-23'
% -------------------------------------------------------------------------
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
% -------------------------------------------------------------------------
    % CASE 3: CS25
    case '-cs-25'
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
% -------------------------------------------------------------------------
    % CASE CS22
    case '-cs-22'
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
% -------------------------------------------------------------------------
    % CASE CSLSA
    case '-cslsa'
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
% -------------------------------------------------------------------------
end
% add(ch,sec);

%% END chapter
%Adding chapters
add(rpt,ch);