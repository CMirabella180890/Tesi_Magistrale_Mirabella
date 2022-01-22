%% CHAPTER 1 - INTRODUCTION
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

ch = Chapter();
ch.Title = 'Introduction';
disp(['Chapter 1', (' "'), ch.Title,('" ') ,'writing...' ])
para = Paragraph();

append(para, Text(strcat('This document defines the SUBPART C - Structure - Flight Loads of the:',' ', char(Aircraft.Certification.Aircraft_Name.value),'.',...
    'The boundaries of the flight envelope will be defined within this document. All speeds are calibrated airspeeds (CAS) (requirement 4.4 [1])',...
    'and given in knots if not stated otherwise.',...
    'All other units used are metric (SI units).',...
    'The weights are given in mass units (kg) but the formulas require force units as input,',...
    'therefore these are calculated in place wherever they are used.', ...
    'Note: The speeds defined within this document should be used for the placards,',...
    'speed markings, aeroplane flight manual (limitations), load calculations and need to be verified by flight test.')));
para.WhiteSpace = 'preserve';
para.Style = {HAlign('justify')};
add(ch,para)
%% END chapter
%Adding chapters
add(rpt,ch);