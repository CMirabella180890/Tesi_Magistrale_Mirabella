%% CHAPTER 14 - Loads on the control surfaces
%Ref EASA/Tecnam p2006 report

import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

% CONSTANTS
requirement            = Aircraft.Certification.Regulation.value;
aileron_airworth_rules = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Aileron_load_SI.Attributes.cs;
% chapter_number = 12;
% ch = strcat('ch' , num2str(chapter_number)); 
ch = Chapter();
ch.Title = 'Loads on the control surfaces';
disp(['Chapter 14', (' "'), ch.Title,('" ') ,'writing...' ])
str = ['According to ' ...
       (' ') ...
       char(requirement) ...
       (' ') ...
       char(aileron_airworth_rules) ...
       (' ') ...
       ', the flight control system and its supporting structure must be' ...
       ' designed for loads corresponding to 125 % of the computed hinge' ...
       ' moments of the movable control surface.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(ch,para);

% +++++++
% AILERON
% +++++++
aileron_load_SI        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Aileron_load_SI.value;
aileron_load_SI_unit   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Aileron_load_SI.Attributes.unit;
aileron_airworth_rules = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Aileron_load_SI.Attributes.cs;
%sec
sec = Section();
sec.Title = 'Ailerons';
str = ['According to ' ...
       (' ') ...
       char(requirement) ...
       (' ') ...
       char(aileron_airworth_rules) ...
       (' ') ...
       ', the total aileron load is equal to' ...
       (' ') ...
       strcat(num2str(aileron_load_SI)) ...
       (' ') ...
       char(aileron_load_SI_unit) ...
       '. The hinge moment is calculated' ...
       ' by the following equation '];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);
% -------------------------------------------------------------------------
C_h_total_aileron = Aircraft.Certification.Aerodynamic_data.Hinge_moments.Aileron.C_h_total_deg.value;
qA                = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
S_aileron         = Aircraft.Geometry.Aileron.S.value;
cf                = Aircraft.Geometry.Aileron.cf.value;
HA                = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA.value;
HA_unit           = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.HA.Attributes.unit;
% -------------------------------------------------------------------------
% HA = C_h_total_deg * qA * S_aileron * cf;
        myNumEq = strcat( ' = ', ...
                                num2str(C_h_total_aileron,4), '*', ...
                                num2str(qA,4), '*', ...
                                num2str(S_aileron,4), '*', ...
                                num2str(cf,4), ' = ', ...
                                num2str(HA,4) , '\,\,', ...
                                HA_unit);
        myEq = "$ H_{aileron} = q*S_{aileron}*c_f*C_{h_{total}} ";
        eq = Equation(strcat(myEq, myNumEq));
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg);
% ------------------------------------------------------------------------- 
str = ['where: '];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);
% -------------------------------------------------------------------------
        % HA 
        myEq = "$ H_{aileron} = \mathrm{aileron hinge moment} (N * m) \mathrm{;} ";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg1 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref1 = eqImg1;         

        % q
        myEq = "$ q = \mathrm{dynamic pressure at point A} (Pa) \mathrm{;} ";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg2 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref2 = eqImg2; 
        
        % S_aileron 
        myEq = "$ S_{aileron} = \mathrm{aileron surface} (m^2) \mathrm{;} ";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg3 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref3 = eqImg3;         

        % cf
        myEq = "$ c_f = \mathrm{reference chord} (m) \mathrm{;}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg4 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref4 = eqImg4;      
        
        % C_h_total
        myEq = "$ C_{h_{total}} = \mathrm{total hinge moment coefficient.}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg5 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref5 = eqImg5; 
        
        ol = UnorderedList({ref1, ref2, ref3, ref4, ref5});
%         ol = UnorderedList({ref1,ref2,ref3,...
%             ref4,ref5,ref6, ref7,ref8});
%         ol = UnorderedList({ref1, ref2, ref3,...
%             ref4,ref5,ref6, ref7, ref8, ref9});
        append(sec,ol);
% -------------------------------------------------------------------------   
str = ['This is the formula used in all the following calculations.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);
% -------------------------------------------------------------------------
add(ch,sec);
% ++++++++
% ELEVATOR
% ++++++++
elevator_load_SI        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Total_elevator_SI.value;
elevator_load_SI_unit   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Total_elevator_SI.Attributes.unit;
%sec
sec = Section();
sec.Title = 'Elevator';
str = ['According to ' ...
       (' ') ...
       char(requirement) ...
       (' ') ...
       char(aileron_airworth_rules) ...
       (' ') ...
       ', the total elevator load is equal to' ...
       (' ') ...
       strcat(num2str(elevator_load_SI)) ...
       (' ') ...
       char(elevator_load_SI_unit) ...
       '.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);

add(ch,sec);
% ++++++
% RUDDER
% ++++++
rudder_load_SI        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.total_rudder_loads.value;
rudder_load_SI_unit   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.total_rudder_loads.Attributes.unit;
%sec
sec = Section();
sec.Title = 'Rudder';
str = ['According to ' ...
       (' ') ...
       char(requirement) ...
       (' ') ...
       char(aileron_airworth_rules) ...
       (' ') ...
       ', the total rudder load is equal to' ...
       (' ') ...
       strcat(num2str(rudder_load_SI)) ...
       (' ') ...
       char(rudder_load_SI_unit) ...
       '.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);

add(ch,sec);



%moving to another path for figure
% cd ..
% cd ..
%  regulation = Aircraft.Certification.Regulation.value;
%  results_path = [pwd '\' regulation '\Output\'];
% 
%  cd (RepDir);
% 
% fig = FormalImage([results_path,'Finalenvelope.pdf']);
%          fig.Caption = 'Maneuver and Gust load factors and diagram';
%          fig.Height = '5in';
%          fig.LinkTarget='maneuver_ref';
%          add(ch,fig);

%add(ch7,sec1);


%% END chapter
%Adding chapters
add(rpt,ch);