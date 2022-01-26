%% CHAPTER 13 - Loads on the wing flaps
%Ref EASA/Tecnam p2006 report

import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

requirement = Aircraft.Certification.Regulation.value;
% CASE CS-VLA 441(a)(1)
case_a1 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.Attributes.cs;
% CASE CS-VLA 441(a)(2)
case_a2 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.Attributes.cs;
% CASE CS-VLA 441(a)(3)
case_a3 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.Attributes.cs;
% MAXIMUM RUDDER DEFLECTION - CONTROL STOPS
rudder_control_stops      = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.value;
rudder_control_stops_unit = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.Attributes.unit;
% LATERAL FORCE COEFFICIENT SLOPE PER DEGREE
dCYddeltar_per_degree = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CYdr.value;
% MAXIMU CY
CY_max_case_a1 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY.value;
% DYNAMIC PRESSURE AT VA
qA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
% VERTICAL TAIL WING SURFACE
S_vertical_total = 2*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value;
% VERTICAL TAIL WING SPAN 
b_vertical = Aircraft.Geometry.Vertical.b.value;
% MAX LATERAL FORCE
Y_max      = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_decanewton.value;
Y_max_unit = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_decanewton.Attributes.unit;
% CALCULATION OF (1.3)*(YAW_ANGLE)
Overswing_sideslip_angle = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.yaw_angle.value;
% LATERAL FORCE COEFFICIENT 
CY_VTP_case_a2 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_VTP.value;
% LATERAL FORCE CASE (a)(2) 
Y_case_a2 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.value;


% chapter_number = 12;
% ch = strcat('ch' , num2str(chapter_number)); 
ch = Chapter();
ch.Title = 'Loads on the wing flaps';
disp(['Chapter 13', (' "'), ch.Title,('" ') ,'writing...' ])

str = ['According to ' ...
    (' ') ...
    char(requirement) ...
    (' ') ...
    'the vertical tail must withstand several manoeuvring loads.' ...
    ' In this chapter, all these load case will be illustrated.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(ch,para);

%sec
sec = Section();
sec.Title = 'Manouevring load';
str = ['At speeds up to VA, the vertical tail surfaces must be designed to '...
       'withstand the following condition. In computing the tail loads, ' ...
       'the yawing velocity may be assumed to be zero. '];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);

%sub
subsec = Section();
subsec.Title = [char(requirement) ...
                char(case_a1)];
str = ['With the aeroplane in unaccelerated flight at zero yaw, ' ...
       'it is assumed that the rudder control is suddenly displaced ' ...
       'to the maximum deflection, as limited by the control stops or by limit pilot forces.' ...
       'The control stops are +/- ' ...
       (' ') ...
       strcat(num2str(rudder_control_stops)) ...
       (' ') ...
       char(rudder_control_stops_unit) ...
       '. The lateral force coefficient acting on the rudder when at '...
       'maximum deflection angle is given by the following simple ' ...
       ' equation: '];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------
        % LATERAL FORCE COEFFICIENT ON THE RUDDER AT MAX DEFLECTION
        % latex interprete with $ simbol
        myEq = "$ C_Y = C_{Y,0} + \frac{d C_Y}{d \delta_r} * \delta_{r,max} ";
        eq = Equation(strcat(myEq));
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(subsec,eqImg);
% -------------------------------------------------------------------------      
str = ['where: '];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------

        % CY
        myEq = "$ C_Y = \mathrm{lateral force coefficient;}";
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

        % C_Y_0
        myEq = "$ C_{Y, 0} = \mathrm{lateral force coefficient at } \beta=\delta=0 \mathrm{, equal to zero for symmetrical airfoil;}";
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
        
        % d C_Y / d delta
        myEq = "$ \frac{d C_Y }{d \delta_r } = \mathrm{lateral force curve slope per deg of rudder deflection;}";
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

        % delta max
        myEq = "$ \delta_{r, max} = \mathrm{rudder control stop.}";
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
        
        ol = UnorderedList({ref1, ref2, ref3, ref4});
%         ol = UnorderedList({ref1,ref2,ref3,...
%             ref4,ref5,ref6, ref7,ref8});
%         ol = UnorderedList({ref1, ref2, ref3,...
%             ref4,ref5,ref6, ref7, ref8, ref9});
        append(subsec,ol);
% -------------------------------------------------------------------------      
str = ['Assuming no deflection of the control cable (the control system ' ...
       'is infinitely rigid), the maximum value of the lateral force coefficient ' ...
       'is: '];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------
        % -----------------------------------------------------------------
        myNumEq = strcat( ' = ', num2str(dCYddeltar_per_degree,4), ' * ', ...
                          num2str(rudder_control_stops,4), ...
                          ' = ', num2str(CY_max_case_a1));
        % latex interprete with $ simbol
        myEq = "$ (C_Y)_{\delta_{r} = 30} ";
        eq = Equation(strcat(myEq, myNumEq));
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(subsec,eqImg);
        % -----------------------------------------------------------------
% -------------------------------------------------------------------------
str = ['The lateral force is calculated as follow: '];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------
        % -----------------------------------------------------------------
        myNumEq = strcat( ' = (1/10)*', num2str(qA,4), '*', ...
                          num2str(S_vertical_total,4), '*', ...
                          num2str(CY_max_case_a1,4), '/', ...
                          num2str(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value), ' = ', ...
                          num2str(Y_max,4), '\,\,', ...
                          Y_max_unit);
        % latex interprete with $ simbol
        myEq = "$ Y = \frac{1}{10}*q_A*S_{vertical}*\frac{(C_Y)_{\delta_{r} = 30}}{S_{ratio}} ";
        eq = Equation(strcat(myEq, myNumEq));
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(subsec,eqImg);
        % -----------------------------------------------------------------
% -------------------------------------------------------------------------
str = ['The lateral force acting on a single fin of the vertical tail plain is '...
    (' ') ...
    strcat(num2str(Y_max,4)) ...
    '/2 = ' ...
    (' ') ...
    strcat(num2str(Y_max*0.5,4)) ...
    (' ') ...
    char(Y_max_unit) ...
    '.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------
% SUB-SECTION CASE (a)(1)
add(sec,subsec);

%sub
subsec = Section();
subsec.Title = [char(requirement) ...
                char(case_a2)];
str = ['With the rudder deflected as specified in sub-paragraph '...
    (' ') ...
    char(requirement) ...
    (' ') ...
    char(case_a1) ...
    (' ') ...
    'of this paragraph, it is assumed that the aeroplane yaws to the ' ...
    'resulting sideslip angle. In lieu of a rational analysis, an ' ...
    'overswing angle equal to 1.3 times the static sideslip angle ' ...
    'of sub-paragraph ' ...
    (' ') ...
    char(requirement) ...
    (' ') ...
    char(case_a3) ...
    (' ') ...
    'of this paragraph may be assumed. The overswing sideslip angle is ' ...
    '1.3 * ' ...
    (' ') ...
    strcat(num2str(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.value)) ...
    (' ') ...
    ' = ' ...
    (' ') ...
    strcat(num2str(Overswing_sideslip_angle)) ...
    (' ') ...
    char(rudder_control_stops_unit) ...
    '. The total lateral force acting on the vertical tail in this case is: ' ];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------
        % -----------------------------------------------------------------
        myNumEq = strcat( ' = (1/10)*', num2str(qA,4), '*', ...
                          num2str(S_vertical_total,4), '*', ...
                          num2str(CY_VTP_case_a2,4), '/', ...
                          num2str(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value), ' = ', ...
                          num2str(Y_case_a2,4), '\,\,', ...
                          Y_max_unit);
        % latex interprete with $ simbol
        myEq = "$ Y =\frac{1}{10}*q_A*S_{vertical}*\frac{C_Y}{S_{ratio}} ";
        eq = Equation(strcat(myEq, myNumEq));
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(subsec,eqImg);
        % -----------------------------------------------------------------
% -------------------------------------------------------------------------
str = ['The lateral force acting on a single fin of the vertical tail plain is '...
    (' ') ...
    strcat(num2str(Y_case_a2,4)) ...
    '/2 = ' ...
    (' ') ...
    strcat(num2str(Y_case_a2*0.5,4)) ...
    (' ') ...
    char(Y_max_unit) ...
    '.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------        
% -------------------------------------------------------------------------
% SUB-SECTION CASE (a)(2)
add(sec,subsec);

%sub
subsec = Section();
subsec.Title = [char(requirement) ...
                char(case_a3)];
str = ['A yaw angle of 15 degrees with the rudder control maintained ' ...
       'in the neutral position (except as limited by pilot strength).' ...
       ' The total lateral force in this case is: ' ];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
        % -----------------------------------------------------------------
        myNumEq = strcat( ' = (1/10)*', num2str(qA,4), '*', ...
                          num2str(S_vertical_total,4), '*', ...
                          num2str(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_VTP.value,4), '/', ...
                          num2str(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value), ' = ', ...
                          num2str(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value,4), '\,\,', ...
                          Y_max_unit);
        % latex interprete with $ simbol
        myEq = "$ Y =\frac{1}{10}*q_A*S_{vertical}*\frac{C_Y}{S_{ratio}} ";
        eq = Equation(strcat(myEq, myNumEq));
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(subsec,eqImg);
        % -----------------------------------------------------------------
% -------------------------------------------------------------------------
str = ['The lateral force acting on a single fin of the vertical tail plain is '...
    (' ') ...
    strcat(num2str(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value,4)) ...
    '/2 = ' ...
    (' ') ...
    strcat(num2str(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value*0.5,4)) ...
    (' ') ...
    char(Y_max_unit) ...
    '.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------    
% SUB-SECTION CASE (a)(3)
add(sec,subsec);

add(ch,sec);

% %sec
% sec = Section();
% sec.Title = 'Balancing loads';
% str = ['ADD HERE details '];
% para = Paragraph(str);
% add(ch,para);
% add(ch,sec);

%sec
sec = Section();
sec.Title = 'Manouevring and gust envelope';
str = ['ADD HERE details '];
para = Paragraph(str);
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