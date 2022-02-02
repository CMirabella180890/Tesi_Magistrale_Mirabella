%% CHAPTER 9 - Loads on the aeroplane
%Ref EASA: WING LOAD CALCULATION
%(Example document for LSA applicants – v1 of 08.03.16)
%Date of issue: DD/MM/YYYY
%Document reference: ABCD-FL-57-00

import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

% (a) Strength requirements are specified in terms of limit loads (the maximum loads to be expected
% in service) and ultimate loads (limit loads multiplied by prescribed factors of safety). Unless
% otherwise provided, prescribed loads are limit loads.

% (b) Unless otherwise provided, the air, ground, and water loads must be placed in equilibrium with
% inertia forces, considering each item of mass in the aeroplane. These loads must be distributed
% to conservatively approximate or closely represent actual conditions.

% (c) If deflections under load would significantly change the distribution ofexternal or internal loads,
% this redistribution must be taken into account.

% (d) Simplified structural design criteria given in this Subpart C and its appendices may be used only
% for aeroplanes with conventional configurations. IfAppendixA is used, the entire appendix must
% be substituted for the corresponding paragraphs of this subpart, i.e. CS-VLA 321 to 459

% chapter_number = 9;
% ch = strcat('ch' , num2str(chapter_number)); 
ch = Chapter();
ch.Title = 'Loads on the aeroplane';
disp(['Chapter 9', (' "'), ch.Title,('" ') ,'writing...' ])

str = ['Strength requirements are specified in terms of limit loads' ...
    ' (the maximum loads to be expected in service) and ultimate loads' ...
    ' (limit loads multiplied by prescribed factors of safety). Unless' ...
    ' otherwise provided, prescribed loads are limit loads.' ...
    ' Unless otherwise provided, the air, ground, and water loads' ...
    ' must be placed in equilibrium with inertia forces, considering' ...
    ' each item of mass in the aeroplane. These loads must be' ...
    ' distributed to conservatively approximate or closely represent actual conditions.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(ch,para);

%sec
sec = Section();
sec.Title = 'Reference axes and sign convention';
str = ['In the figure is represented the reference frame used' ...
    ' to project forces and moment acting on the aircraft structures.' ...
    ' The origin is located at the airplane axis of symmetry' ...
    ' (x axis) with the y axis passing through the leading edge' ...
    ' of the mean aerodynamic chord section of the wing.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);
% -------------------------------------------------------------------------
%moving to another path for figure
cd ..
cd ..
 regulation = Aircraft.Certification.Regulation.value;
 results_path = [pwd '\' regulation '\Output\'];

 cd (RepDir);
 
 fig = FormalImage([results_path,'reference_axis.png']);
 fig.Caption = 'Reference axis';
 fig.Height = '5in';
 fig.LinkTarget='reference_axis';
 add(sec,fig);
% reference_axis
% -------------------------------------------------------------------------
% subsec
subsec = Section();
subsec.Title = 'Sign conventions and symbols';
str = ['Sign conventions and symbols used are summarized as follows:'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------
        % x
        myEq = "$ x = \mathrm{longitudinal axis of the aircraft;}";
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

        % y
        myEq = "$ y = \mathrm{lateral axis of the aircraft;}";
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
        
        % z
        myEq = "$ z = \mathrm{vertical axis of the aircraft;}";
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
        
        % Mx
        myEq = "$ M_x = \mathrm{total rolling moment;}";
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
        
        % My
        myEq = "$ M_y = \mathrm{total pitching moment;}";
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
        
        % Mz
        myEq = "$ M_z = \mathrm{total yawing moment;}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg6 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref6 = eqImg6;
        
        % Fx
        myEq = "$ F_x = \mathrm{toal axial force;}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg7 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref7 = eqImg7; 
        
        % Fy
        myEq = "$ F_y = \mathrm{total lateral force;}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg8 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref8 = eqImg8; 
        
        % Fz
        myEq = "$ F_z = \mathrm{total normal force;}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg9 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref9 = eqImg9;
        
        ol = UnorderedList({ref1, ref2, ref3, ref4, ref5, ref6, ref7, ref8, ref9});
%         ol = UnorderedList({ref1,ref2,ref3,...
%             ref4,ref5,ref6, ref7,ref8});
%         ol = UnorderedList({ref1, ref2, ref3,...
%             ref4,ref5,ref6, ref7, ref8, ref9});
        append(subsec,ol);
% ------------------------------------------------------------------------- 
add(sec,subsec);
% -------------------------------------------------------------------------
add(ch,sec);

%sec
sec = Section();
sec.Title = 'Symmetrical flight conditions';
str = ['The external forces and moments acting on the aeroplane in' ...
    ' a balanced flight condition have been determined. The' ...
    ' simplified scheme in figure is considered. The aeroplane' ...
    ' is reduced to the wing and the horizontal tail. The symbols' ...
    ' used are:'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%moving to another path for figure
cd ..
cd ..
 regulation = Aircraft.Certification.Regulation.value;
 results_path = [pwd '\' regulation '\Output\'];

 cd (RepDir);
 
 fig = FormalImage([results_path,'balance_reference.png']);
 fig.Caption = 'Simplified equilibrium of the aircraft';
 fig.Height = '3in';
 fig.LinkTarget='balance_reference';
 add(sec,fig);
% balance_reference
% -------------------------------------------------------------------------
        % xcg
        myEq = "$ x_{CG} = \mathrm{distance to aircraft centre of gravity;}";
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

        % xacwingfuselage
        myEq = "$ x_{AC_{f+w}} = \mathrm{distance to wing fuselage combination aerodynamic centre;}";
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
        
        % xPf+w
        myEq = "$ x_{P_{f+w}} = \mathrm{distance to wing fuselage combination centre of pressure;}";
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
        
        % nW
        myEq = "$ nW = \mathrm{aircraft total weight force;}";
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
        
        % xHT
        myEq = "$ x_{HT} = \mathrm{distance to HT quarter chord line;}";
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
% -------------------------------------------------------------------------
add(ch,sec);

%sec
sec = Section();
sec.Title = 'Aerodynamic centre';
str = ['ADD HERE details for balancing Equation'];
para = Paragraph(str);
add(ch,para);
add(ch,sec);

%sec
sec = Section();
sec.Title = 'Pitching moment of the wing-body';
str = ['ADD HERE details for balancing Equation'];
para = Paragraph(str);
add(ch,para);
add(ch,sec);

%moving to another path for figure
cd ..
cd ..
 regulation = Aircraft.Certification.Regulation.value;
 results_path = [pwd '\' regulation '\Output\'];

 cd (RepDir);
 
 fig = FormalImage([results_path,'Wingairloads.png']);
 fig.Caption = 'Wing airloads';
 fig.Height = '5in';
 fig.LinkTarget='wing_loads';
 add(sec,fig);

  fig = FormalImage([results_path,'Balancingloads.png']);
 fig.Caption = 'Balancing loads';
 fig.Height = '5in';
 fig.LinkTarget='bala_loads';
 add(sec,fig);

%sec
sec = Section();
sec.Title = 'Complete aircraft balancing loads';
str = ['ADD HERE table of results - LWB and LTAIL'];
para = Paragraph(str);
add(ch,para);
add(ch,sec);

%% END chapter
%Adding chapters
add(rpt,ch);