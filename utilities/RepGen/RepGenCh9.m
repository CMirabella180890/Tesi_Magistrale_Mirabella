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
str = ['Assuming positive forces and moment as depicted in the figures' ...
    ' it is possible to write the following equilibrium equations:'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);
% -------------------------------------------------------------------------
        % PITCHING ANGLE DIFFERENTIAL EQUATION 
        %
        myEq = "$ \mathrm{Z equilibrium: } \quad L_{wb}\cdot\cos{\alpha} + D_{wb}\cdot\sin{\alpha} + \frac{L_{t}}{\cos{\alpha}} - nW\cdot\cos{\alpha} = 0 ";
        eq = Equation(strcat(myEq));
        eq.DisplayInline = true;
        eq.FontSize = 10;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg); 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
        % -----------------------------------------------------------------
        % PITCHING ANGLE DIFFERENTIAL EQUATION 
        %
        myEq = "$ \mathrm{Rotation about } Y\mathrm{ : } ";
        eq = Equation(strcat(myEq));
        eq.DisplayInline = true;
        eq.FontSize = 10;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg); 
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        % PITCHING ANGLE DIFFERENTIAL EQUATION 
        %
        myEq = "$ [-L_{wb}\cdot\cos{\alpha}\cdot(x_{AC} - x_{CG}) - L_{wb}\cdot\sin{\alpha}\cdot b_{CG}] + [D_{wb}\cdot\cos{\alpha}\cdot b_{CG} - D_{wb}\cdot\sin{\alpha}\cdot c\cdot(x_{AC} - x_{CG})] ";
        eq = Equation(strcat(myEq));
        eq.DisplayInline = true;
        eq.FontSize = 10;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg); 
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        % PITCHING ANGLE DIFFERENTIAL EQUATION 
        %
        myEq = "$ + \biggl[-\frac{L_{t}}{\cos{\alpha}}\cdot(l_{t} + c\cdot(x_{AC} - x_{CG})) \biggr] - \\ T\cdot h + M_{wb} = 0 ";
        eq = Equation(strcat(myEq));
        eq.DisplayInline = true;
        eq.FontSize = 10;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg); 
        % -----------------------------------------------------------------
% -------------------------------------------------------------------------
add(ch,sec);
% -------------------------------------------------------------------------   
%sec
sec = Section();
sec.Title = 'Aerodynamic centre';
str = ['The aerodynamic centre of the wing-body aircraft was fixed at 25%' ...
    ' of the mean aerodynamic chord. This assumptions is use in all' ...
    ' calculations relative to forces and moments. Also, the maximum' ...
    ' takeoff weight of the aircraft has been assumed.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);
add(ch,sec);

%sec
sec = Section();
sec.Title = 'Pitching moment of the wing-body';
str = ['The pitching moment coefficient of the wing-body aircraft' ...
    ' was determined by high fidelity computational fluidynamics.' ...
    ' The results are depicted in the following figure.'];
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
 
 fig = FormalImage([results_path,'cfd_cm_results.png']);
 fig.Caption = 'Pitching moment coefficient from CFD calculations.';
 fig.Height = '4in';
 fig.LinkTarget='cfd_cm_results';
 add(sec,fig);
% reference_axis
% -------------------------------------------------------------------------
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
str = ['Here, the aircraft balancing loads are collected inside' ...
    ' table.'];
para = Paragraph(str);
add(sec,para);
% -------------------------------------------------------------------------
% SWITCH CASE
% -------------------------------------------------------------------------
switch (Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Straight_flight.value)
    % ---------------------------------------------------------------------
    % CASE 1
    % ---------------------------------------------------------------------
    case 'Case 1'
        n_gust_cruise_plus = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.value;
        nmax               = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value;
        % -----------------------------------------------------------------
        if max(n_gust_cruise_plus) > nmax
            % -----------------------------------------------------------------
            n1         = 1.0; 
            point_A    = "Point A";
            va         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
            va_unit    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.Attributes.unit;
            na         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
            na_unit    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.Attributes.unit;
            alfaA      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.value;
            alpha_unit = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.Attributes.unit;
            CLA        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value;
            LA         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA_new.value;
            L_unit     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA_new.Attributes.unit;
            LtailA     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.value;

            point_C    = "Point C";
            vc         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
            nc         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value;
            alfaC      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value;
            CLC        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value;
            LC         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.value;
            LtailC     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.value;

            point_D    = "Point D";
            vd         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
            nd         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value;
            alfaD      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value;
            CLD        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value;
            LD         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.value;
            LtailD     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value;

            point_G    = "Point G";
            vg         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG.value;
            ng         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.nG.value;
            alfaG      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.value;
            CLG        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value;
            LG         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG_new.value;
            LtailG     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.value;

            point_F    = "Point F";
            vf         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.value;
            nf         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.value;
            alfaF      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.value;
            CLF        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value;
            LF         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF_new.value;
            LtailF     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.value;

            point_E    = "Point E";
            ve         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.VE.value;
            ne         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.nE.value;
            alfaE      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.value;
            CLE        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.value;
            LE         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE_new.value;
            LtailE     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.value;
            % -----------------------------------------------------------------
            header = {'Point', strcat('V(',va_unit,')'), strcat('n(',na_unit,')'), ...
                strcat('alpha(',alpha_unit,')'), strcat('CL'), strcat('L(',L_unit,')'), strcat('L tail(',L_unit,')')};
            %each table row needs of a fieldValue
            %1
            name       = {char(point_A); char(point_C); char(point_D); ...
                          char(point_G); char(point_F); char(point_E)};
            speeds     = {num2str(va,4); num2str(vc,4); num2str(vd,4); ...
                          num2str(vg,4); num2str(vf,4); num2str(ve,4)};
            load_fact  = {num2str(na,4); num2str(nc,4); num2str(nd,4); ...
                          num2str(ng,4); num2str(nf,4); num2str(ne,4)};
            alfa       = {num2str(alfaA, 4); num2str(alfaC, 4); num2str(alfaD, 4); ...
                          num2str(alfaG, 4); num2str(alfaF, 4); num2str(alfaE, 4)};
            lift_coeff = {num2str(CLA, 4); num2str(CLC, 4); num2str(CLD, 4); ...
                          num2str(CLG, 4); num2str(CLF, 4); num2str(CLE, 4)};
            wing_lift  = {num2str(LA, 4); num2str(LC, 4); num2str(LD, 4); ...
                          num2str(LG, 4); num2str(LF, 4); num2str(LE, 4)};
            tail_lift  = {num2str(LtailA, 4); num2str(LtailC, 4); num2str(LtailD, 4); ...
                          num2str(LtailG, 4); num2str(LtailF, 4); num2str(LtailE, 4)};
            fieldValue = [name, speeds, load_fact, alfa, lift_coeff, wing_lift, tail_lift];


            tbl = FormalTable(header,fieldValue);
            % In order to put a table with a caption, the API Report denomination should
            % be used, the other options are from API DOM. In order to solve the problem,
            % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
            tbl = BaseTable(tbl);
            tbl.Title = strcat('Flight envelope points.');
            tbl.LinkTarget = 'flightenvelopepoints';
            add(sec,tbl);
            % -----------------------------------------------------------------            
        elseif max(n_gust_cruise_plus) < nmax 
             % -----------------------------------------------------------------
            n1         = 1.0; 
            point_A    = "Point A";
            va         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
            va_unit    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.Attributes.unit;
            na         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
            na_unit    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.Attributes.unit;
            alfaA      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.value;
            alpha_unit = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.Attributes.unit;
            CLA        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value;
            LA         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA_new.value;
            L_unit     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA_new.Attributes.unit;
            LtailA     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.value;

            point_C    = "Point C";
            vc         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
            nc         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value;
            alfaC      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value;
            CLC        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value;
            LC         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.value;
            LtailC     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.value;

            point_D    = "Point D";
            vd         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
            nd         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value;
            alfaD      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value;
            CLD        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value;
            LD         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.value;
            LtailD     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value;

            point_G    = "Point G";
            vg         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG.value;
            ng         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.nG.value;
            alfaG      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.value;
            CLG        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value;
            LG         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG_new.value;
            LtailG     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.value;

            point_F    = "Point F";
            vf         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.value;
            nf         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.value;
            alfaF      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.value;
            CLF        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value;
            LF         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF_new.value;
            LtailF     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.value;

            point_E    = "Point E";
            ve         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.VE.value;
            ne         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.nE.value;
            alfaE      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.value;
            CLE        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.value;
            LE         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE_new.value;
            LtailE     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.value;
            % -----------------------------------------------------------------
            header = {'Point', strcat('V(',va_unit,')'), strcat('n(',na_unit,')'), ...
                strcat('alpha(',alpha_unit,')'), strcat('CL'), strcat('L(',L_unit,')'), strcat('L tail(',L_unit,')')};
            %each table row needs of a fieldValue
            %1
            name       = {char(point_A); char(point_C); char(point_D); ...
                          char(point_G); char(point_F); char(point_E)};
            speeds     = {num2str(va,4); num2str(vc,4); num2str(vd,4); ...
                          num2str(vg,4); num2str(vf,4); num2str(ve,4)};
            load_fact  = {num2str(na,4); num2str(nc,4); num2str(nd,4); ...
                          num2str(ng,4); num2str(nf,4); num2str(ne,4)};
            alfa       = {num2str(alfaA, 4); num2str(alfaC, 4); num2str(alfaD, 4); ...
                          num2str(alfaG, 4); num2str(alfaF, 4); num2str(alfaE, 4)};
            lift_coeff = {num2str(CLA, 4); num2str(CLC, 4); num2str(CLD, 4); ...
                          num2str(CLG, 4); num2str(CLF, 4); num2str(CLE, 4)};
            wing_lift  = {num2str(LA, 4); num2str(LC, 4); num2str(LD, 4); ...
                          num2str(LG, 4); num2str(LF, 4); num2str(LE, 4)};
            tail_lift  = {num2str(LtailA, 4); num2str(LtailC, 4); num2str(LtailD, 4); ...
                          num2str(LtailG, 4); num2str(LtailF, 4); num2str(LtailE, 4)};
            fieldValue = [name, speeds, load_fact, alfa, lift_coeff, wing_lift, tail_lift];


            tbl = FormalTable(header,fieldValue);
            % In order to put a table with a caption, the API Report denomination should
            % be used, the other options are from API DOM. In order to solve the problem,
            % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
            tbl = BaseTable(tbl);
            tbl.Title = strcat('Flight envelope points.');
            tbl.LinkTarget = 'flightenvelopepoints';
            add(sec,tbl);
            % -----------------------------------------------------------------    
        end
    % ---------------------------------------------------------------------
    % CASE 2
    % ---------------------------------------------------------------------
    case 'Case 2'
        % -----------------------------------------------------------------
        n1         = 1.0; 
        point_A    = "Point A";
        va         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
        va_unit    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.Attributes.unit;
        na         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
        na_unit    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.Attributes.unit;
        alfaA      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.value;
        alpha_unit = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.Attributes.unit;
        CLA        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value;
        LA         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA_new.value;
        L_unit     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA_new.Attributes.unit;
        LtailA     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.value;

        point_C    = "Point C";
        vc         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
        nc         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value;
        alfaC      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value;
        CLC        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value;
        LC         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.value;
        LtailC     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.value;
        
        point_D    = "Point D";
        vd         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
        nd         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value;
        alfaD      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value;
        CLD        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value;
        LD         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.value;
        LtailD     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value;
        
        point_G    = "Point G";
        vg         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG.value;
        ng         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.nG.value;
        alfaG      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.value;
        CLG        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value;
        LG         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG_new.value;
        LtailG     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.value;
        
        point_F    = "Point F";
        vf         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.value;
        nf         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.value;
        alfaF      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.value;
        CLF        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value;
        LF         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF_new.value;
        LtailF     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.value;
        
        point_E    = "Point E";
        ve         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.VE.value;
        ne         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.nE.value;
        alfaE      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.value;
        CLE        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.value;
        LE         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE_new.value;
        LtailE     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.value;
        
        % -----------------------------------------------------------------
        header = {'Point', strcat('V(',va_unit,')'), strcat('n(',na_unit,')'), ...
            strcat('alpha(',alpha_unit,')'), strcat('CL'), strcat('L(',L_unit,')'), strcat('L tail(',L_unit,')')};
        %each table row needs of a fieldValue
        %1
        name       = {char(point_A); char(point_C); char(point_D); ...
                      char(point_G); char(point_F); char(point_E)};
        speeds     = {num2str(va,4); num2str(vc,4); num2str(vd,4); ...
                      num2str(vg,4); num2str(vf,4); num2str(ve,4)};
        load_fact  = {num2str(na,4); num2str(nc,4); num2str(nd,4); ...
                      num2str(ng,4); num2str(nf,4); num2str(ne,4)};
        alfa       = {num2str(alfaA, 4); num2str(alfaC, 4); num2str(alfaD, 4); ...
                      num2str(alfaG, 4); num2str(alfaF, 4); num2str(alfaE, 4)};
        lift_coeff = {num2str(CLA, 4); num2str(CLC, 4); num2str(CLD, 4); ...
                      num2str(CLG, 4); num2str(CLF, 4); num2str(CLE, 4)};
        wing_lift  = {num2str(LA, 4); num2str(LC, 4); num2str(LD, 4); ...
                      num2str(LG, 4); num2str(LF, 4); num2str(LE, 4)};
        tail_lift  = {num2str(LtailA, 4); num2str(LtailC, 4); num2str(LtailD, 4); ...
                      num2str(LtailG, 4); num2str(LtailF, 4); num2str(LtailE, 4)};
        fieldValue = [name, speeds, load_fact, alfa, lift_coeff, wing_lift, tail_lift];
    
          
        tbl = FormalTable(header,fieldValue);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
        tbl = BaseTable(tbl);
        tbl.Title = strcat('Flight envelope points.');
        tbl.LinkTarget = 'flightenvelopepoints';
        add(sec,tbl);
        % -----------------------------------------------------------------
% -------------------------------------------------------------------------
end
% -------------------------------------------------------------------------
% END SWITCH CASE
% -------------------------------------------------------------------------
add(ch,sec);

%% END chapter
%Adding chapters
add(rpt,ch);