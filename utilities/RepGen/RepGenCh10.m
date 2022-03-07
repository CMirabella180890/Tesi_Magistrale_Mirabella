%% CHAPTER 9 - Loads on the wing
%Ref EASA: WING LOAD CALCULATION
%(Example document for LSA applicants – v1 of 08.03.16)
%Date of issue: DD/MM/YYYY
%Document reference: ABCD-FL-57-00

import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*
% -------------------------------------------------------------------------
ch = Chapter();
ch.Title = 'Loads on the wing';
disp(['Chapter 10', (' "'), ch.Title,('" ') ,'writing...' ])
% -------------------------------------------------------------------------
% moving to another path for figure
cd ..
cd ..
%  regulation = Aircraft.Certification.Regulation.value;
% results_path = [pwd '\utilities\Geometry\DroneVLA_results\'];
Check_name = convertCharsToStrings(Aircraft.Certification.Aircraft_Name.value);
if Check_name == "Drone_VLA"
    results_path =  [pwd '\utilities\Geometry\DroneVLA_results\'];
elseif Check_name == "TecnamP92"
    results_path =  [pwd '\utilities\Geometry\TecnamP92_results\'];
end
cd(RepDir);

fig = FormalImage([results_path,'Wing-Top-View.png']);
         fig.Caption = 'Wing, top view.';
         fig.Height = '5in';
         fig.LinkTarget='Wingtopview';
         add(ch,fig);
% -------------------------------------------------------------------------
str = ['In this section will be shown all the resulting internal' ...
    ' forces acting on the wing structural elements; having calculated' ...
    ' lift, drag and pitching moment coefficient distribution on the wing' ...
    ' with a panel method and the geometrical chord distribution, it is' ...
    ' possible to evaluate normal and shear forces and pitching moment' ...
    ' distributions along the wing span.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(ch,para);

% -------------------------------------------------------------------------
%sec
sec = Section();
sec.Title = 'Influence of the fuselage';
str = ['The effects of the fuselage on the wing span lift distribution' ...
    ' cause a reduction of lift at stations near the wing root; this lift' ...
    ' reduction can be discounted because is often negligible, leading to' ...
    ' a more conservative design loads. On the other hand, its influence on' ...
    ' the aeroplane equilibrium is accounted for, in particular on the' ...
    ' pitching moment distribution.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);

add(ch,sec);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%sec
sec = Section();
sec.Title = 'Forces and moments acting on the wings';
str = ['Numerical and graphical results from internal forces' ...
    ' and moments calculations will be shown in this section.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%sub
subsec = Section();
subsec.Title = 'SpanWise Airloads Distribution';
str = ['Spanwise airloads distributions along the wing semi-span' ...
    ' are obtained from a panel method; then, an interpolation through' ...
    ' all the values of the angle of attack is performed. Results are' ...
    ' represented in the following figures.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);

% moving to another path for figure
cd ..
cd ..
%  regulation = Aircraft.Certification.Regulation.value;
 results_path = [pwd '\csvla\Output\'];
cd(RepDir);

%to be checked
fig = FormalImage([results_path,'ClInterpolation3dplot.png']);
 fig.Caption = 'Wing lift coefficient spanwise distribution';
 fig.Height = '5in';
 fig.LinkTarget='wing_lift';
 add(subsec,fig);

 fig = FormalImage([results_path,'CdInterpolation3dplot.png']);
 fig.Caption = 'Wing drag coefficient spanwise distribution';
 fig.Height = '5in';
 fig.LinkTarget='wing_drag';
 add(subsec,fig);
 
 fig = FormalImage([results_path,'CmInterpolation3dplot.png']);
 fig.Caption = 'Wing pitching moment coefficient (0.25mac) spanwise distribution';
 fig.Height = '5in';
 fig.LinkTarget='wing_pitch';
 add(subsec,fig);
 
add(sec,subsec);

% -------------------------------------------------------------------------
%sub
subsec = Section();
subsec.Title = 'Normal and parallel component';
% -------------------------------------------------------------------------
%moving to another path for figure
cd ..
cd ..
 regulation = Aircraft.Certification.Regulation.value;
 results_path = [pwd '\' regulation '\Output\'];

 cd (RepDir);
 
 fig = FormalImage([results_path,'normalandaxial.png']);
 fig.Caption = 'Normal and axial forces decomposition.';
%  fig.Height = '4in';
 fig.LinkTarget='normalandaxial';
 add(subsec,fig);
 % ------------------------------------------------------------------------
 str = ['In order to determine the wing structural loads, the aerodynamic forces' ...
     ', Lift and Drag, are divided into their Normal and Parallel components in' ...
     ' respect of reference longitudinal axis, as depicted in the figure. The wing' ...
     ' spars and all structural component are normal or parallel to this axes which' ...
     ' the wing angle of attack is referred too.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------
        %
        myEq = "$ \mathrm{Normal and Axial decomposition: } ";
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
% -------------------------------------------------------------------------
        %
        myEq = "$ N_{wb} = L \cdot \cos{\alpha} + D \cdot \sin{\alpha} ";
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
% -------------------------------------------------------------------------
        %
        myEq = "$ A_{wb} = D \cdot \cos{\alpha} - L \cdot \sin{\alpha} ";
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
% -------------------------------------------------------------------------
 str = ['In non dimensionalized form, the decomposition become the following: '];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------
        %
        myEq = "$ C_{N} = C_{L} \cdot \cos{\alpha} + C_{D} \cdot \sin{\alpha} ";
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
        %
        myEq = "$ C_{A} = C_{D} \cdot \cos{\alpha} - C_{L} \cdot \sin{\alpha} ";
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
add(sec,subsec);

% ------------------------------------------------------------------------
%sub
subsec = Section();
subsec.Title = 'Shear, Bending and Torsion';
str = ['Shear, bending and torsion along the wing semi-span are' ...
    ' shown in the following figures; these distributions are also' ...
    ' reported inside a table, for each flight condition.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
%A
Check_name = convertCharsToStrings(Aircraft.Certification.Aircraft_Name.value);
if Check_name == "Drone_VLA"
    fig = FormalImage([results_path,'ShearBendingTorsionDiagramPointA.png']);
    fig.Caption = 'Shear, Bending and Torsion due to airloads - POINT A';
    fig.LinkTarget='A_distribution';
elseif Check_name == "TecnamP92"
    fig = FormalImage([results_path,'ShearBendingTorsionDiagramPointA1.png']);
    fig.Caption = 'Shear, Bending and Torsion due to airloads - POINT A1';
    fig.LinkTarget='A_distribution';
end
% fig = FormalImage([results_path,'ShearBendingTorsionDiagramPointA.png']);
%  fig.Caption = 'Shear, Bending and Torsion due to airloads - POINT A';
%  fig.Height = '4in';
%  fig.LinkTarget='A_distribution';
 add(subsec,fig);

 %C
fig = FormalImage([results_path,'ShearBendingTorsionDiagramPointC.png']);
 fig.Caption = 'Shear, Bending and Torsion due to airloads - POINT C';
 fig.Height = '4in';
 fig.LinkTarget='C_distribution';
 add(subsec,fig);

  %D
fig = FormalImage([results_path,'ShearBendingTorsionDiagramPointD.png']);
 fig.Caption = 'Shear, Bending and Torsion due to airloads - POINT D';
 fig.Height = '4in';
 fig.LinkTarget='D_distribution';
 add(subsec,fig);
% ------------------------------------------------------------------------
 str = ['To evaluate the internal forces needed to the structural' ...
     ' sizing of the aircraft wing, a summation along the wing semi' ...
     ' span is performed. First, axial and normal components of the' ...
     ' aerodynamic forces are calculated from the semi-span distributions' ...
     ' of lift and drag coefficients, which are used to evaluate shear,' ...
     ' bending and torsion dimensional distributions.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------
%                 SHEAR
%                 a = 0.5*(FZ(i-1)+FZ(i))*(y(i-1) - y(i));
%                 S(i) = S(i-1) + a; 
% -------------------------------------------------------------------------
                    %
                    myEq = "$ S_{i} = S_{i-1} + s ";
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
 str = ['where'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% ------------------------------------------------------------------------- 
                    %
                    myEq = "$ s = \frac{1}{2} \cdot (N_{i-1} + N_{i})\cdot(y_{i-1} - y_{i}) ";
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
 str = ['where N is the normal force acting on the wing at the span station y.' ...
     ' For the bending moment distribution, the calculations are very similar:' ...
     ' Torsion is calculated via the same summation, but with torsion couple' ...
     ' along the wing semi-span:'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------
%                 BENDING
                    %
                    myEq = "$ BM_{i} = BM_{i-1} + b ";
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
                    %                 BENDING
                    % a = 0.5*(S(i-1)+S(i))*(y(i-1) - y(i));
                    myEq = "$ b = \frac{1}{2} \cdot (S_{i-1} + S_{i}) \cdot (y_{i-1} - y_{i}) ";
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
% -------------------------------------------------------------------------
str = ['Torsion is calculated via the same summation, but with torsion couple' ...
     ' along the wing semi-span:'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------
%                 TORSION
%                 x         = c_of_y(i);
%                 M_of_y(i) = CM(i)*Dyn_Press*(x^2);
                    %
                    myEq = "$ M_{y_{i}} = C_{M_{i}}\cdot q \cdot (c_{y_{i}})^{2} ";
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
 str = ['where the chords are evaluated at every semi-span station.' ...
     ' The summation is performed in the same way as did with shear and' ...
     ' bending distributions.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% ------------------------------------------------------------------------- 
        % -----------------------------------------------------------------
        n1           = 1.0; 
        
        point_S    = "Point S";
        vs         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.value;
        ns         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.nS.value;
        shearS     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_distr.value(end);
        bendingS   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Bend_mom_distr.value(end);
        torsionS   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Tors_mom_distr.value(end);
        
        point_A      = "Point A";
        va           = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
        va_unit      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.Attributes.unit;
        na           = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
        na_unit      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.Attributes.unit;
        if isfield(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA, 'Shear_distr') == 1
            shearA       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_distr.value(end);
            shear_unit   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_distr.Attributes.unit;
            bendingA     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Bend_mom_distr.value(end);
            bending_unit = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Bend_mom_distr.Attributes.unit;
            torsionA     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Tors_mom_distr.value(end);
            torsion_unit = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Tors_mom_distr.Attributes.unit;
        elseif isfield(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope, 'PointA1') == 1
            shearA       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Shear_distr.value(end);
            shear_unit   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Shear_distr.Attributes.unit;
            bendingA     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Bend_mom_distr.value(end);
            bending_unit = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Bend_mom_distr.Attributes.unit;
            torsionA     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Tors_mom_distr.value(end);
            torsion_unit = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Tors_mom_distr.Attributes.unit;
        end

        point_C    = "Point C";
        vc         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
        nc         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value;
        shearC     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_distr.value(end);
        bendingC   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Bend_mom_distr.value(end);
        torsionC   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Tors_mom_distr.value(end);
        
        point_D    = "Point D";
        vd         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
        nd         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value;
        shearD     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_distr.value(end);
        bendingD   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Bend_mom_distr.value(end);
        torsionD   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Tors_mom_distr.value(end);
        
        point_Sinv  = "Point Sinv";
        vs_inv      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.VS_inv.value;
        ns_inv      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.nS_inv.value;
        shearSinv   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.Shear_distr.value(end);
        bendingSinv = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.Bend_mom_distr.value(end);
        torsionSinv = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.Tors_mom_distr.value(end);
        
        point_G    = "Point G";
        vg         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG.value;
        ng         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.nG.value;
        alfaG      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alfaG.value;
        CLG        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value;
        LG         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG_new.value;
        LtailG     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.value;
        shearG     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Shear_distr.value(end);
        bendingG   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Bend_mom_distr.value(end);
        torsionG   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Tors_mom_distr.value(end);
        
        point_F    = "Point F";
        vf         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.value;
        nf         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.value;
        alfaF      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alfaF.value;
        CLF        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value;
        LF         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF_new.value;
        LtailF     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.value;
        shearF     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Shear_distr.value(end);
        bendingF   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Bend_mom_distr.value(end);
        torsionF   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Tors_mom_distr.value(end);
        
        point_E    = "Point E";
        ve         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.VE.value;
        ne         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.nE.value;
        alfaE      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alfaE.value;
        CLE        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.value;
        LE         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE_new.value;
        LtailE     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.value;
        shearE     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Shear_distr.value(end);
        bendingE   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Bend_mom_distr.value(end);
        torsionE   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Tors_mom_distr.value(end);
        % -----------------------------------------------------------------
        check_shear = [abs(shearS),   abs(shearA),   abs(shearC),   abs(shearD),   abs(shearSinv),   abs(shearG),   abs(shearF),   abs(shearE)  ];
        check_bend  = [abs(bendingS), abs(bendingA), abs(bendingC), abs(bendingD), abs(bendingSinv), abs(bendingG), abs(bendingF), abs(bendingE)];
        check_tors  = [abs(torsionS), abs(torsionA), abs(torsionC), abs(torsionD), abs(torsionSinv), abs(torsionG), abs(torsionF), abs(torsionE)];
      % -----------------------------------------------------------------
        if max(check_shear) == abs(shearS)
            critical_point_shear   = 'Point S';
            critical_shear         = shearS;
            string_paragraph_shear = [' The critical shear is at point S and is' ...
                ' equal to' ...
                (' ') ...
                num2str(shearS, 4) ...
                (' ') ...
                char(shear_unit) ...
                '.']; 
        elseif max(check_shear) == abs(shearA)
            critical_point_shear   = 'Point A';
            critical_shear         = shearA;
            string_paragraph_shear = [' The critical shear is at point A and is' ...
                ' equal to' ...
                (' ') ...
                num2str(shearA, 4) ...
                (' ') ...
                char(shear_unit) ...
                '.']; 
        elseif max(check_shear) == abs(shearC)
            critical_point_shear   = 'Point C';
            critical_shear         = shearC;
            string_paragraph_shear = [' The critical shear is at point C and is' ...
                ' equal to' ...
                (' ') ...
                num2str(shearC, 4) ...
                (' ') ...
                char(shear_unit) ...
                '.'];
        elseif max(check_shear) == abs(shearD)
            critical_point_shear   = 'Point D';
            critical_shear         = shearD;
            string_paragraph_shear = [' The critical shear is at point D and is' ...
                ' equal to' ...
                (' ') ...
                num2str(shearC, 4) ...
                (' ') ...
                char(shear_unit) ...
                '.'];
        elseif max(check_shear) == abs(shearSinv)
            critical_point_shear   = 'Point Sinv';
            critical_shear         = shearSinv;
            string_paragraph_shear = [' The critical shear is at point Sinv and is' ...
                ' equal to' ...
                (' ') ...
                num2str(shearSinv, 4) ...
                (' ') ...
                char(shear_unit) ...
                '.'];
        elseif max(check_shear) == abs(shearG)
            critical_point_shear   = 'Point G';
            critical_shear         = shearG;
            string_paragraph_shear = [' The critical shear is at point G and is' ...
                ' equal to' ...
                (' ') ...
                num2str(shearG, 4) ...
                (' ') ...
                char(shear_unit) ...
                '.'];
        elseif max(check_shear) == abs(shearF)
            critical_point_shear   = 'Point F';
            critical_shear         = shearF;
            string_paragraph_shear = [' The critical shear is at point F and is' ...
                ' equal to' ...
                (' ') ...
                num2str(shearF, 4) ...
                (' ') ...
                char(shear_unit) ...
                '.'];
        elseif max(check_shear) == abs(shearE)
            critical_point_shear   = 'Point E';
            critical_shear         = shearE;
            string_paragraph_shear = [' The critical shear is at point E and is' ...
                ' equal to' ...
                (' ') ...
                num2str(shearE, 4) ...
                (' ') ...
                char(shear_unit) ...
                '.'];
        end
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        if max(check_bend) == abs(bendingS)
            critical_point_bend   = 'Point S';
            critical_bend         = bendingS;
            string_paragraph_bend = [' The critical bending is at point S and is' ...
                ' equal to' ...
                (' ') ...
                num2str(bendingS, 4) ...
                (' ') ...
                char(bending_unit) ...
                '.']; 
        elseif max(check_bend) == abs(bendingA)
            critical_point_bend   = 'Point A';
            critical_bend         = bendingA;
            string_paragraph_bend = [' The critical bending is at point A and is' ...
                ' equal to' ...
                (' ') ...
                num2str(bendingA, 4) ...
                (' ') ...
                char(bending_unit) ...
                '.']; 
        elseif max(check_bend) == abs(bendingC)
            critical_point_bend   = 'Point C';
            critical_bend         = bendingC;
            string_paragraph_bend = [' The critical bending is at point C and is' ...
                ' equal to' ...
                (' ') ...
                num2str(bendingC, 4) ...
                (' ') ...
                char(bending_unit) ...
                '.'];
        elseif max(check_bend) == abs(bendingD)
            critical_point_bend   = 'Point D';
            critical_bend         = bendingD;
            string_paragraph_bend = [' The critical bending is at point D and is' ...
                ' equal to' ...
                (' ') ...
                num2str(bendingD, 4) ...
                (' ') ...
                char(bending_unit) ...
                '.'];
        elseif max(check_bend) == abs(bendingSinv)
            critical_point_bend   = 'Point Sinv';
            critical_bend         = bendingSinv;
            string_paragraph_bend = [' The critical bending is at point Sinv and is' ...
                ' equal to' ...
                (' ') ...
                num2str(bendingSinv, 4) ...
                (' ') ...
                char(bending_unit) ...
                '.'];
        elseif max(check_bend) == abs(bendingG)
            critical_point_bend   = 'Point G';
            critical_bend         = bendingG;
            string_paragraph_bend = [' The critical bending is at point G and is' ...
                ' equal to' ...
                (' ') ...
                num2str(bendingG, 4) ...
                (' ') ...
                char(torsion_unit) ...
                '.'];
        elseif max(check_bend) == abs(bendingF)
            critical_point_bend   = 'Point F';
            critical_bend         = bendingF;
            string_paragraph_bend = [' The critical bending is at point F and is' ...
                ' equal to' ...
                (' ') ...
                num2str(bendingF, 4) ...
                (' ') ...
                char(bending_unit) ...
                '.'];
        elseif max(check_bend) == abs(bendingE)
            critical_point_bend   = 'Point E';
            critical_bend         = bendingE;
            string_paragraph_bend = [' The critical bending is at point E and is' ...
                ' equal to' ...
                (' ') ...
                num2str(bendingE, 4) ...
                (' ') ...
                char(bending_unit) ...
                '.'];
        end
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        if max(check_tors) == abs(torsionS)
            critical_point_tors   = 'Point S';
            critical_tors         = torsionS;
            string_paragraph_tors = [' The critical torsion is at point S and is' ...
                ' equal to' ...
                (' ') ...
                num2str(torsionS, 4) ...
                (' ') ...
                char(torsion_unit) ...
                '.']; 
        elseif max(check_tors) == abs(torsionA)
            critical_point_tors   = 'Point A';
            critical_tors         = torsionA;
            string_paragraph_tors = [' The critical torsion is at point A and is' ...
                ' equal to' ...
                (' ') ...
                num2str(torsionA, 4) ...
                (' ') ...
                char(torsion_unit) ...
                '.']; 
        elseif max(check_tors) == abs(torsionC)
            critical_point_tors   = 'Point C';
            critical_tors         = torsionC;
            string_paragraph_tors = [' The critical torsion is at point C and is' ...
                ' equal to' ...
                (' ') ...
                num2str(torsionC, 4) ...
                (' ') ...
                char(torsion_unit) ...
                '.'];
        elseif max(check_tors) == abs(torsionD)
            critical_point_tors   = 'Point D';
            critical_tors         = torsionD;
            string_paragraph_tors = [' The critical torsion is at point D and is' ...
                ' equal to' ...
                (' ') ...
                num2str(torsionD, 4) ...
                (' ') ...
                char(torsion_unit) ...
                '.'];
        elseif max(check_tors) == abs(torsionSinv)
            critical_point_tors   = 'Point Sinv';
            critical_tors         = torsionSinv;
            string_paragraph_tors = [' The critical torsion is at point Sinv and is' ...
                ' equal to' ...
                (' ') ...
                num2str(torsionSinv, 4) ...
                (' ') ...
                char(torsion_unit) ...
                '.'];
        elseif max(check_tors) == abs(torsionG)
            critical_point_tors   = 'Point G';
            critical_tors         = torsionG;
            string_paragraph_tors = [' The critical torsion is at point G and is' ...
                ' equal to' ...
                (' ') ...
                num2str(torsionG, 4) ...
                (' ') ...
                char(torsion_unit) ...
                '.'];
        elseif max(check_tors) == abs(torsionF)
            critical_point_tors   = 'Point F';
            critical_tors         = torsionF;
            string_paragraph_tors = [' The critical torsion is at point F and is' ...
                ' equal to' ...
                (' ') ...
                num2str(torsionF, 4) ...
                (' ') ...
                char(torsion_unit) ...
                '.'];
        elseif max(check_tors) == abs(torsionE)
            critical_point_tors   = 'Point E';
            critical_tors         = torsionE;
            string_paragraph_tors = [' The critical torsion is at point E and is' ...
                ' equal to' ...
                (' ') ...
                num2str(torsionE, 4) ...
                (' ') ...
                char(torsion_unit) ...
                '.'];
        end
        % -----------------------------------------------------------------        
        header = {'Point', strcat('V (',va_unit,')'), strcat('n (',na_unit,')'), ...
            strcat('S (',shear_unit,')'), strcat('M (',bending_unit,')'), strcat('T (',torsion_unit,')')};
        %each table row needs of a fieldValue
        %1
        name       = {char(point_S); char(point_A); char(point_C); ...
                      char(point_D); char(point_Sinv); char(point_G); ...
                      char(point_F); char(point_E)};
        speeds     = {num2str(vs,4); num2str(va,4); num2str(vc,4); num2str(vd,4); ...
                      num2str(vs_inv,4); num2str(vg,4); num2str(vf,4); num2str(ve,4)};
        load_fact  = {num2str(ns,4); num2str(na,4); num2str(nc,4); num2str(nd,4); ...
                      num2str(ns_inv,4); num2str(ng,4); num2str(nf,4); num2str(ne,4)};
        shear      = {num2str(shearS, 4); num2str(shearA, 4); num2str(shearC, 4); num2str(shearD, 4); ... 
                      num2str(shearSinv, 4); num2str(shearG, 4); num2str(shearF, 4); num2str(shearE, 4)};
        bending    = {num2str(bendingS, 4); num2str(bendingA, 4); num2str(bendingC, 4); num2str(bendingD, 4); ...
                      num2str(bendingSinv, 4); num2str(bendingG, 4); num2str(bendingF, 4); num2str(bendingE, 4)};
        torsion    = {num2str(torsionS, 4); num2str(torsionA, 4); num2str(torsionC, 4); num2str(torsionD, 4); ...
                      num2str(torsionSinv, 4); num2str(torsionG, 4); num2str(torsionF, 4); num2str(torsionE, 4)};
        fieldValue = [name, speeds, load_fact, shear, bending, torsion];
          
        tbl = FormalTable(header,fieldValue);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
        tbl = BaseTable(tbl);
        tbl.Title = strcat('Shear, bending and torsion distribution along the semi-span.');
        tbl.LinkTarget = 'shearbendingtors_distr';
        add(subsec,tbl);
        % -----------------------------------------------------------------
        str = ['From the table is evident that Point A is critical for shear and' ...
               ' bending, while torsion is critical at points D and E.'];
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(subsec,para);
% ------------------------------------------------------------------------
add(sec,subsec);

% -------------------------------------------------------------------------
%sub
subsec = Section();
subsec.Title = 'Critical load condition';
 str = ['Critical values of internal forces and moments are relative to' ...
     ' the root station of the wing.'];
%  str = ['Critical values of internal forces and moments are relative to' ...
%      ' the root station of the wing. Shear is critical at point C, where is equal to' ... 
%      (' ') ...
%      num2str(shearC, 4) ...
%      (' ') ...
%      char(shear_unit) ...
%      '; the corresponding inverted flight point F is also' ...
%      ' a point where shear is higher than other flight conditions in inverted flight.' ...
%      ' At point C the bending moment at the root is equal to' ...
%      (' ') ...
%      num2str(bendingC, 4) ...
%      (' ') ...
%      char(bending_unit) ...
%      ', while bending is not critical for inverted flight.' ...
%      ' Torsion is critical at points D and E, where we have' ...
%      (' ') ...
%      num2str(torsionD, 4) ...
%      (' ') ...
%      char(torsion_unit) ...
%      ' and ' ...
%      (' ') ...
%      num2str(torsionE, 4) ...
%      (' ') ...
%      char(torsion_unit) ...
%      ', respectively.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %ordered list
        ref1 = string_paragraph_shear;
        ref2 = string_paragraph_bend;
        ref3 = string_paragraph_tors;
        ol = UnorderedList({ref1, ref2, ref3});
        append(subsec,ol);
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        % -----------------------------------------------------------------
        header = {'Point', 'Value', 'Units', 'Load'};
        %each table row needs of a fieldValue
        %1
        name       = {char(critical_point_shear); char(critical_point_bend); char(critical_point_tors)};
        values     = {num2str(critical_shear, 4); num2str(critical_bend, 4); num2str(critical_tors, 4)};
        units      = {char(shear_unit); char(bending_unit); char(torsion_unit)};
        load       = {'Shear'; 'Bending'; 'Torsion'};
        fieldValue = [name, values, units, load];
          
        tbl = FormalTable(header,fieldValue);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
        tbl = BaseTable(tbl);
        tbl.Title = strcat('Critical shear, bending and torsion.');
        tbl.LinkTarget = 'shearbendingtors_crit';
        add(subsec,tbl);
        % -----------------------------------------------------------------
 str = ['In the table these values are summarized, for convenience.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);        
% ------------------------------------------------------------------------- 

%SHEAR
fig = FormalImage([results_path,'ShearComparison.png']);
fig.Caption = 'Shear comparison';
fig.Height = '5in';
fig.LinkTarget='Bending';
add(subsec,fig);

%Bending
fig = FormalImage([results_path,'BendingComparison.png']);
fig.Caption = 'Bending comparison';
fig.Height = '5in';
fig.LinkTarget='Bending';
add(subsec,fig);

%Torsion
fig = FormalImage([results_path,'TorsionComparison.png']);
fig.Caption = 'Torsion comparison';
fig.Height = '5in';
fig.LinkTarget='Torsion';
add(subsec,fig);
add(sec,subsec);
% -------------------------------------------------------------------------
add(ch,sec);

% -------------------------------------------------------------------------
%sec
requirement         = Aircraft.Certification.Regulation.value;
Unsymm_req_airworth = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.cs;
%   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%   CS-VLA 349 Rolling conditions 
%   The wing and wing bracing must be designed for the following loading
%   conditions:   
%   (a) Unsymmetrical wing loads. Unless the following values result in 
%       unrealistic loads, the rolling accelerations may be obtained by
%       modifying the symmetrical flight conditions in CS-VLA 333(d) as
%       follows: In condition A, assume that 100% of the semispan wing
%       airload acts on one side of the aeroplane and 70% of this load
%       acts on the other side.  
%   (b) The  loads  resulting  from  the  aileron  deflections  and  
%       speeds  specified  in  CS-VLA 455, in combination with an aero- 
%       plane load factor of at least two thirds of the positive 
%       manoeuvring load factor used for design. Unless the following 
%       values result in unrealistic loads, the effect of aileron
%       displacement on wing torsion may be accounted for by adding the 
%       following increment to  the basic aerofoil moment coefficient
%       over the aileron portion of  the span in  the critical condition 
%       determined in CS-VLA 333(d):
%       
%                      DELTA_CM = (-0.01)*DELTA_AILERON
%
%       with 
%      
%       DELTA_CM      --> Moment coefficient increment
%       DELTA_AILERON --> Down aileron deflection in degrees in the 
%                         critical condition
%       
%       NOTE: The angle at critical condition DELTA_AILERON must be given
%             in degrees. 
%   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sec = Section();
sec.Title = 'Unsymmetrical loads';
str = ['According to' ...
    (' ') ...
    char(requirement) ...
    (' ') ...
    char(Unsymm_req_airworth) ...
    ', the wing and wing bracing must be designed for the following' ...
    ' loading conditions:'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);
% -------------------------------------------------------------------------
flight_envelope_cs   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.Attributes.cs;
aileron_req_airworth = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.cs;
        %ordered list
        ref1 = ['Unsymmetrical wing loads. Unless the following values result in' ...
      ' unrealistic loads, the rolling accelerations may be obtained by'  ...
      ' modifying the symmetrical flight conditions in' ...
      (' ') ...
      char(requirement) ...
      (' ') ...
      char(flight_envelope_cs) ...
      (' ') ...
      ' follows: in condition A, assume that 100% of the semispan wing' ...
      ' airload acts on one side of the aeroplane and 70% of this load' ...
      ' acts on the other side.' ];
        ref2 = ['Aileron deflection. The  loads  resulting  from  the  aileron  deflections' ...
            ' and speeds  specified  in'...
            (' ') ...
            char(requirement) ...
            (' ') ...
            char(aileron_req_airworth) ...
            ', in combination with an aeroplane load factor of at least' ...
            ' two thirds of the positive manoeuvring load factor used for' ...
            ' design. Unless the following values result in unrealistic' ...
            ' loads, the effect of aileron displacement on wing torsion' ...
            ' may be accounted for by adding the following increment to' ...
            ' the basic aerofoil moment coefficient over the aileron' ...
            ' portion of  the span in  the critical condition determined in' ...
            (' ') ...
            char(requirement) ...
            (' ') ...
            char(flight_envelope_cs) ...
            '.'];
        ol = OrderedList({ref1, ref2});
        append(sec,ol);
        % -----------------------------------------------------------------
% -------------------------------------------------------------------------
Global_torsion_full_load_airwort_reg = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.cs;
Flight_envelope_reg_airworth         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.Attributes.cs;
%sub
subsec = Section();
subsec.Title = 'Rolling condition';
str = ['According to' ...
    (' ') ...
    char(requirement) ...
    (' ') ...
    char(Global_torsion_full_load_airwort_reg) ...
    ', the aileron displacement cause a significant change of' ...
    ' pitching moment distribution along the wing span; these changes are' ...
    ' shown in the following diagrams for different fligh condition' ...
    ' in red. Unless the following values result in unrealistic loads,' ...
    ' the effect of aileron displacement on wing torsion may be' ...
    ' accounted for by adding the following increment to  the basic' ...
    ' aerofoil moment coefficient over the aileron portion of  the' ...
    ' span in  the critical condition determined in' ...
    (' ') ...
    char(requirement) ...
    (' ') ...
    char(Flight_envelope_reg_airworth) ...
    ':'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
        % PITCHING ANGLE DIFFERENTIAL EQUATION   
%       
%                      DELTA_CM = (-0.01)*DELTA_AILERON
        %
        myEq = "$ \Delta C_{m} = (-0.01)*\delta_{aileron} ";
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
%sec
str = ['where: '];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------
%       with 
%      
%       DELTA_CM      --> Moment coefficient increment
%       DELTA_AILERON --> Down aileron deflection in degrees in the 
%                         critical condition
%       
%       NOTE: The angle at critical condition DELTA_AILERON must be given
%             in degrees. 
        % Ude
        myEq = "$ \Delta C_{m} = \mathrm{pitching moment coefficient increment;}";
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

        % Lvt
        myEq = "$ \delta_{aileron} = \mathrm{down aileron deflection in degrees at critical condition.}";
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
        
        ol = UnorderedList({ref1, ref2});
%         ol = UnorderedList({ref1,ref2,ref3,...
%             ref4,ref5,ref6, ref7,ref8});
%         ol = UnorderedList({ref1, ref2, ref3,...
%             ref4,ref5,ref6, ref7, ref8, ref9});
        append(subsec,ol);
% ------------------------------------------------------------------------- 
%sec
str = ['The aileron deflection at critical condition must be given' ...
       ' in degrees.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------
        % -----------------------------------------------------------------
        n1         = 1.0; 
        point_A    = "Point A";
        if isfield(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA, 'Unsymmetrical_loads') == 1
            cm_A       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.cm_A.value(1:3:end);
            yA         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.yA.value(1:3:end)';
            yA_unit    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.yA.Attributes.unit;
        elseif isfield(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA, 'Unsymmetrical_loads') == 0
            cm_A       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.cm_A1.value(1:3:end);
            yA         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.yA1.value(1:3:end)';
            yA_unit    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.yA1.Attributes.unit;
        end
%         va         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
%         va_unit    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.Attributes.unit;
%         na         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
%         na_unit    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.Attributes.unit;
%         alfaA      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.value;
%         alpha_unit = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.Attributes.unit;
%         CLA        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value;
%         LA         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA_new.value;
%         L_unit     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA_new.Attributes.unit;
%         LtailA     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.value;

        point_C    = "Point C";
        cm_C       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.cm_C.value(1:3:end);
%         vc         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
%         nc         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value;
%         alfaC      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value;
%         CLC        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value;
%         LC         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.value;
%         LtailC     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.value;
        
        point_D    = "Point D";
        cm_D       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.cm_D.value(1:3:end);
%         vd         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
%         nd         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value;
%         alfaD      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value;
%         CLD        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value;
%         LD         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.value;
%         LtailD     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value;
        
        % -----------------------------------------------------------------
        header = {strcat('y(', yA_unit,')'), 'cmA', 'cmC', 'cmD'};
        %each table row needs of a fieldValue
        %1
        value        = {num2str(yA,4), num2str(cm_A,4), num2str(cm_C,4), ...
                        num2str(cm_D,4)};
%         speeds     = {num2str(va,4); num2str(vc,4); num2str(vd,4); ...
%                       num2str(vg,4); num2str(vf,4); num2str(ve,4)};
%         load_fact  = {num2str(na,4); num2str(nc,4); num2str(nd,4); ...
%                       num2str(ng,4); num2str(nf,4); num2str(ne,4)};
%         alfa       = {num2str(alfaA, 4); num2str(alfaC, 4); num2str(alfaD, 4); ...
%                       num2str(alfaG, 4); num2str(alfaF, 4); num2str(alfaE, 4)};
%         lift_coeff = {num2str(CLA, 4); num2str(CLC, 4); num2str(CLD, 4); ...
%                       num2str(CLG, 4); num2str(CLF, 4); num2str(CLE, 4)};
%         wing_lift  = {num2str(LA, 4); num2str(LC, 4); num2str(LD, 4); ...
%                       num2str(LG, 4); num2str(LF, 4); num2str(LE, 4)};
%         tail_lift  = {num2str(LtailA, 4); num2str(LtailC, 4); num2str(LtailD, 4); ...
%                       num2str(LtailG, 4); num2str(LtailF, 4); num2str(LtailE, 4)};
        fieldValue = [value];
    
          
        tbl = FormalTable(header,fieldValue);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
        tbl = BaseTable(tbl);
        tbl.Title = strcat('Unsymmetrical flight conditions. Pitching moment at various flight conditions.');
        tbl.LinkTarget = 'unsymmtable1';
        add(subsec,tbl);
        % -----------------------------------------------------------------
% -------------------------------------------------------------------------
%cm_A
Check_name = convertCharsToStrings(Aircraft.Certification.Aircraft_Name.value);
if Check_name == "Drone_VLA"
    fig = FormalImage([results_path,'CmComparisonPointA.png']);
    fig.Caption = 'Pithcing moment coefficient - POINT A';
    fig.Height = '5in';
    fig.LinkTarget='cm_unsimm_A';
    add(subsec,fig);
elseif Check_name == "TecnamP92"
    fig = FormalImage([results_path,'CmComparisonPointA1.png']);
    fig.Caption = 'Pithcing moment coefficient - POINT A1';
    fig.Height = '5in';
    fig.LinkTarget='cm_unsimm_A';
    add(subsec,fig);
end

% fig = FormalImage([results_path,'CmComparisonPointA.png']);
% fig.Caption = 'Pithcing moment coefficient - POINT A';
% fig.Height = '5in';
% fig.LinkTarget='cm_unsimm_A';
% add(subsec,fig);

%cm_C
fig = FormalImage([results_path,'CmComparisonPointC.png']);
fig.Caption = 'Pithcing moment coefficient - POINT C';
fig.Height = '5in';
fig.LinkTarget='cm_unsimm_C';
add(subsec,fig);

%cm_D
fig = FormalImage([results_path,'CmComparisonPointD.png']);
fig.Caption = 'Pithcing moment coefficient - POINT D';
fig.Height = '5in';
fig.LinkTarget='cm_unsimm_D';
add(subsec,fig);
add(sec,subsec);

%sub
subsec = Section();
subsec.Title = 'Effect of aileron displacement on the wing torsion';
str = ['According to' ...
    (' ') ...
    char(requirement) ...
    (' ') ...
    char(Global_torsion_full_load_airwort_reg) ...
    ', the aileron displacement cause a significant change of' ...
    ' applied wing torsion. The following diagram show this increment' ...
    ' in red.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------
        % -----------------------------------------------------------------
        n1         = 1.0; 
        point_A          = "Point A";
        if isfield(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA, 'Unsymmetrical_loads') == 1 
            TA_full_airloads = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.TA_full_airloads.value(1:3:end);
            TA_70_airloads   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.TA_70_airloads.value(1:3:end);
            TA_unit          = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.TA_70_airloads.Attributes.unit;
            yA               = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.yA.value(1:3:end)';
            yA_unit          = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.yA.Attributes.unit;
        elseif isfield(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA, 'Unsymmetrical_loads') == 0 
            TA_full_airloads = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.TA1_full_airloads.value(1:3:end);
            TA_70_airloads   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.TA1_70_airloads.value(1:3:end);
            TA_unit          = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.TA1_70_airloads.Attributes.unit;
            yA               = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.yA1.value(1:3:end)';
            yA_unit          = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.yA1.Attributes.unit;
        end
%         va         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
%         va_unit    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.Attributes.unit;
%         na         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
%         na_unit    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.Attributes.unit;
%         alfaA      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.value;
%         alpha_unit = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.Attributes.unit;
%         CLA        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value;
%         LA         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA_new.value;
%         L_unit     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA_new.Attributes.unit;
%         LtailA     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.value;

        point_C    = "Point C";
        TC_full_airloads = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.TC_full_airloads.value(1:3:end);
        TC_70_airloads   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.TC_70_airloads.value(1:3:end);
        TC_unit          = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.TC_70_airloads.Attributes.unit;
%         vc         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
%         nc         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value;
%         alfaC      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value;
%         CLC        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value;
%         LC         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.value;
%         LtailC     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.value;
        
        point_D    = "Point D";
        TD_full_airloads = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.TD_full_airloads.value(1:3:end);
        TD_70_airloads   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.TD_70_airloads.value(1:3:end);
        TD_unit          = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.TD_70_airloads.Attributes.unit;
%         vd         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
%         nd         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value;
%         alfaD      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value;
%         CLD        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value;
%         LD         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.value;
%         LtailD     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value;
        
        % -----------------------------------------------------------------
        header = {strcat('y(', yA_unit,')'), strcat('A_100(', TA_unit,')'), ...
                  strcat('A_70(', TA_unit,')'), strcat('C_100(', TA_unit,')'), ...
                  strcat('C_70(', TA_unit,')'), strcat('D_100(', TA_unit,')'), ...
                  strcat('D_70(', TA_unit,')')};
        %each table row needs of a fieldValue
        %1
        value        = {num2str(yA,4), num2str(TA_full_airloads,4), num2str(TA_70_airloads,4), ...
                        num2str(TC_full_airloads,4), num2str(TC_70_airloads,4), ...
                        num2str(TD_full_airloads,4), num2str(TD_70_airloads,4)};
%         speeds     = {num2str(va,4); num2str(vc,4); num2str(vd,4); ...
%                       num2str(vg,4); num2str(vf,4); num2str(ve,4)};
%         load_fact  = {num2str(na,4); num2str(nc,4); num2str(nd,4); ...
%                       num2str(ng,4); num2str(nf,4); num2str(ne,4)};
%         alfa       = {num2str(alfaA, 4); num2str(alfaC, 4); num2str(alfaD, 4); ...
%                       num2str(alfaG, 4); num2str(alfaF, 4); num2str(alfaE, 4)};
%         lift_coeff = {num2str(CLA, 4); num2str(CLC, 4); num2str(CLD, 4); ...
%                       num2str(CLG, 4); num2str(CLF, 4); num2str(CLE, 4)};
%         wing_lift  = {num2str(LA, 4); num2str(LC, 4); num2str(LD, 4); ...
%                       num2str(LG, 4); num2str(LF, 4); num2str(LE, 4)};
%         tail_lift  = {num2str(LtailA, 4); num2str(LtailC, 4); num2str(LtailD, 4); ...
%                       num2str(LtailG, 4); num2str(LtailF, 4); num2str(LtailE, 4)};
        fieldValue = [value];
    
          
        tbl = FormalTable(header,fieldValue);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
        tbl = BaseTable(tbl);
        tbl.Title = strcat('Unsymmetrical flight conditions. Full and partial torsion at various flight conditions.');
        tbl.LinkTarget = 'unsymmtable2';
        add(subsec,tbl);
        % -----------------------------------------------------------------
% -------------------------------------------------------------------------
%cm_A
Check_name = convertCharsToStrings(Aircraft.Certification.Aircraft_Name.value);
if Check_name == "Drone_VLA"
    fig = FormalImage([results_path,'UnsymmetricalTorsionFullPointA.png']);
    fig.Caption = 'Torsion distribution full loads - POINT A';
    fig.Height = '5in';
    fig.LinkTarget='tors_unsimm_A';
    add(subsec,fig);
elseif Check_name == "TecnamP92"
    fig = FormalImage([results_path,'UnsymmetricalTorsionFullPointA1.png']);
    fig.Caption = 'Torsion distribution full loads - POINT A1';
    fig.Height = '5in';
    fig.LinkTarget='tors_unsimm_A';
    add(subsec,fig);
end

% fig = FormalImage([results_path,'UnsymmetricalTorsionFullPointA.png']);
% fig.Caption = 'Torsion distribution full loads - POINT A';
% fig.Height = '5in';
% fig.LinkTarget='tors_unsimm_A';
% add(subsec,fig);

%cm_C
fig = FormalImage([results_path,'UnsymmetricalTorsionFullPointC.png']);
fig.Caption = 'Torsion distribution full loads - POINT C';
fig.Height = '5in';
fig.LinkTarget='tors_unsimm_C';
add(subsec,fig);

%cm_D
fig = FormalImage([results_path,'UnsymmetricalTorsionFullPointD.png']);
fig.Caption = 'Torsion distribution full loads - POINT D';
fig.Height = '5in';
fig.LinkTarget='torsion_unsimm_D';
add(subsec,fig);
add(sec,subsec);
% -------------------------------------------------------------------------
%         % -----------------------------------------------------------------
%         n1         = 1.0; 
%         point_A    = "Point A";
%         va         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
%         va_unit    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.Attributes.unit;
%         na         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
%         na_unit    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.Attributes.unit;
%         alfaA      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.value;
%         alpha_unit = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alfaA.Attributes.unit;
%         CLA        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value;
%         LA         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA_new.value;
%         L_unit     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA_new.Attributes.unit;
%         LtailA     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.value;
% 
%         point_C    = "Point C";
%         vc         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
%         nc         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value;
%         alfaC      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alfaC.value;
%         CLC        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value;
%         LC         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC_new.value;
%         LtailC     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.value;
%         
%         point_D    = "Point D";
%         vd         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
%         nd         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value;
%         alfaD      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alfaD.value;
%         CLD        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value;
%         LD         = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD_new.value;
%         LtailD     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value;
%         
%         % -----------------------------------------------------------------
%         header = {'Point', strcat('V(',va_unit,')'), strcat('n(',na_unit,')'), ...
%             strcat('alpha(',alpha_unit,')'), strcat('CL'), strcat('L(',L_unit,')'), strcat('L tail(',L_unit,')')};
%         %each table row needs of a fieldValue
%         %1
%         name       = {char(point_A); char(point_C); char(point_D); ...
%                       char(point_G); char(point_F); char(point_E)};
%         speeds     = {num2str(va,4); num2str(vc,4); num2str(vd,4); ...
%                       num2str(vg,4); num2str(vf,4); num2str(ve,4)};
%         load_fact  = {num2str(na,4); num2str(nc,4); num2str(nd,4); ...
%                       num2str(ng,4); num2str(nf,4); num2str(ne,4)};
%         alfa       = {num2str(alfaA, 4); num2str(alfaC, 4); num2str(alfaD, 4); ...
%                       num2str(alfaG, 4); num2str(alfaF, 4); num2str(alfaE, 4)};
%         lift_coeff = {num2str(CLA, 4); num2str(CLC, 4); num2str(CLD, 4); ...
%                       num2str(CLG, 4); num2str(CLF, 4); num2str(CLE, 4)};
%         wing_lift  = {num2str(LA, 4); num2str(LC, 4); num2str(LD, 4); ...
%                       num2str(LG, 4); num2str(LF, 4); num2str(LE, 4)};
%         tail_lift  = {num2str(LtailA, 4); num2str(LtailC, 4); num2str(LtailD, 4); ...
%                       num2str(LtailG, 4); num2str(LtailF, 4); num2str(LtailE, 4)};
%         fieldValue = [name, speeds, load_fact, alfa, lift_coeff, wing_lift, tail_lift];
%     
%           
%         tbl = FormalTable(header,fieldValue);
%         % In order to put a table with a caption, the API Report denomination should
%         % be used, the other options are from API DOM. In order to solve the problem,
%         % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
%         tbl = BaseTable(tbl);
%         tbl.Title = strcat('Unsymmetrical flight conditions.');
%         tbl.LinkTarget = 'unsymmtable';
%         add(sec,tbl);
%         % -----------------------------------------------------------------
% -------------------------------------------------------------------------
add(ch,sec);


%% END chapter
%Adding chapters
add(rpt,ch);