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

add(sec,subsec);

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
fig = FormalImage([results_path,'ShearBendingTorsionDiagramPointA.png']);
 fig.Caption = 'Shear, Bending and Torsion due to airloads - POINT A';
 fig.Height = '4in';
 fig.LinkTarget='A_distribution';
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
 
        % -----------------------------------------------------------------
        n1           = 1.0; 
        point_A      = "Point A";
        va           = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
        va_unit      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.Attributes.unit;
        na           = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
        na_unit      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.Attributes.unit;
        shearA       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_distr.value(end);
        shear_unit   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_distr.Attributes.unit;
        bendingA     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Bend_mom_distr.value(end);
        bending_unit = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Bend_mom_distr.Attributes.unit;
        torsionA     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Tors_mom_distr.value(end);
        torsion_unit = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Tors_mom_distr.Attributes.unit;

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
        header = {'Point', strcat('V (',va_unit,')'), strcat('n (',na_unit,')'), ...
            strcat('S (',shear_unit,')'), strcat('M (',bending_unit,')'), strcat('T (',torsion_unit,')')};
        %each table row needs of a fieldValue
        %1
        name       = {char(point_A); char(point_C); char(point_D); ...
                      char(point_G); char(point_F); char(point_E)};
        speeds     = {num2str(va,4); num2str(vc,4); num2str(vd,4); ...
                      num2str(vg,4); num2str(vf,4); num2str(ve,4)};
        load_fact  = {num2str(na,4); num2str(nc,4); num2str(nd,4); ...
                      num2str(ng,4); num2str(nf,4); num2str(ne,4)};
        shear      = {num2str(shearA, 4); num2str(shearC, 4); num2str(shearD, 4); ... 
                      num2str(shearG, 4); num2str(shearF, 4); num2str(shearE, 4)};
        bending    = {num2str(bendingA, 4); num2str(bendingC, 4); num2str(bendingD, 4); ...
                      num2str(bendingG, 4); num2str(bendingF, 4); num2str(bendingE, 4)};
        torsion    = {num2str(torsionA, 4); num2str(torsionC, 4); num2str(torsionD, 4); ...
                      num2str(torsionG, 4); num2str(torsionF, 4); num2str(torsionE, 4)};
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

% -------------------------------------------------------------------------
%cm_A
fig = FormalImage([results_path,'CmComparisonPointA.png']);
fig.Caption = 'Pithcing moment coefficient - POINT A';
fig.Height = '5in';
fig.LinkTarget='cm_unsimm_a';
add(subsec,fig);

%cm_C
fig = FormalImage([results_path,'CmComparisonPointC.png']);
fig.Caption = 'Pithcing moment coefficient - POINT C';
fig.Height = '5in';
fig.LinkTarget='cm_unsimm_C';
add(subsec,fig);

%cm_D
fig = FormalImage([results_path,'CmComparisonPointC.png']);
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

%cm_A
fig = FormalImage([results_path,'UnsymmetricalTorsionFullPointA.png']);
fig.Caption = 'Torsion distribution full loads - POINT A';
fig.Height = '5in';
fig.LinkTarget='cm_unsimm_a';
add(subsec,fig);

%cm_C
fig = FormalImage([results_path,'UnsymmetricalTorsionFullPointC.png']);
fig.Caption = 'Torsion distribution full loads - POINT C';
fig.Height = '5in';
fig.LinkTarget='cm_unsimm_C';
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