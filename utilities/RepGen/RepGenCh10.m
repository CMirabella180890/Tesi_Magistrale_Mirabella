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

ch = Chapter();
ch.Title = 'Loads on the wing';
disp(['Chapter 10', (' "'), ch.Title,('" ') ,'writing...' ])


str = ['ADD HERE details for balancing Equation'];
para = Paragraph(str);
add(ch,para);


%sec
sec = Section();
sec.Title = 'Influence of the fuselage';
str = ['ADD HERE details for fuselage effect how are they accounted?'];
para = Paragraph(str);
add(ch,para);

add(ch,sec);


%sec
sec = Section();
sec.Title = 'Forces and moments acting on the wings';
str = ['ADD HERE details for balancing Equation'];
para = Paragraph(str);
add(ch,para);
%sub
subsec = Section();
subsec.Title = 'SpanWise Airloads Distribution';

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

%sub
subsec = Section();
subsec.Title = 'Normal and parallel component';

add(sec,subsec);

%sub
subsec = Section();
subsec.Title = 'Shear, Bending and Torsion';

%A
fig = FormalImage([results_path,'ShearBendingTorsionDiagramPointA.png']);
 fig.Caption = 'Shear, Bending and Torsion due to airloads - POINT A';
 fig.Height = '5in';
 fig.LinkTarget='A_distribution';
 add(subsec,fig);

 %C
fig = FormalImage([results_path,'ShearBendingTorsionDiagramPointC.png']);
 fig.Caption = 'Shear, Bending and Torsion due to airloads - POINT C';
 fig.Height = '5in';
 fig.LinkTarget='C_distribution';
 add(subsec,fig);

  %D
fig = FormalImage([results_path,'ShearBendingTorsionDiagramPointD.png']);
 fig.Caption = 'Shear, Bending and Torsion due to airloads - POINT D';
 fig.Height = '5in';
 fig.LinkTarget='D_distribution';
 add(subsec,fig);
 
add(sec,subsec);

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
add(ch,sec);

%sec
sec = Section();
sec.Title = 'Unsymmetrical loads';
str = ['ADD HERE details for uns loads'];
para = Paragraph(str);
add(ch,para);

%sub
subsec = Section();
subsec.Title = 'Rolling condition';

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


add(ch,sec);


%% END chapter
%Adding chapters
add(rpt,ch);