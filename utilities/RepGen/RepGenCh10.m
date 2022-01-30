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
str = ['In this section will be shown all the resulting internal' ...
    ' forces acting on the wing structural elements; having calculated' ...
    ' lift, drag and pitching moment coefficient distribution on the wing' ...
    ' with a panel method and the geometrical chord distribution, it is' ...
    ' possible to evaluate normal and shear forces and pitching moment' ...
    ' distributions along the wing span.'];
para = Paragraph(str);
add(ch,para);


%sec
sec = Section();
sec.Title = 'Influence of the fuselage';
str = ['The effects of the fuselage on the wing span lift distribution' ...
    ' cause a reduction of lift at stations near the wing root; this lift' ...
    ' reduction can be discounted because is often negligible, leading to' ...
    ' a more conservative design loads. On the other hand, its influence on' ...
    ' the aeroplane equilibrium is accounted for, in particular on the' ...
    ' pitching moment distribution.'];
% On the other hand, its influence on the aeroplane equilibrium is accounted for 
% as a contribution to the pitching moment. In particular the following formula
% (taken from ref [9]) provides the shift of the aerodynamic centre of the wing
% due to fuselage pitching moment.

para = Paragraph(str);
add(sec,para);

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
        ref2 = ['the  loads  resulting  from  the  aileron  deflections' ...
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