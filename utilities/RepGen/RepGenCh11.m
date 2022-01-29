%% CHAPTER 11 - Loads on the htail
%Ref EASA: WING LOAD CALCULATION
%(Example document for LSA applicants – v1 of 08.03.16)
%Date of issue: DD/MM/YYYY
%Document reference: ABCD-FL-57-00

import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

% chapter_number = 11;
% ch = strcat('ch' , num2str(chapter_number)); 
ch = Chapter();
ch.Title = 'Loads on the horizontal tail';
disp(['Chapter 11', (' "'), ch.Title,('" ') ,'writing...' ])

str = ['ADD HERE details '];
para = Paragraph(str);
add(ch,para);

% BALANCING LOADS
requirement            = Aircraft.Certification.Regulation.value;
% CS-VLA 421 Balancing loads
balancing_airworth_req = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Attributes.cs;
flaps_load_req         = char(Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.Attributes.cs);
%sec
sec = Section();
sec.Title = 'Balancing loads';
str = ['According to ' ...
      (' ') ...
      char(requirement) ...
      (' ') ...
      char(balancing_airworth_req) ...
      (' ') ...
      ', a horizontal tail balancing load is a load necessary to maintain' ...
      ' equilibrium in any specified flight condition with no pitching' ...
      ' acceleration. Horizontal tail surfaces must also be designed for'...
      ' the balancing loads occuring at any point on the limit manoeuvring' ...
      ' envelope and in the flap conditions specified in' ...
      (' ') ...
      char(requirement) ...
      (' ') ...
      flaps_load_req(1:4) ...
      '. The distribution in figure B6 of Appendix B may be used. '];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);

fig = FormalImage([results_path,'Balancingloads.png']);
fig.Caption = 'Balancing loads';
fig.Height = '5in';
fig.LinkTarget='bala_loads';
add(sec,fig);

add(ch,sec);

% CS-VLA 423 Manouevring loads
manouevring_load_req = char(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.DeltaLimitLTail.Attributes.cs);
%sec
sec = Section();
sec.Title = 'Manouevring loads';
str = ['According to ' ...
       (' ') ...
       char(requirement) ...
       (' ') ...
       manouevring_load_req(1:4) ...
       ', each horizontal tail surface must be designed for manoeuvring' ...
       ' loads imposed by one of the cited load conditions.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);
% -------------------------------------------------------------------------
%         % -----------------------------------------------------------------
%         %ordered list
%         ref1 = ['a sudden deflection of the elevator control at VA, to' ... 
%                 ' (1) the maximum updward deflection and (2) the maximum' ...
%                 ' downward deflection, as limited by the control stops or' ...
%                 ' pilot effort, whichever is critical. The average loading' ...
%                 ' of B11 of Appendix B and the distribution in figure B7' ...
%                 ' of Appendix B may be used;'];
%         ref2 = ['75 % of the loads according to ' ...
%             (' ') ...
%             char(requirement) ...
%             (' ') ... 
%             char(global_airworth_rules_horiz) ...
%             (' ') ...
%             'for the horizontal tail and' ...
%             (' ') ...
%             char(requirement) ...
%             (' ') ...
%             char(global_airworth_rules) ...
%             (' ') ...
%             'for the vertical tail must be assumed acting simultaneously;' ...
%             ' this prescription results in a combined load equal to' ...
%             (' ') ...
%             strcat(num2str(Tailloads_subpar_b,4)) ...
%             (' ') ...
%             char(Comb_tailloads_unit) ...
%             '.'];
%         ol = OrderedList({ref1, ref2});
%         append(sec,ol);
%         % -----------------------------------------------------------------
% -------------------------------------------------------------------------        

% UNCHECKED MANOEUVRE DATA
pitchup_control_stops   = Aircraft.Certification.Aerodynamic_data.Elevator.Max_deflection.value;
pitchup_control_unit    = Aircraft.Certification.Aerodynamic_data.Elevator.Max_deflection.Attributes.unit;
pitchdown_control_stops = Aircraft.Certification.Aerodynamic_data.Elevator.Max_deflection.value;  

%sub
subsec = Section();
subsec.Title = 'Unchecked manoeuvre';
str = ['At speed VA the pitching control is suddenly displaced to the' ...
    ' maximum deflection as limited by the control stops. The control' ...
    ' stops are' ...
    (' ') ...
    strcat(num2str(pitchup_control_stops,4)) ...
    (' ') ...
    char(pitchup_control_unit) ...
    (' ') ...
    ' pitch up and ' ...
    (' ') ...
    strcat(num2str(pitchdown_control_stops,4)) ...
    (' ') ...
    char(pitchup_control_unit) ...
    (' ') ...
    ' pitch down. Assuming a linear increment of deflection angle,' ...
    ' the tail lift and its moment about the center of gravity pitching' ...
    ' axis grow accordingly. The aircraft angular pitching acceleration' ...
    ' is the consequence, which, at the tail station, leads to a' ...
    ' tangential acceleration nearly normal to tail plane. ' ...
    ' In the time interval delta t a relative speed delta v develops,' ...
    ' which, in composition with the aircraft speed VA causes a decrement' ...
    ' of the tail incidence angle equal to the delta v divided by VA.' ...
    ' This damping effect is the major relevant fact of the control finite' ...
    ' time and its consequence is a less unchecked manoeuvring load.' ...
    ' Taking into account the drag forces, which are opposed to the body' ...
    ' rotation, and other minus occurrencies, a conservative damping' ...
    ' reduction factor of about 0.3 is introduced. This is a standard' ...
    ' assumption for the sudden manoeuvring deflection from neutral' ...
    ' position to stops. Assuming the direction and the intensity of' ...
    ' airspeed at the center of gravity constant during the control time,' ...
    ' the differential equation representing the motion is:'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para)
% -------------------------------------------------------------------------
        % PITCHING ANGLE DIFFERENTIAL EQUATION 
        %
        myEq = "$ \frac{d^{2}\theta}{dt} = \frac{q*S_{tail}*a_{tail}*d}{I_{y}}*\Biggl(\omega*dt - \frac{\Delta v}{V_A}*DF\Biggr) ";
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

        % Ude
        myEq = "$ \theta = \mathrm{rotation pitching angle}";
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
        myEq = "$ q = \mathrm{dynamic pressure (Pa)}";
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
        
        % Kgt
        myEq = "$ S_{tail} = \frac{horizontal tail area (m^2)}";
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
        
%                          2 * M                   K^2
%      mu_gt = ------------------------------- * ------- = lat. mass ratio;
%              rho * c_bar_t * g * a_vt * S_vt   (l_t^2)
%              
        % mu_gt
        myEq = "$ d = \mathrm{C.G.-tail-A.C. distance (m) with } x_{C.G.} = 0.25*MAC ";
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
%              M    = aeroplane mass (kg); 
%              rho  = air density (kg/m^3); 
%              S_vt = area of vertical tail (m^2);
%              l_t  = distance from aeroplane c.g. to lift centre of
%                     vertical surface (m); 
%              a_vt = lift curve slope of vertical tail (1/rad); 
%              V    = aeroplane equivalent speed (m/s); 
%              K    = radius of gyration in yaw (m);
%              g    = acceleration due to gravity (m/s^2).
        % M
        myEq = "$ a_{tail} = \mathrm{tail lift curve slope (1/deg)}; ";
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
        
        % rho
        myEq = "$ I_{y} = \mathrm{airplane pitching inertia moment (kg*m^2)} ";
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
        
        % lt
        myEq = "$ \omega = \mathrm{control angular speed of plane deflection (1/sec)}";
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
        
        % Svt
        myEq = "$ \frac{\Delta v}{V_A} = \mathrm{damping angle} ";
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
        
        % avt
        myEq = "$ DF = \mathrm{damping effect reduction factor} = 0.3 ";
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
        
        ol = UnorderedList({ref1, ref2, ref3, ref4, ...
            ref5, ref6, ref7, ref8, ref9});
%         ol = UnorderedList({ref1,ref2,ref3,...
%             ref4,ref5,ref6, ref7,ref8});
%         ol = UnorderedList({ref1, ref2, ref3,...
%             ref4,ref5,ref6, ref7, ref8, ref9});
        append(subsec,ol);
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------    
%sec
% PITCH UP CASE 423(a) 
DeltaLimitLTail_cs_airworth        = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.DeltaLimitLTail.Attributes.cs;
DeltaLimitLTail_pitch_up           = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.DeltaLimitLTail.value;
DeltaLimitLTail_pitch_up_unit      = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.DeltaLimitLTail.Attributes.unit;
alpha_new_horiz_grad_pitch_up      = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.alpha_new_horiz_grad.value(end);
alpha_new_horiz_grad_pitch_up_unit = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.alpha_new_horiz_grad.Attributes.unit;
% PITCH DOWN CASE 423(a) 
DeltaLimitLTail_pitch_down      = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.DeltaLimitLTail.value;
alpha_new_horiz_grad_pitch_down = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.alpha_new_horiz_grad.value(end);
str = ['It is possible to solve this equation using simple and reliable' ...
    ' numerical methods. According to ' ...
    (' ') ...
    char(requirement) ...
    (' ') ...
    char(DeltaLimitLTail_cs_airworth) ...
    (' ') ...
    ', the following results are presented: '];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------
        %ordered list
        ref1 = ['a pitch up case where the resultant tailplane angle of' ...
            ' attack is ' ...
            (' ') ...
            strcat(num2str(alpha_new_horiz_grad_pitch_up,4)) ...
            (' ') ...
            char(alpha_new_horiz_grad_pitch_up_unit) ...
            (' ') ...
            ' and a corresponding limit tail load of ' ...
            (' ') ...
            strcat(num2str(DeltaLimitLTail_pitch_up,4)) ...
            (' ') ...
            char(DeltaLimitLTail_pitch_up_unit) ...
            ';'];
        ref2 = ['a pitch down case where the resultant tailplane angle of' ...
            ' attack is ' ...
            (' ') ...
            strcat(num2str(alpha_new_horiz_grad_pitch_down,4)) ...
            (' ') ...
            char(alpha_new_horiz_grad_pitch_up_unit) ...
            (' ') ...
            ' and a corresponding limit tail load of ' ...
            (' ') ...
            strcat(num2str(DeltaLimitLTail_pitch_down,4)) ...
            (' ') ...
            char(DeltaLimitLTail_pitch_up_unit) ...
            '.'];
        ol = OrderedList({ref1, ref2});
        append(subsec,ol);
        % -----------------------------------------------------------------
add(sec,subsec);
% -------------------------------------------------------------------------  

%sub

% CS - VLA 423 (b)
Airworth_csvla_methodb = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.critical_tail_airloads.Attributes.cs;
% PITCHING INERTIA MOMENT 
pitching_inertia_moment      = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.IY.value; 
pitching_inertia_moment_unit = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.IY.Attributes.unit;
% METHOD B RESULTS 
Total_critical_loads_methodb      = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.Total_critical_loads.value;
Total_critical_loads_methodb_unit = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.Total_critical_loads.Attributes.unit;

subsec = Section();
subsec.Title = 'Checked manoeuvre';
% -------------------------------------------------------------------------
str = ['According to ' ...
    (' ') ...
    char(requirement) ...
    (' ') ...
    char(Airworth_csvla_methodb) ...
    (' ') ...
    ' a sudden upward deflection of the elevator must be studied,' ...
    ' at speeds abobe VA, followed by a downward deflection of the' ...
    ' elevator, resulting in specified combinations of normal and angular acceleration.' ...
    ' The airplane pitching inertia moment is estimated equal to' ...
    (' ') ...
    strcat(num2str(pitching_inertia_moment,4)) ...
    (' ' ) ...
    char(pitching_inertia_moment_unit) ...
    ' at maximum takeoff weight and center of gravity at 25% of' ...
    ' the mean aerodynamic chord. The maximum limit load in the checked' ...
    ' manoeuvre is ' ...
    (' ') ...
    strcat(num2str(Total_critical_loads_methodb,4)) ...
    (' ' ) ...
    char(Total_critical_loads_methodb_unit) ...
    '.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(subsec,para);
% -------------------------------------------------------------------------
add(sec,subsec);

%sub
subsec = Section();
subsec.Title = 'Gust loads';

add(sec,subsec);


add(ch,sec);

%sec
sec = Section();
sec.Title = 'Horizontal tail loads summary';
str = ['ADD HERE details '];
para = Paragraph(str);
add(ch,para);

add(ch,sec);


%sec
sec = Section();
sec.Title = 'Unsysmmetrical loads';
str = ['ADD HERE details '];
para = Paragraph(str);
add(ch,para);

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