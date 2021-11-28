
%% SUPPLEMENTARY CONDITIONS FOR TAIL SURFACES 
%  ------------------------------------------
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%  CS - VLA 447 COMBINED LOADS ON TAIL SURFACES 
%  
%  (a) With the aeroplane in a loading condition correspondint to point A
%      or point D in the V - n diagram (whichever condition leads to the
%      higher balance load) the loads on the horizontal tail must be
%      combined with those on the vertical tail as specified in CS - VLA
%      441. 
%  (b) 75% of the loads according to CS - VLA 423 for the horizontal tail
%      and CS - VLA 441 for the vertical tail must be assumed acting
%      simultaneously. 
%
%  AMC S 447
%   For aeroplanes where the horizontal stabilising surfaces are arranged
%   considerably above or below the centre of area of the vertical
%   stabilising surfaces, the stabilising surfaces and their supporting
%   structure including the rear portion of the fuselage should be designed
%   to whitstand combined horizontal and vertical loads. In this case, the
%   prescribed loadings on the vertical stabilising surfaces and the roll
%   moments induced at the horizontal stabilising surfaces should be
%   accounted for. 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% SUBPARAGRAPH (a)
% Horizontal tail loads
TailLoadsD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTail_D.value;
TailLoadsA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTail_A.value;

% Choosing design condition 
if abs(TailLoadsD) > abs(TailLoadsA) 
    Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LHTail.value = TailLoadsD;
elseif abs(TailLoadsA) > abs(TailLoadsD) 
    Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LHTail.value = TailLoadsA;
end
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LHTail.Attributes.unit = "daN";

% Vertical tail loads 
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LVTail.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.value;
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LVTail.Attributes.unit = "daN";

% Total loads acting on the tail 
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LTail.value = sqrt( Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LHTail.value^2 + Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LVTail.value^2);
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LTail.Attributes.unit = "daN";
x = Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LHTail.value/Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LVTail.value;
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_rad.value = pi - atan(x);
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_rad.Attributes.unit = "rad"; 
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_rad.value);
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_deg.Attributes.unit = "degrees"; 

% SUBPARAGRAPH (b)
% Horizontal tail loads 
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LHTail.value = 0.75*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.value;
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LHTail.Attributes.unit = "daN";

% Vertical tail loads 
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LVTail.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.value;
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LVTail.Attributes.unit = "daN";

% Total loads acting on the tail 
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LTail.value = sqrt( Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LHTail.value^2 + Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LVTail.value^2);
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LTail.Attributes.unit = "daN";
x = Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LHTail.value/Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LVTail.value;
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_rad.value = pi - atan(x);
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_rad.Attributes.unit = "rad"; 
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_rad.value);
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_deg.Attributes.unit = "degrees"; 

%% CRITICAL COMBINED LOADS 
Tailloads_subpar_a = Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LTail.value;
Tailloads_subpar_b = Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LTail.value;

if abs(Tailloads_subpar_a) > abs(Tailloads_subpar_b) 
    Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.Total_loads.value = Tailloads_subpar_a;
    Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.theta_rad.value = Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_rad.value;
    Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.theta_deg.value = Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_deg.value;
elseif abs(Tailloads_subpar_b) > abs(Tailloads_subpar_a) 
    Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.Total_loads.value = Tailloads_subpar_b;
    Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.theta_rad.value = Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_rad.value;
    Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.theta_deg.value = Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_deg.value;
end
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.Total_loads.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.theta_rad.Attributes.unit = "rad";
Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.theta_deg.Attributes.unit = "degrees";

disp(" ")
% Horizontal tail loads increments
Increment = [Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LTail.value, ...
             Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LTail.value, ...
             Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.Total_loads.value];
disp(" ++++++++++ COMBINED TAIL LOADS - [daN] ++++++++++ ")
format = ' %6.6f           %6.6f           %6.6f\n';
label  = ' L_Tail Subpar. (a) L_Tail Subpar. (b) L_Tail Critical\n';
fprintf(label);
fprintf(format, Increment.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
