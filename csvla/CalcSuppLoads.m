
%% SUPPLEMENTARY CONDITIONS FOR TAIL SURFACES 
%  ------------------------------------------
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%  CS - VLA 447 COMBINED LOADS ON TAIL SURFACES 
%  
%  (a) With the aeroplane in a loading condition corresponding to point A
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

switch (Straight_flight_Case)
    % CASE 1: VA greater than the intercept
    case 'Case 1'

        disp(" ")
        disp(" ++++ CS - VLA 447 COMBINED LOADS ON TAIL SURFACES ++++ ")
        % SUBPARAGRAPH (a)
        % Horizontal tail loads
        TailLoadsD  = 0.5*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value;
        TailLoadsA1 = 0.5*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.value;

        % Choosing design condition 
        if abs(TailLoadsD) > abs(TailLoadsA1) 
            Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LHTail.value = TailLoadsD;
        elseif abs(TailLoadsA1) > abs(TailLoadsD) 
            Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LHTail.value = TailLoadsA1;
        end
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LHTail.Attributes.unit = "daN";

        % Vertical tail loads 
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LVTail.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.value;
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LVTail.Attributes.unit = "daN";

        % Total loads acting on the tail 
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LTail.value = sqrt( Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LHTail.value^2 + Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LVTail.value^2);
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LTail.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LTail.Attributes.cs = " 447(a) ";
        x = abs(Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LHTail.value)/abs(Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LVTail.value);
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_rad.value = 2*pi - atan(x);
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_rad.Attributes.unit = "rad"; 
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_rad.value);
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_deg.Attributes.unit = "deg"; 

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
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LTail.Attributes.cs = " 447(b) ";
        x = abs(Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LHTail.value)/abs(Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LVTail.value);
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_rad.value = 2*pi - atan(x);
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_rad.Attributes.unit = "rad"; 
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_rad.value);
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_deg.Attributes.unit = "deg"; 

        %% CRITICAL COMBINED LOADS 
        Tailloads_subpar_a = Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LTail.value;
        Tailloads_subpar_b = Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LTail.value;

        if abs(Tailloads_subpar_a) > abs(Tailloads_subpar_b) 
            Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.Total_loads.value = Tailloads_subpar_a;
            Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.Total_loads.Attributes.cs = " 447(a) ";
            Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.theta_rad.value = Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_rad.value;
            Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.theta_deg.value = Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_deg.value;
        elseif abs(Tailloads_subpar_b) > abs(Tailloads_subpar_a) 
            Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.Total_loads.value = Tailloads_subpar_b;
            Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.Total_loads.Attributes.cs = " 447(b) ";
            Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.theta_rad.value = Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_rad.value;
            Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.theta_deg.value = Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_deg.value;
        end
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.Total_loads.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.theta_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.theta_deg.Attributes.unit = "deg";

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

    % CASE 1: VA lower than the intercept
    case 'Case 2'
        
        disp(" ")
        disp(" ++++ CS - VLA 447 COMBINED LOADS ON TAIL SURFACES ++++ ")
        % SUBPARAGRAPH (a)
        % Horizontal tail loads
        TailLoadsD = 0.5*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value;
        TailLoadsA = 0.5*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.value;

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
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LTail.Attributes.cs = " 447(a) ";
        x = abs(Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LHTail.value)/abs(Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.LVTail.value);
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_rad.value = 2*pi - atan(x);
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_rad.Attributes.unit = "rad"; 
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_rad.value);
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_a.theta_deg.Attributes.unit = "deg"; 

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
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LTail.Attributes.cs = " 447(b) ";
        x = abs(Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LHTail.value)/abs(Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.LVTail.value);
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_rad.value = 2*pi - atan(x);
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_rad.Attributes.unit = "rad"; 
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_rad.value);
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.subparagraph_b.theta_deg.Attributes.unit = "deg"; 

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
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.Total_loads.Attributes.cs = " 447 ";
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.theta_rad.Attributes.unit = "rad";
        Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.theta_deg.Attributes.unit = "deg";

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
        
end
        
%% GRAPHICAL REPRESENTATION OF THE SUPPLEMENTARY LOAD CONDITIONS 

Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.Graphical_representation.value = figure;
hold on; grid on; grid minor; 

% SWITCH CASE TO TAKE INTO ACCOUNT THE EMPENNAGE CONFIGURATION 
switch(Aircraft.Geometry.Vertical.empennage_flag.value)
    case 'Double fin'
        combined_load_figure = figure(160);
        hold on ;
        grid on; grid minor;
        plot([0 Aircraft.Geometry.Horizontal.b.value], [0 0], '-k.', 'LineWidth', 2.5, 'DisplayName', 'Horiz. tail');
        plot([Aircraft.Geometry.Horizontal.b.value Aircraft.Geometry.Horizontal.b.value], [0 Aircraft.Geometry.Vertical.b.value], '-k.', 'LineWidth', 2.5, 'DisplayName', 'Vert. tail');
        critical_load = Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.Total_loads.value;
        theta_rad     = Aircraft.Certification.Regulation.SubpartC.CombinedLoads.Critical_condition.theta_rad.value;
        HTail_comp    = critical_load*(sin(theta_rad))*0.0125;
        VTail_comp    = critical_load*(cos(theta_rad))*0.0125;
        X_vertical = [Aircraft.Geometry.Horizontal.b.value Aircraft.Geometry.Horizontal.b.value+VTail_comp];
        Y_vertical = [Aircraft.Geometry.Vertical.b.value*0.75 Aircraft.Geometry.Vertical.b.value*0.75];
        plot(X_vertical, Y_vertical, '-g.', 'MarkerSize', 15, 'DisplayName', 'Vert. loads')
        
        X_horizontal = [Aircraft.Geometry.Horizontal.b.value*0.75 Aircraft.Geometry.Horizontal.b.value*0.75];
        Y_horizontal = [0.0 HTail_comp];
        plot(X_horizontal, Y_horizontal, '-b.', 'MarkerSize', 15, 'DisplayName', 'Horiz. loads')
        
        plot([Aircraft.Geometry.Horizontal.b.value*0.75 Aircraft.Geometry.Horizontal.b.value+VTail_comp], ...
             [HTail_comp HTail_comp], '--b.', 'LineWidth', 0.5, 'MarkerSize', 15, 'DisplayName', 'Projection')
        plot([Aircraft.Geometry.Horizontal.b.value+VTail_comp Aircraft.Geometry.Horizontal.b.value+VTail_comp], ...
             [Aircraft.Geometry.Vertical.b.value*0.75 HTail_comp], '--g.', 'LineWidth', 0.5, 'MarkerSize', 15, 'DisplayName', 'Projection') 
        plot([Aircraft.Geometry.Horizontal.b.value*0.75 Aircraft.Geometry.Horizontal.b.value+VTail_comp], ...
             [Aircraft.Geometry.Vertical.b.value*0.75 HTail_comp], '-r.', 'MarkerSize', 15, 'DisplayName', 'Crit. loads')
         
        xlabel("$Y_b$ - $(m)$", "Interpreter", "latex")
        ylabel("$Z_b$ - $(m)$", "Interpreter", "latex")
        title("Combined loads on the empennage", "Interpreter", "latex")
        legend('Interpreter', 'latex', 'Location', 'northwestoutside')
        % legend({'Horiz. empennage', 'Vertical empennage', 'H load', 'V load', 'Critical load'}, 'Interpreter', 'latex', 'Location', 'northwestoutside')
        
        xlim([0 (Aircraft.Geometry.Horizontal.b.value+0.5)])
        ylim([-(Aircraft.Geometry.Vertical.b.value+0.5) Aircraft.Geometry.Vertical.b.value+0.25]);

        % EXPORT FIGURE
        exportgraphics(combined_load_figure, 'Combinedload.pdf', 'ContentType', 'vector')
        exportgraphics(combined_load_figure, 'Combinedload.png', 'ContentType', 'vector')

        % Saving figures inside correct folder
        fprintf('Saving Combinedload.pdf in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile Combinedload.pdf Output
        movefile Combinedload.png Output 
        % -----------------------------------------------------------------
    case 'Single fin'
        
    case 'T tail'
        
    case 'Others'
end