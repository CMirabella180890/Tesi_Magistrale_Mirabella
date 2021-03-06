close all;
%% CalcUnsymmLoads 
% Script to evaluate unsymmetrical load conditions associated with aileron
% deflection.
% =========================================================================
%   DESCRIPTION
%   It is useful to remember the following airworthiness rules: 
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
%   CS-VLA 455 Ailerons 
%   (a) The ailerons must be designed for the loads to which they are
%       subjected  
%       (1) In the neutral position during symmetrical flight conditions; 
%           and  
%       (2) By the following deflections (except as limited by pilot 
%           effort), during unsymmetrical flight conditions; and 
%           (i) Sudden maximum displacement of the aileron control at VA. 
%               Suitable allowance may be made for control system defl-
%               ections. 
%          (ii) Sufficient deflection at VC, where VC is more than VA, to
%               produce a rate of roll not less than obtained in sub-
%               paragraph (a)(2)(i) of this paragraph. 
%         (iii) Sufficient deflection at VD to produce a rate of roll not
%               less than one-third of that obtained in subparagraph 
%               (a)(2)(i) of this paragraph.  
%   (b) The average loading in Appendix B, B11  and figure B1  of Appendix
%       B and the distribution in figure B9 of Appendix B may be used.
%
% =========================================================================
% MAXIMUM AILERON DEFLECTION;
delta_max  = Aircraft.Geometry.Aileron.Max_deflection.value;
y_iniziale = Aircraft.Geometry.Aileron.y_inner.value;
y_finale   = Aircraft.Geometry.Aileron.y_outer.value;

%% STRAIGHT FLIGHT
switch (Straight_flight_Case)
    % CASE 1: VA greater than the intercept
    case 'Case 1'
        % =================================================================
        if max(n_gust_cruise_plus) > nmax
            % =================================================================
            %% AILERON DEFLECTION AT POINT A1
%             if isfield(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope, 'PointA1') == 1
%                 cm_A1      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.cm_A1.value;
%             elseif isfield(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope, 'PointA1') == 0
%                 cm_A1      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cm_A.value;
%             end
            delta_A1   = delta_max;

            % Aileron deflection at Point A1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.delta_A1.value = delta_A1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.delta_A1.Attributes.unit = "deg";

            % Dynamic pressure times delta_A1
            qA1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.value;
            F_Aileron_A1 = qA1 * delta_A1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.F_Aileron_A1.value = F_Aileron_A1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.F_Aileron_A1.Attributes.unit = "N/m^2";

            % Initialization of the Unsymmetrical_load field - Point A
            % half_span = flip(half_span);
            cm_A1_unsymm = zeros(length(half_span), 1);

            for i = 1:length(half_span')
                cm_A1_unsymm(i) = cm_A1(i);
                if (abs(half_span(i) - y_iniziale) < 1e-1) && (abs(half_span(i) - y_finale) > 1e-1)  
                    j = i;
                end
            end
            for i = 1:length(half_span')
                if (abs(half_span(i) - y_finale) < 1e-1) && (abs(half_span(i) - y_iniziale) > 1e-1)
                    k = i;
                end
            end
            for z = j:k 
                Aircraft.Geometry.Aileron.y_span.value(z) = half_span(z);
                Aircraft.Geometry.Aileron.y_span.Attributes.unit = "meters";
                cm_A1_unsymm(z) = cm_A1(z) - 0.01*delta_A1;
            end  

            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.cm_A1.value = cm_A1_unsymm;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.cm_A1.Attributes.unit = "Non dimensional";

            % Pitching moment coefficient diagram 
            disp(" ++++ FIGURE 27 - POINT A1 SYMM. AND UNSYMM. PITCH MOM. COEFFICIENTS ++++ ");
            Pitching_moment_diagram_comparison = Pitching_moment_coefficients_diagram(half_span, cm_A1, cm_A1_unsymm, PointA1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value = Pitching_moment_diagram_comparison;

            exportgraphics(Pitching_moment_diagram_comparison, 'CmComparisonPointA1.pdf', 'ContentType', 'vector')
            exportgraphics(Pitching_moment_diagram_comparison, 'CmComparisonPointA1.png', 'ContentType', 'vector')

            % Saving figures inside correct folder
            fprintf('Saving CmComparisonPointA1.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile CmComparisonPointA1.pdf Output
            movefile CmComparisonPointA1.png Output   

            % Applied torsion associated with maximum aileron deflection at Point A1 on
            % the wing with 100% of the symmetric airloads
            mdistr_full_A1 = zeros(length(cm_A1_unsymm), 1);
            for i = 1:length(cm_A1_unsymm)
                mdistr_full_A1(i) = cm_A1_unsymm(i) * qA1 *((chord_distr(i))^2);
            end
            TA1_full_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_full_A1))*(1e-1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.TA1_full_airloads.value = TA1_full_airloads;
            Aircraft.Certification.Regulation.SubpartC.Final_envelope.PointA1.Unsymmetrical_loads.TA1_full_airloads.Attributes.unit = "daN*m";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.mdistr_full_A1.value = mdistr_full_A1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.mdistr_full_A1.Attributes.unit = "N"; 

            % Applied torsion associated with maximum aileron deflection at Point A1 on
            % the wing with 70% of the symmetric airloads
            mdistr_seventy_A1 = zeros(length(cm_A1_unsymm), 1);
            for i = 1:length(cm_A1_unsymm)
                mdistr_seventy_A1(i) = 0.7*cm_A1_unsymm(i) * qA1 *((chord_distr(i))^2);
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.mdistr_seventy_A1.value = mdistr_seventy_A1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.mdistr_seventy_A1.Attributes.unit = "N";
            TA1_70_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_seventy_A1))*(1e-1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.TA1_70_airloads.value = TA1_70_airloads;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.TA1_70_airloads.Attributes.unit = "daN*m";                

            % Global torsion on the wing        
            Global_torsion_full_load = trapz(half_span, TA1_full_airloads);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_full_load.value = Global_torsion_full_load;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_full_load.Attributes.unit = "daN*m";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_full_load.Attributes.cs = " 349 ";
            Global_torsion_seventypercent_load = trapz(half_span, TA1_70_airloads);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_seventypercent_load.value = Global_torsion_seventypercent_load;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.unit = "daN*m";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.cs = " 349 ";
            
            % Unsymmetrical torsion load due to aileron deflection diagram - 70% AIRLOADS - Point A1 
            disp(" ++++ FIGURE 28 - POINT A1 PARTIAL LOAD - UNSYMM. TORSION ++++ ");
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load = Unsymm_load_diagram(half_span, flip(TA1_70_airloads), flip(Tors_mom_distr_A1), PointA1);

            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointA1.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointA1.png', 'ContentType', 'vector')

            % Saving figures inside correct folder
            fprintf('Saving UnsymmetricalTorsionSeventyPerCentPointA1.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile UnsymmetricalTorsionSeventyPerCentPointA1.pdf Output
            movefile UnsymmetricalTorsionSeventyPerCentPointA1.png Output 

            % Unsymmetrical torsion load due to aileron deflection diagram - FULL AIRLOADS - Point A1 
            disp(" ++++ FIGURE 29 - POINT A1 FULL LOAD - UNSYMM. TORSION ++++ ");
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load = Unsymm_load_diagram(half_span, flip(TA1_full_airloads), flip(Tors_mom_distr_A1), PointA1);

            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointA1.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointA1.png', 'ContentType', 'vector')

            % Saving figures inside correct folder
            fprintf('Saving in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile UnsymmetricalTorsionFullPointA1.pdf Output
            movefile UnsymmetricalTorsionFullPointA1.png Output

            %% AILERON DEFLECTION AT POINT C

            % Aileron deflection at Point C
            VA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
            VC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
            v_ratio = VA / VC;
            delta_C = delta_max * v_ratio;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.delta_C.value = delta_C;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.delta_C.Attributes.unit = "deg";

            % Dynamic pressure times delta_C
            qC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value;
            F_Aileron_C = qC * delta_C;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.F_Aileron_C.value = F_Aileron_C;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.F_Aileron_C.Attributes.unit = "N/m^2";

            % Initialization of the Unsymmetrical_load field - Point C
            cm_C_unsymm = zeros(length(half_span), 1);
            cm_C        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.value;
            for i = 1:length(half_span')
                cm_C_unsymm(i) = cm_C(i);
                if (abs(half_span(i) - y_iniziale) < 1e-1) && (abs(half_span(i) - y_finale) > 1e-1)  
                    j = i;
                end
            end
            for i = 1:length(half_span')
                if (abs(half_span(i) - y_finale) < 1e-1) && (abs(half_span(i) - y_iniziale) > 1e-1)
                    k = i;
                end
            end
            for z = j:k 
                Aircraft.Geometry.Aileron.y_span.value(z) = half_span(z);
                cm_C_unsymm(z) = cm_C(z) - 0.01 * delta_C;
            end        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.cm_C.value = cm_C_unsymm;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.cm_C.Attributes.unit = "Non dimensional";

            % Pitching moment coefficient diagram 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value = Pitching_moment_coefficients_diagram(half_span, cm_C, cm_C_unsymm, PointC);

            disp(" ++++ FIGURE 30 -POINT C SYMM. AND UNSYMM. PITCH MOM. COEFFICIENTS ++++ ");
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value, 'CmComparisonPointC.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value, 'CmComparisonPointC.png', 'ContentType', 'vector')

            % Saving figures inside correct folder
            fprintf('Saving CmComparisonPointC.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile CmComparisonPointC.pdf Output
            movefile CmComparisonPointC.png Output        

            % Applied torsion associated with maximum aileron deflection at Point C on
            % the wing with 100% of the symmetric airloads
            mdistr_full_C = zeros(length(cm_C_unsymm), 1);

            for i = 1:length(cm_C_unsymm)
                mdistr_full_C(i) = cm_C_unsymm(i) * qC *((chord_distr(i))^2);
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.mdistr_full_C.value = mdistr_full_C;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.mdistr_full_C.Attributes.unit = "N";
            TC_full_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_full_C))*(1e-1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.TC_full_airloads.value = TC_full_airloads;
            Aircraft.Certification.Regulation.SubpartC.Final_envelope.PointC.Unsymmetrical_loads.TC_full_airloads.Attributes.unit = "daN*m";

            % Applied torsion associated with maximum aileron deflection at Point C on
            % the wing with 70% of the symmetric airloads
            mdistr_seventy_C = zeros(length(cm_C_unsymm), 1);

            for i = 1:length(cm_C_unsymm)
                mdistr_seventy_C(i) = 0.7 * cm_C_unsymm(i) * qC*((chord_distr(i))^2);
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.mdistr_seventy_C.value = mdistr_seventy_C;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.mdistr_seventy_C.Attributes.unit = "N";
            TC_70_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_seventy_C))*(1e-1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.TC_70_airloads.value = TC_70_airloads;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.TC_70_airloads.Attributes.unit = "daN*m";

            % Global torsion on the wing    
            Global_torsion_full_load = trapz(half_span, TC_full_airloads);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value = Global_torsion_full_load;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.Attributes.unit = "daN*m";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.Attributes.cs = " 349 ";
            Global_torsion_seventypercent_load = trapz(half_span, TC_70_airloads);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_seventypercent_load.value = Global_torsion_seventypercent_load;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.unit = "daN*m";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.cs = " 349 ";

            % Unsymmetrical torsion load due to aileron deflection diagram - 70% AIRLOADS - Point C 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load = Unsymm_load_diagram(half_span, flip(TC_70_airloads), flip(Tors_mom_distr_C), PointC);

            disp(" ++++ FIGURE 31 - POINT C PARTIAL LOAD - UNSYMM. TORSION ++++ ");
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointC.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointC.png', 'ContentType', 'vector')

            % Saving figures inside correct folder
            fprintf('Saving UnsymmetricalTorsionSeventyPerCentPointC.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile UnsymmetricalTorsionSeventyPerCentPointC.pdf Output
            movefile UnsymmetricalTorsionSeventyPerCentPointC.png Output   

            % Unsymmetrical torsion load due to aileron deflection diagram - FULL AIRLOADS - Point C 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load = Unsymm_load_diagram(half_span, flip(TC_full_airloads), flip(Tors_mom_distr_C), PointC);

            disp(" ++++ FIGURE 32 - POINT C FULL LOAD - UNSYMM. TORSION ++++ ");
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointC.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointC.png', 'ContentType', 'vector')

            % Saving figures inside correct folder
            fprintf('Saving UnsymmetricalTorsionSeventyPerCentPointC.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile UnsymmetricalTorsionFullPointC.pdf Output
            movefile UnsymmetricalTorsionFullPointC.png Output     

            %% AILERON DEFLECTION AT POINT D

            % Aileron deflection at Point D
            VA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
            VD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
            v_ratio = (VA) / (3 * VD);
            delta_D = delta_max * v_ratio;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.delta_D.value = delta_D;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.delta_D.Attributes.unit = "deg";       

            % Dynamic pressure times delta_D
            qD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value;
            F_Aileron_D = qD * delta_D;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.F_Aileron_D.value = F_Aileron_D;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.F_Aileron_D.Attributes.unit = "N/m^2";

            % Initialization of the Unsymmetrical_load field - Point D
            cm_D_unsymm = zeros(length(half_span), 1);
            cm_D        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.value;

            for i = 1:length(half_span')
                cm_D_unsymm(i) = cm_D(i);
                if (abs(half_span(i) - y_iniziale) < 1e-2) && (abs(half_span(i) - y_finale) > 1e-2)  
                    j = i;
                end
            end
            for i = 1:length(half_span')
                if (abs(half_span(i) - y_finale) < 1e-2) && (abs(half_span(i) - y_iniziale) > 1e-2)
                    k = i;
                end
            end
            for z = j:k 
                Aircraft.Geometry.Aileron.y_span.value(z) = half_span(z);
                cm_D_unsymm(z) = cm_D(z) - 0.01 * delta_D;
            end    

            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.cm_D.value = cm_D_unsymm;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.cm_D.Attributes.unit = "Non dimensional";

            % Pitching moment coefficient diagram 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value = Pitching_moment_coefficients_diagram(half_span, cm_D, cm_D_unsymm, PointD);

            disp(" ++++ FIGURE 33 -POINT D SYMM. AND UNSYMM. PITCH MOM. COEFFICIENTS ++++ ");
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value, 'CmComparisonPointD.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value, 'CmComparisonPointD.png', 'ContentType', 'vector')

            % Saving figures inside correct folder
            fprintf('Saving CmComparisonPointD.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile CmComparisonPointD.pdf Output
            movefile CmComparisonPointD.png Output  

            % Applied torsion associated with maximum aileron deflection at Point D on
            % the wing with 100% of the symmetric airloads
            mdistr_full_D = zeros(length(cm_D_unsymm), 1);

            for i = 1:length(cm_D_unsymm)
                mdistr_full_D(i) = cm_D_unsymm(i) * qD *((chord_distr(i))^2);
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.mdistr_full_D.value = mdistr_full_D;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.mdistr_full_D.Attributes.unit = "N";
            TD_full_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_full_D))*(1e-1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.TD_full_airloads.value = TD_full_airloads;
            Aircraft.Certification.Regulation.SubpartC.Final_envelope.PointD.Unsymmetrical_loads.TD_full_airloads.Attributes.unit = "daN*m";

            % Applied torsion associated with maximum aileron deflection at Point D on
            % the wing with 70% of the symmetric airloads
            mdistr_seventy_D = zeros(length(cm_D_unsymm), 1);

            for i = 1:length(cm_D_unsymm)
                mdistr_seventy_D(i) = 0.7 * cm_D_unsymm(i) * qD *((chord_distr(i))^2);
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.mdistr_seventy_D.value = mdistr_seventy_D;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.mdistr_seventy_D.Attributes.unit = "N";
            TD_70_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_seventy_D))*(1e-1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.TD_70_airloads.value = TD_70_airloads;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.TD_70_airloads.Attributes.unit = "daN*m";

            % Global torsion on the wing                                                                                                                                      
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value = trapz(half_span, TD_full_airloads);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.Attributes.unit = "daN*m";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.Attributes.cs = " 349 "; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_seventypercent_load.value = trapz(half_span, TD_70_airloads);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.unit = "daN*m";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.cs = " 349 ";

            % Unsymmetrical torsion load due to aileron deflection diagram - 70% AIRLOADS - Point D 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load = Unsymm_load_diagram(half_span, flip(TD_70_airloads), flip(Tors_mom_distr_D), PointD);

            disp(" ++++ FIGURE 34 - POINT D PARTIAL LOAD - UNSYMM. TORSION ++++ ");
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointD.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointD.png', 'ContentType', 'vector')

            % Saving figures inside correct folder
            fprintf('Saving UnsymmetricalTorsionSeventyPerCentPointD.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile UnsymmetricalTorsionSeventyPerCentPointD.pdf Output
            movefile UnsymmetricalTorsionSeventyPerCentPointD.png Output

            % Unsymmetrical torsion load due to aileron deflection diagram - FULL AIRLOADS - Point D 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load = Unsymm_load_diagram(half_span, flip(TD_full_airloads), flip(Tors_mom_distr_D), PointD);

            disp(" ++++ FIGURE 35 - POINT D FULL LOAD - UNSYMM. TORSION ++++ ");
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointD.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointD.png', 'ContentType', 'vector')


            % Saving figures inside correct folder
            fprintf('Saving UnsymmetricalTorsionFullPointD.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile UnsymmetricalTorsionFullPointD.pdf Output
            movefile UnsymmetricalTorsionFullPointD.png Output  

            %% COMPARISON BETWEEN TORSION DISTRIBUTIONS 
            % Full load torsion loads
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value = Compare_unsymm_load(half_span, ...
                                                        flip(TA1_full_airloads), ...
                                                        flip(TC_full_airloads), ...
                                                        flip(TD_full_airloads), ... 
                                                        PointA1, ...
                                                        PointC, ...
                                                        PointD, ...
                                                        'Full load');

            disp(" ++++ FIGURE 36 - FULL LOAD UNSYMM. TORSION COMPARISON ++++ ");
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value, 'UnsymmetricalLoadsComparisonFullLoad.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value, 'UnsymmetricalLoadsComparisonFullLoad.png', 'ContentType', 'vector')


            % Saving figures inside correct folder
            fprintf('Saving UnsymmetricalLoadsComparisonFullLoad.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile UnsymmetricalLoadsComparisonFullLoad.pdf Output
            movefile UnsymmetricalLoadsComparisonFullLoad.png Output

            % 70% torsion loads
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value = Compare_unsymm_load(half_span, ...
                                                        flip(TA1_70_airloads), ...
                                                        flip(TC_70_airloads), ...
                                                        flip(TD_70_airloads), ... 
                                                        PointA1, ...
                                                        PointC, ...
                                                        PointD, ...
                                                        '$70\%$ load');

            disp(" ++++ FIGURE 37 - PARTIAL LOAD UNSYMM. TORSION COMPARISON ++++ ");
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value, 'UnsymmetricalLoadsComparisonSeventyPerCentLoad.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value, 'UnsymmetricalLoadsComparisonSeventyPerCentLoad.png', 'ContentType', 'vector')
            % Saving figures inside correct folder
            fprintf('Saving UnsymmetricalLoadsComparisonSeventyPerCentLoad.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile UnsymmetricalLoadsComparisonSeventyPerCentLoad.pdf Output
            movefile UnsymmetricalLoadsComparisonSeventyPerCentLoad.png Output

            %% CRITICAL CONDITION FOR AILERON 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.value = NaN;
            Condition1 = F_Aileron_C / F_Aileron_A1;
            Condition2 = F_Aileron_D / F_Aileron_A1;
            if Condition1 > Condition2 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.value = Condition1; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.expression = "FC/FA1";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.description = 'Critical condition for aileron at Point C';
            elseif Condition2 > Condition1 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.value = Condition2; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.expression = "FD/FA1";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.description = 'Critical condition for aileron at Point D'; 
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.cs   = " 455 ";

            %% CRITICAL CONDITION FOR TORSION DUE TO AILERON DEFLECTION
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.value = NaN;
            if (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_full_load.value)) && (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value))
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.Attributes.description = 'Critical condition for torsion at Point D';
            elseif (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value)) && (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value))
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_full_load.value; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.Attributes.description = 'Critical condition for torsion at Point A';
            elseif (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_full_load.value)) && (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value))
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.Attributes.description = 'Critical condition for torsion at Point C';
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.Attributes.unit = "N*m";        

        elseif max(n_gust_cruise_plus) < nmax
            % =================================================================  
            %% AILERON DEFLECTION AT POINT A1

            cm_A1      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.cm_A1.value;
            delta_A1   = delta_max;
            y_iniziale = Aircraft.Geometry.Aileron.y_inner.value;
            % Aileron deflection at Point A1
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.delta_A1.value = delta_A1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.delta_A1.Attributes.unit = "deg";

            % Dynamic pressure times delta_A1
            qA1 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.qA1.value;
            F_Aileron_A1 = qA1 * delta_A1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.F_Aileron_A1.value = F_Aileron_A1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.F_Aileron_A1.Attributes.unit = "N/m^2";

            % Initialization of the Unsymmetrical_load field - Point A1
            % half_span = flip(half_span);
            cm_A1_unsymm = zeros(length(half_span), 1);

            for i = 1:length(half_span')
                cm_A1_unsymm(i) = cm_A1(i);
                if (abs(half_span(i) - y_iniziale) < 1e-1) && (abs(half_span(i) - y_finale) > 1e-1)  
                    j = i;
                end
            end
            for i = 1:length(half_span')
                if (abs(half_span(i) - y_finale) < 1e-1) && (abs(half_span(i) - y_iniziale) > 1e-1)
                    k = i;
                end
            end
            for z = j:k 
                Aircraft.Geometry.Aileron.y_span.value(z) = half_span(z);
                Aircraft.Geometry.Aileron.y_span.Attributes.unit = "meters";
                cm_A1_unsymm(z) = cm_A1(z) - 0.01 * delta_A1;
            end  

            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.cm_A1.value = cm_A1_unsymm;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.cm_A1.Attributes.unit = "Non dimensional";

            % Pitching moment coefficient diagram 
            disp(" ++++ FIGURE 27 - POINT A1 SYMM. AND UNSYMM. PITCH MOM. COEFFICIENTS ++++ ");
            Pitching_moment_diagram_comparison = Pitching_moment_coefficients_diagram(half_span, cm_A1, cm_A1_unsymm, PointA1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value = Pitching_moment_diagram_comparison;

            exportgraphics(Pitching_moment_diagram_comparison, 'CmComparisonPointA1.pdf', 'ContentType', 'vector')
            exportgraphics(Pitching_moment_diagram_comparison, 'CmComparisonPointA1.png', 'ContentType', 'vector')

            % Saving figures inside correct folder
            fprintf('Saving CmComparisonPointA1.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile CmComparisonPointA1.pdf Output
            movefile CmComparisonPointA1.png Output   

            % Applied torsion associated with maximum aileron deflection at Point A on
            % the wing with 100% of the symmetric airloads
            mdistr_full_A1 = zeros(length(cm_A1_unsymm), 1);
            for i = 1:length(cm_A1_unsymm)
                mdistr_full_A1(i) = cm_A1_unsymm(i) * qA1 *((chord_distr(i))^2);
            end
            TA1_full_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_full_A1))*(1e-1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.TA1_full_airloads.value = TA1_full_airloads;
            Aircraft.Certification.Regulation.SubpartC.Final_envelope.PointA1.Unsymmetrical_loads.TA1_full_airloads.Attributes.unit = "daN*m";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.mdistr_full_A1.value = mdistr_full_A1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.mdistr_full_A1.Attributes.unit = "N"; 

            % Applied torsion associated with maximum aileron deflection at Point A1 on
            % the wing with 70% of the symmetric airloads
            mdistr_seventy_A1 = zeros(length(cm_A1_unsymm), 1);
            for i = 1:length(cm_A1_unsymm)
                mdistr_seventy_A1(i) = 0.7*cm_A1_unsymm(i) * qA1 *((chord_distr(i))^2);
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.mdistr_seventy_A1.value = mdistr_seventy_A1;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.mdistr_seventy_A1.Attributes.unit = "N";
            TA1_70_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_seventy_A1))*(1e-1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.TA1_70_airloads.value = TA1_70_airloads;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.TA1_70_airloads.Attributes.unit = "daN*m";                

            % Global torsion on the wing        
            Global_torsion_full_load = trapz(half_span, TA1_full_airloads);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_full_load.value = Global_torsion_full_load;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_full_load.Attributes.unit = "daN*m";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_full_load.Attributes.cs = " 349 ";
            Global_torsion_seventypercent_load = trapz(half_span, TA1_70_airloads);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_seventypercent_load.value = Global_torsion_seventypercent_load;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.unit = "daN*m";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_full_load.Attributes.cs = " 349 ";
            
            % Unsymmetrical torsion load due to aileron deflection diagram - 70% AIRLOADS - Point A1 
            disp(" ++++ FIGURE 28 - POINT A1 PARTIAL LOAD - UNSYMM. TORSION ++++ ");
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load = Unsymm_load_diagram(half_span, flip(TA1_70_airloads), flip(Tors_mom_distr_A1), PointA1);

            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointA1.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointA1.png', 'ContentType', 'vector')

            % Saving figures inside correct folder
            fprintf('Saving UnsymmetricalTorsionSeventyPerCentPointA1.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile UnsymmetricalTorsionSeventyPerCentPointA1.pdf Output
            movefile UnsymmetricalTorsionSeventyPerCentPointA1.png Output 

            % Unsymmetrical torsion load due to aileron deflection diagram - FULL AIRLOADS - Point A1 
            disp(" ++++ FIGURE 29 - POINT A1 FULL LOAD - UNSYMM. TORSION ++++ ");
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load = Unsymm_load_diagram(half_span, flip(TA1_full_airloads), flip(Tors_mom_distr_A1), PointA1);

            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointA1.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointA1.png', 'ContentType', 'vector')

            % Saving figures inside correct folder
            fprintf('Saving in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile UnsymmetricalTorsionFullPointA1.pdf Output
            movefile UnsymmetricalTorsionFullPointA1.png Output

            %% AILERON DEFLECTION AT POINT C

            % Aileron deflection at Point C
            VA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
            VC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
            v_ratio = VA / VC;
            delta_C = delta_max * v_ratio;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.delta_C.value = delta_C;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.delta_C.Attributes.unit = "deg";

            % Dynamic pressure times delta_C
            qC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value;
            F_Aileron_C = qC * delta_C;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.F_Aileron_C.value = F_Aileron_C;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.F_Aileron_C.Attributes.unit = "N/m^2";

            % Initialization of the Unsymmetrical_load field - Point C
            cm_C_unsymm = zeros(length(half_span), 1);
            cm_C        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.value;
            for i = 1:length(half_span')
                cm_C_unsymm(i) = cm_C(i);
                if (abs(half_span(i) - y_iniziale) < 1e-1) && (abs(half_span(i) - y_finale) > 1e-1) 
                    j = i;
                end
            end
            for i = 1:length(half_span')
                if (abs(half_span(i) - y_finale) < 1e-2) && (abs(half_span(i) - y_iniziale) > 1e-2)
                    k = i;
                end
            end
            for z = j:k 
                Aircraft.Geometry.Aileron.y_span.value(z) = half_span(z);
                cm_C_unsymm(z) = cm_C(z) - 0.01 * delta_C;
            end        
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.cm_C.value = cm_C_unsymm;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.cm_C.Attributes.unit = "Non dimensional";

            % Pitching moment coefficient diagram 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value = Pitching_moment_coefficients_diagram(half_span, cm_C, cm_C_unsymm, PointC);

            disp(" ++++ FIGURE 30 -POINT C SYMM. AND UNSYMM. PITCH MOM. COEFFICIENTS ++++ ");
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value, 'CmComparisonPointC.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value, 'CmComparisonPointC.png', 'ContentType', 'vector')

            % Saving figures inside correct folder
            fprintf('Saving CmComparisonPointC.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile CmComparisonPointC.pdf Output
            movefile CmComparisonPointC.png Output        

            % Applied torsion associated with maximum aileron deflection at Point C on
            % the wing with 100% of the symmetric airloads
            mdistr_full_C = zeros(length(cm_C_unsymm), 1);

            for i = 1:length(cm_C_unsymm)
                mdistr_full_C(i) = cm_C_unsymm(i) * qC *((chord_distr(i))^2);
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.mdistr_full_C.value = mdistr_full_C;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.mdistr_full_C.Attributes.unit = "N";
            TC_full_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_full_C))*(1e-1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.TC_full_airloads.value = TC_full_airloads;
            Aircraft.Certification.Regulation.SubpartC.Final_envelope.PointC.Unsymmetrical_loads.TC_full_airloads.Attributes.unit = "daN*m";

            % Applied torsion associated with maximum aileron deflection at Point C on
            % the wing with 70% of the symmetric airloads
            mdistr_seventy_C = zeros(length(cm_C_unsymm), 1);

            for i = 1:length(cm_C_unsymm)
                mdistr_seventy_C(i) = 0.7 * cm_C_unsymm(i) * qC*((chord_distr(i))^2);
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.mdistr_seventy_C.value = mdistr_seventy_C;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.mdistr_seventy_C.Attributes.unit = "N";
            TC_70_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_seventy_C))*(1e-1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.TC_70_airloads.value = TC_70_airloads;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.TC_70_airloads.Attributes.unit = "daN*m";

            % Global torsion on the wing    
            Global_torsion_full_load = trapz(half_span, TC_full_airloads);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value = Global_torsion_full_load;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.Attributes.unit = "daN*m";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.Attributes.cs = " 349 ";
            Global_torsion_seventypercent_load = trapz(half_span, TC_70_airloads);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_seventypercent_load.value = Global_torsion_seventypercent_load;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.unit = "daN*m";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.cs = " 349 ";

            % Unsymmetrical torsion load due to aileron deflection diagram - 70% AIRLOADS - Point C 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load = Unsymm_load_diagram(half_span, flip(TC_70_airloads), flip(Tors_mom_distr_C), PointC);

            disp(" ++++ FIGURE 31 - POINT C PARTIAL LOAD - UNSYMM. TORSION ++++ ");
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointC.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointC.png', 'ContentType', 'vector')

            % Saving figures inside correct folder
            fprintf('Saving UnsymmetricalTorsionSeventyPerCentPointC.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile UnsymmetricalTorsionSeventyPerCentPointC.pdf Output
            movefile UnsymmetricalTorsionSeventyPerCentPointC.png Output   

            % Unsymmetrical torsion load due to aileron deflection diagram - FULL AIRLOADS - Point C 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load = Unsymm_load_diagram(half_span, flip(TC_full_airloads), flip(Tors_mom_distr_C), PointC);

            disp(" ++++ FIGURE 32 - POINT C FULL LOAD - UNSYMM. TORSION ++++ ");
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointC.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointC.png', 'ContentType', 'vector')

            % Saving figures inside correct folder
            fprintf('Saving UnsymmetricalTorsionSeventyPerCentPointC.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile UnsymmetricalTorsionFullPointC.pdf Output
            movefile UnsymmetricalTorsionFullPointC.png Output     

            %% AILERON DEFLECTION AT POINT D

            % Aileron deflection at Point D
            VA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
            VD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
            v_ratio = (VA) / (3 * VD);
            delta_D = delta_max * v_ratio;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.delta_D.value = delta_D;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.delta_D.Attributes.unit = "deg";       

            % Dynamic pressure times delta_D
            qD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value;
            F_Aileron_D = qD * delta_D;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.F_Aileron_D.value = F_Aileron_D;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.F_Aileron_D.Attributes.unit = "N/m^2";

            % Initialization of the Unsymmetrical_load field - Point D
            cm_D_unsymm = zeros(length(half_span), 1);
            cm_D        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.value;

            for i = 1:length(half_span')
                cm_D_unsymm(i) = cm_D(i);
                if (abs(half_span(i) - y_iniziale) < 1e-2) && (abs(half_span(i) - y_finale) > 1e-2)  
                    j = i;
                end
            end
            for i = 1:length(half_span')
                if (abs(half_span(i) - y_finale) < 1e-2) && (abs(half_span(i) - y_iniziale) > 1e-2)
                    k = i;
                end
            end
            for z = j:k 
                Aircraft.Geometry.Aileron.y_span.value(z) = half_span(z);
                cm_D_unsymm(z) = cm_D(z) - 0.01 * delta_D;
            end    

            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.cm_D.value = cm_D_unsymm;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.cm_D.Attributes.unit = "Non dimensional";

            % Pitching moment coefficient diagram 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value = Pitching_moment_coefficients_diagram(half_span, cm_D, cm_D_unsymm, PointD);

            disp(" ++++ FIGURE 33 -POINT D SYMM. AND UNSYMM. PITCH MOM. COEFFICIENTS ++++ ");
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value, 'CmComparisonPointD.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value, 'CmComparisonPointD.png', 'ContentType', 'vector')

            % Saving figures inside correct folder
            fprintf('Saving CmComparisonPointD.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile CmComparisonPointD.pdf Output
            movefile CmComparisonPointD.png Output  

            % Applied torsion associated with maximum aileron deflection at Point D on
            % the wing with 100% of the symmetric airloads
            mdistr_full_D = zeros(length(cm_D_unsymm), 1);

            for i = 1:length(cm_D_unsymm)
                mdistr_full_D(i) = cm_D_unsymm(i) * qD *((chord_distr(i))^2);
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.mdistr_full_D.value = mdistr_full_D;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.mdistr_full_D.Attributes.unit = "N";
            TD_full_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_full_D))*(1e-1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.TD_full_airloads.value = TD_full_airloads;
            Aircraft.Certification.Regulation.SubpartC.Final_envelope.PointD.Unsymmetrical_loads.TD_full_airloads.Attributes.unit = "daN*m";

            % Applied torsion associated with maximum aileron deflection at Point D on
            % the wing with 70% of the symmetric airloads
            mdistr_seventy_D = zeros(length(cm_D_unsymm), 1);

            for i = 1:length(cm_D_unsymm)
                mdistr_seventy_D(i) = 0.7 * cm_D_unsymm(i) * qD *((chord_distr(i))^2);
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.mdistr_seventy_D.value = mdistr_seventy_D;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.mdistr_seventy_D.Attributes.unit = "N";
            TD_70_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_seventy_D))*(1e-1);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.TD_70_airloads.value = TD_70_airloads;
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.TD_70_airloads.Attributes.unit = "daN*m";

            % Global torsion on the wing                                                                                                                                      
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value = trapz(half_span, TD_full_airloads);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.Attributes.unit = "daN*m";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.Attributes.cs = " 349 ";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_seventypercent_load.value = trapz(half_span, TD_70_airloads);
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.unit = "daN*m";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.cs = " 349 ";

            % Unsymmetrical torsion load due to aileron deflection diagram - 70% AIRLOADS - Point D 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load = Unsymm_load_diagram(half_span, flip(TD_70_airloads), flip(Tors_mom_distr_D), PointD);

            disp(" ++++ FIGURE 34 - POINT D PARTIAL LOAD - UNSYMM. TORSION ++++ ");
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointD.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointD.png', 'ContentType', 'vector')

            % Saving figures inside correct folder
            fprintf('Saving UnsymmetricalTorsionSeventyPerCentPointD.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile UnsymmetricalTorsionSeventyPerCentPointD.pdf Output
            movefile UnsymmetricalTorsionSeventyPerCentPointD.png Output

            % Unsymmetrical torsion load due to aileron deflection diagram - FULL AIRLOADS - Point D 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load = Unsymm_load_diagram(half_span, flip(TD_full_airloads), flip(Tors_mom_distr_D), PointD);

            disp(" ++++ FIGURE 35 - POINT D FULL LOAD - UNSYMM. TORSION ++++ ");
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointD.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointD.png', 'ContentType', 'vector')


            % Saving figures inside correct folder
            fprintf('Saving UnsymmetricalTorsionFullPointD.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile UnsymmetricalTorsionFullPointD.pdf Output
            movefile UnsymmetricalTorsionFullPointD.png Output  

            %% COMPARISON BETWEEN TORSION DISTRIBUTIONS 
            % Full load torsion loads
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value = Compare_unsymm_load(half_span, ...
                                                        flip(TA1_full_airloads), ...
                                                        flip(TC_full_airloads), ...
                                                        flip(TD_full_airloads), ... 
                                                        PointA1, ...
                                                        PointC, ...
                                                        PointD, ...
                                                        'Full load');

            disp(" ++++ FIGURE 36 - FULL LOAD UNSYMM. TORSION COMPARISON ++++ ");
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value, 'UnsymmetricalLoadsComparisonFullLoad.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value, 'UnsymmetricalLoadsComparisonFullLoad.png', 'ContentType', 'vector')


            % Saving figures inside correct folder
            fprintf('Saving UnsymmetricalLoadsComparisonFullLoad.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile UnsymmetricalLoadsComparisonFullLoad.pdf Output
            movefile UnsymmetricalLoadsComparisonFullLoad.png Output

            % 70% torsion loads
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value = Compare_unsymm_load(half_span, ...
                                                        flip(TA1_70_airloads), ...
                                                        flip(TC_70_airloads), ...
                                                        flip(TD_70_airloads), ... 
                                                        PointA1, ...
                                                        PointC, ...
                                                        PointD, ...
                                                        '$70\%$ load');

            disp(" ++++ FIGURE 37 - PARTIAL LOAD UNSYMM. TORSION COMPARISON ++++ ");
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value, 'UnsymmetricalLoadsComparisonSeventyPerCentLoad.pdf', 'ContentType', 'vector')
            exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value, 'UnsymmetricalLoadsComparisonSeventyPerCentLoad.png', 'ContentType', 'vector')
            % Saving figures inside correct folder
            fprintf('Saving UnsymmetricalLoadsComparisonSeventyPerCentLoad.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile UnsymmetricalLoadsComparisonSeventyPerCentLoad.pdf Output
            movefile UnsymmetricalLoadsComparisonSeventyPerCentLoad.png Output

            %% CRITICAL CONDITION FOR AILERON 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.value = NaN;
            Condition1 = F_Aileron_C / F_Aileron_A1;
            Condition2 = F_Aileron_D / F_Aileron_A1;
            if Condition1 > Condition2 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.value = Condition1; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.expression = "FC/FA1";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.description = 'Critical condition for aileron at Point C';
            elseif Condition2 > Condition1 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.value = Condition2; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.expression = "FD/FA1";
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.description = 'Critical condition for aileron at Point D'; 
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.unit = "Non dimensional";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.cs   = " 455 ";

            %% CRITICAL CONDITION FOR TORSION DUE TO AILERON DEFLECTION
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.value = NaN;
            if (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_full_load.value)) && (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value))
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.Attributes.description = 'Critical condition for torsion at Point D';
            elseif (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value)) && (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value))
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_full_load.value; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.Attributes.description = 'Critical condition for torsion at Point A';
            elseif (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.Unsymmetrical_loads.Global_torsion_full_load.value)) && (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value))
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value; 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.Attributes.description = 'Critical condition for torsion at Point C';
            end
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.Attributes.unit = "N*m";          
        
        end
    % CASE 2: VA lower than the intercept
    case 'Case 2'
        % =================================================================  
        %% AILERON DEFLECTION AT POINT A
        
        cm_A      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cm_A.value;
        delta_A   = delta_max;
        y_iniziale = Aircraft.Geometry.Aileron.y_inner.value;
        % Aileron deflection at Point A
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.delta_A.value = delta_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.delta_A.Attributes.unit = "deg";
        
        % Dynamic pressure times delta_A
        qA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        F_aileron_A = qA * delta_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.F_Aileron_A.value = F_aileron_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.F_Aileron_A.Attributes.unit = "N/m^2";
        
        % Initialization of the Unsymmetrical_load field - Point A
        % half_span = flip(half_span);
        cm_A_unsymm = zeros(length(half_span), 1);

        for i = 1:length(half_span')
            cm_A_unsymm(i) = cm_A(i);
            if (abs(half_span(i) - y_iniziale) < 1e-1) && (abs(half_span(i) - y_finale) > 1e-1)  
                j = i;
            end
        end
        for i = 1:length(half_span')
            if (abs(half_span(i) - y_finale) < 1e-1) && (abs(half_span(i) - y_iniziale) > 1e-1)
                k = i;
            end
        end
        for z = j:k 
            Aircraft.Geometry.Aileron.y_span.value(z) = half_span(z);
            Aircraft.Geometry.Aileron.y_span.Attributes.unit = "meters";
            cm_A_unsymm(z) = cm_A(z) - 0.01*delta_A;
        end  
        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.cm_A.value = cm_A_unsymm;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.cm_A.Attributes.unit = "Non dimensional";
        
        % Pitching moment coefficient diagram 
        disp(" ++++ FIGURE 27 - POINT A SYMM. AND UNSYMM. PITCH MOM. COEFFICIENTS ++++ ");
        Pitching_moment_diagram_comparison = Pitching_moment_coefficients_diagram(half_span, cm_A, cm_A_unsymm, PointA);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value = Pitching_moment_diagram_comparison;
        
        exportgraphics(Pitching_moment_diagram_comparison, 'CmComparisonPointA.pdf', 'ContentType', 'vector')
        exportgraphics(Pitching_moment_diagram_comparison, 'CmComparisonPointA.png', 'ContentType', 'vector')
        
        % Saving figures inside correct folder
        fprintf('Saving CmComparisonPointA.pdf in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile CmComparisonPointA.pdf Output
        movefile CmComparisonPointA.png Output   
        
        % Applied torsion associated with maximum aileron deflection at Point A on
        % the wing with 100% of the symmetric airloads
        mdistr_full_A = zeros(length(cm_A_unsymm), 1);
        for i = 1:length(cm_A_unsymm)
            mdistr_full_A(i) = cm_A_unsymm(i) * qA *((chord_distr(i))^2);
        end
        TA_full_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_full_A))*(1e-1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.TA_full_airloads.value = TA_full_airloads;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.TA_full_airloads.Attributes.unit = "daN*m";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.mdistr_full_A.value = mdistr_full_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.mdistr_full_A.Attributes.unit = "N"; 
        
        % Applied torsion associated with maximum aileron deflection at Point A on
        % the wing with 70% of the symmetric airloads
        mdistr_seventy_A = zeros(length(cm_A_unsymm), 1);
        for i = 1:length(cm_A_unsymm)
            mdistr_seventy_A(i) = 0.7*cm_A_unsymm(i) * qA *((chord_distr(i))^2);
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.mdistr_seventy_A.value = mdistr_seventy_A;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.mdistr_seventy_A.Attributes.unit = "N";
        TA_70_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_seventy_A))*(1e-1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.TA_70_airloads.value = TA_70_airloads;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.TA_70_airloads.Attributes.unit = "daN*m";                
        
        % Global torsion on the wing        
        Global_torsion_full_load = trapz(half_span, TA_full_airloads);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Global_torsion_full_load.value = Global_torsion_full_load;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Global_torsion_full_load.Attributes.unit = "daN*m";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Global_torsion_full_load.Attributes.cs = " 349 ";
        Global_torsion_seventypercent_load = trapz(half_span, TA_70_airloads);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Global_torsion_seventypercent_load.value = Global_torsion_seventypercent_load;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.unit = "daN*m";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.cs = " 349 ";
                
        % Unsymmetrical torsion load due to aileron deflection diagram - 70% AIRLOADS - Point A 
        disp(" ++++ FIGURE 28 - POINT A PARTIAL LOAD - UNSYMM. TORSION ++++ ");
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load = Unsymm_load_diagram(half_span, flip(TA_70_airloads), flip(Tors_mom_distr_A), PointA);
        
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointA.pdf', 'ContentType', 'vector')
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointA.png', 'ContentType', 'vector')
        
        % Saving figures inside correct folder
        fprintf('Saving UnsymmetricalTorsionSeventyPerCentPointA.pdf in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile UnsymmetricalTorsionSeventyPerCentPointA.pdf Output
        movefile UnsymmetricalTorsionSeventyPerCentPointA.png Output 
        
        % Unsymmetrical torsion load due to aileron deflection diagram - FULL AIRLOADS - Point A 
        disp(" ++++ FIGURE 29 - POINT A FULL LOAD - UNSYMM. TORSION ++++ ");
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load = Unsymm_load_diagram(half_span, flip(TA_full_airloads), flip(Tors_mom_distr_A), PointA);
        
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointA.pdf', 'ContentType', 'vector')
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointA.png', 'ContentType', 'vector')
        
        % Saving figures inside correct folder
        fprintf('Saving in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile UnsymmetricalTorsionFullPointA.pdf Output
        movefile UnsymmetricalTorsionFullPointA.png Output

        %% AILERON DEFLECTION AT POINT C
        
        % Aileron deflection at Point C
        VA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
        VC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
        v_ratio = VA / VC;
        delta_C = delta_max * v_ratio;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.delta_C.value = delta_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.delta_C.Attributes.unit = "deg";
        
        % Dynamic pressure times delta_C
        qC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value;
        F_Aileron_C = qC * delta_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.F_Aileron_C.value = F_Aileron_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.F_Aileron_C.Attributes.unit = "N/m^2";
        
        % Initialization of the Unsymmetrical_load field - Point C
        cm_C_unsymm = zeros(length(half_span), 1);
        cm_C        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.value;
        for i = 1:length(half_span')
            cm_C_unsymm(i) = cm_C(i);
            if (abs(half_span(i) - y_iniziale) < 1e-2) && (abs(half_span(i) - y_iniziale) > 1e-2)  
                j = i;
            end
        end
        for i = 1:length(half_span')
            if (abs(half_span(i) - y_iniziale) < 1e-2) && (abs(half_span(i) - y_iniziale) > 1e-2)
                k = i;
            end
        end
        for z = j:k 
            Aircraft.Geometry.Aileron.y_span.value(z) = half_span(z);
            cm_C_unsymm(z) = cm_C(z) - 0.01 * delta_C;
        end        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.cm_C.value = cm_C_unsymm;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.cm_C.Attributes.unit = "Non dimensional";
        
        % Pitching moment coefficient diagram 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value = Pitching_moment_coefficients_diagram(half_span, cm_C, cm_C_unsymm, PointC);
        
        disp(" ++++ FIGURE 30 -POINT C SYMM. AND UNSYMM. PITCH MOM. COEFFICIENTS ++++ ");
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value, 'CmComparisonPointC.pdf', 'ContentType', 'vector')
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value, 'CmComparisonPointC.png', 'ContentType', 'vector')
        
        % Saving figures inside correct folder
        fprintf('Saving CmComparisonPointC.pdf in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile CmComparisonPointC.pdf Output
        movefile CmComparisonPointC.png Output        

        % Applied torsion associated with maximum aileron deflection at Point C on
        % the wing with 100% of the symmetric airloads
        mdistr_full_C = zeros(length(cm_C_unsymm), 1);
         
        for i = 1:length(cm_C_unsymm)
            mdistr_full_C(i) = cm_C_unsymm(i) * qC *((chord_distr(i))^2);
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.mdistr_full_C.value = mdistr_full_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.mdistr_full_C.Attributes.unit = "N";
        TC_full_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_full_C))*(1e-1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.TC_full_airloads.value = TC_full_airloads;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.TC_full_airloads.Attributes.unit = "daN*m";
        
        % Applied torsion associated with maximum aileron deflection at Point C on
        % the wing with 70% of the symmetric airloads
        mdistr_seventy_C = zeros(length(cm_C_unsymm), 1);
        
        for i = 1:length(cm_C_unsymm)
            mdistr_seventy_C(i) = 0.7 * cm_C_unsymm(i) * qC*((chord_distr(i))^2);
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.mdistr_seventy_C.value = mdistr_seventy_C;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.mdistr_seventy_C.Attributes.unit = "N";
        TC_70_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_seventy_C))*(1e-1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.TC_70_airloads.value = TC_70_airloads;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.TC_70_airloads.Attributes.unit = "daN*m";
                
        % Global torsion on the wing    
        Global_torsion_full_load = trapz(half_span, TC_full_airloads);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value = Global_torsion_full_load;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.Attributes.unit = "daN*m";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.Attributes.cs = " 349 ";
        Global_torsion_seventypercent_load = trapz(half_span, TC_70_airloads);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_seventypercent_load.value = Global_torsion_seventypercent_load;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.unit = "daN*m";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.cs = " 349 ";
          
        % Unsymmetrical torsion load due to aileron deflection diagram - 70% AIRLOADS - Point C 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load = Unsymm_load_diagram(half_span, flip(TC_70_airloads), flip(Tors_mom_distr_C), PointC);
        
        disp(" ++++ FIGURE 31 - POINT C PARTIAL LOAD - UNSYMM. TORSION ++++ ");
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointC.pdf', 'ContentType', 'vector')
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointC.png', 'ContentType', 'vector')
        
        % Saving figures inside correct folder
        fprintf('Saving UnsymmetricalTorsionSeventyPerCentPointC.pdf in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile UnsymmetricalTorsionSeventyPerCentPointC.pdf Output
        movefile UnsymmetricalTorsionSeventyPerCentPointC.png Output   
        
        % Unsymmetrical torsion load due to aileron deflection diagram - FULL AIRLOADS - Point C 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load = Unsymm_load_diagram(half_span, flip(TC_full_airloads), flip(Tors_mom_distr_C), PointC);
        
        disp(" ++++ FIGURE 32 - POINT C FULL LOAD - UNSYMM. TORSION ++++ ");
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointC.pdf', 'ContentType', 'vector')
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointC.png', 'ContentType', 'vector')
        
        % Saving figures inside correct folder
        fprintf('Saving UnsymmetricalTorsionSeventyPerCentPointC.pdf in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile UnsymmetricalTorsionFullPointC.pdf Output
        movefile UnsymmetricalTorsionFullPointC.png Output     
        
        %% AILERON DEFLECTION AT POINT D
        
        % Aileron deflection at Point D
        VA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
        VD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
        v_ratio = (VA) / (3 * VD);
        delta_D = delta_max * v_ratio;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.delta_D.value = delta_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.delta_D.Attributes.unit = "deg";       
        
        % Dynamic pressure times delta_D
        qD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value;
        F_Aileron_D = qD * delta_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.F_Aileron_D.value = F_Aileron_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.F_Aileron_D.Attributes.unit = "N/m^2";
        
        % Initialization of the Unsymmetrical_load field - Point D
        cm_D_unsymm = zeros(length(half_span), 1);
        cm_D        = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.value;
        
        for i = 1:length(half_span')
            cm_D_unsymm(i) = cm_D(i);
            if (abs(half_span(i) - y_iniziale) < 1e-2) && (abs(half_span(i) - y_finale) > 1e-2)  
                j = i;
            end
        end
        for i = 1:length(half_span')
            if (abs(half_span(i) - y_finale) < 1e-2) && (abs(half_span(i) - y_iniziale) > 1e-2)
                k = i;
            end
        end
        for z = j:k 
            Aircraft.Geometry.Aileron.y_span.value(z) = half_span(z);
            cm_D_unsymm(z) = cm_D(z) - 0.01 * delta_D;
        end    
        
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.cm_D.value = cm_D_unsymm;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.cm_D.Attributes.unit = "Non dimensional";
        
        % Pitching moment coefficient diagram 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value = Pitching_moment_coefficients_diagram(half_span, cm_D, cm_D_unsymm, PointD);
        
        disp(" ++++ FIGURE 33 -POINT D SYMM. AND UNSYMM. PITCH MOM. COEFFICIENTS ++++ ");
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value, 'CmComparisonPointD.pdf', 'ContentType', 'vector')
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Pitching_moment_diagram_comparison.value, 'CmComparisonPointD.png', 'ContentType', 'vector')
        
        % Saving figures inside correct folder
        fprintf('Saving CmComparisonPointD.pdf in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile CmComparisonPointD.pdf Output
        movefile CmComparisonPointD.png Output  
        
        % Applied torsion associated with maximum aileron deflection at Point D on
        % the wing with 100% of the symmetric airloads
        mdistr_full_D = zeros(length(cm_D_unsymm), 1);
        
        for i = 1:length(cm_D_unsymm)
            mdistr_full_D(i) = cm_D_unsymm(i) * qD *((chord_distr(i))^2);
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.mdistr_full_D.value = mdistr_full_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.mdistr_full_D.Attributes.unit = "N";
        TD_full_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_full_D))*(1e-1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.TD_full_airloads.value = TD_full_airloads;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.TD_full_airloads.Attributes.unit = "daN*m";
        
        % Applied torsion associated with maximum aileron deflection at Point D on
        % the wing with 70% of the symmetric airloads
        mdistr_seventy_D = zeros(length(cm_D_unsymm), 1);
        
        for i = 1:length(cm_D_unsymm)
            mdistr_seventy_D(i) = 0.7 * cm_D_unsymm(i) * qD *((chord_distr(i))^2);
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.mdistr_seventy_D.value = mdistr_seventy_D;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.mdistr_seventy_D.Attributes.unit = "N";
        TD_70_airloads = calc_tors_mom(obj2, half_span, flip(mdistr_seventy_D))*(1e-1);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.TD_70_airloads.value = TD_70_airloads;
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.TD_70_airloads.Attributes.unit = "daN*m";
              
        % Global torsion on the wing                                                                                                                                      
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value = trapz(half_span, TD_full_airloads);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.Attributes.unit = "daN*m";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.Attributes.cs = " 349 ";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_seventypercent_load.value = trapz(half_span, TD_70_airloads);
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.unit = "daN*m";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_seventypercent_load.Attributes.cs = " 349 ";
        
        % Unsymmetrical torsion load due to aileron deflection diagram - 70% AIRLOADS - Point D 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load = Unsymm_load_diagram(half_span, flip(TD_70_airloads), flip(Tors_mom_distr_D), PointD);
                
        disp(" ++++ FIGURE 34 - POINT D PARTIAL LOAD - UNSYMM. TORSION ++++ ");
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointD.pdf', 'ContentType', 'vector')
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Seventy_percent_load, 'UnsymmetricalTorsionSeventyPerCentPointD.png', 'ContentType', 'vector')
        
        % Saving figures inside correct folder
        fprintf('Saving UnsymmetricalTorsionSeventyPerCentPointD.pdf in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile UnsymmetricalTorsionSeventyPerCentPointD.pdf Output
        movefile UnsymmetricalTorsionSeventyPerCentPointD.png Output
        
        % Unsymmetrical torsion load due to aileron deflection diagram - FULL AIRLOADS - Point D 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load = Unsymm_load_diagram(half_span, flip(TD_full_airloads), flip(Tors_mom_distr_D), PointD);
        
        disp(" ++++ FIGURE 35 - POINT D FULL LOAD - UNSYMM. TORSION ++++ ");
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointD.pdf', 'ContentType', 'vector')
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Torsion_due_to_aileron.Full_load, 'UnsymmetricalTorsionFullPointD.png', 'ContentType', 'vector')
        
        
        % Saving figures inside correct folder
        fprintf('Saving UnsymmetricalTorsionFullPointD.pdf in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile UnsymmetricalTorsionFullPointD.pdf Output
        movefile UnsymmetricalTorsionFullPointD.png Output  
        
        %% COMPARISON BETWEEN TORSION DISTRIBUTIONS 
        % Full load torsion loads
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value = Compare_unsymm_load(half_span, ...
                                                    flip(TA_full_airloads), ...
                                                    flip(TC_full_airloads), ...
                                                    flip(TD_full_airloads), ... 
                                                    PointA, ...
                                                    PointC, ...
                                                    PointD, ...
                                                    'Full load');
        
        disp(" ++++ FIGURE 36 - FULL LOAD UNSYMM. TORSION COMPARISON ++++ ");
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value, 'UnsymmetricalLoadsComparisonFullLoad.pdf', 'ContentType', 'vector')
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value, 'UnsymmetricalLoadsComparisonFullLoad.png', 'ContentType', 'vector')
        
        
        % Saving figures inside correct folder
        fprintf('Saving UnsymmetricalLoadsComparisonFullLoad.pdf in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile UnsymmetricalLoadsComparisonFullLoad.pdf Output
        movefile UnsymmetricalLoadsComparisonFullLoad.png Output
        
        % 70% torsion loads
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value = Compare_unsymm_load(half_span, ...
                                                    flip(TA_70_airloads), ...
                                                    flip(TC_70_airloads), ...
                                                    flip(TD_70_airloads), ... 
                                                    PointA, ...
                                                    PointC, ...
                                                    PointD, ...
                                                    '$70\%$ load');
        
        disp(" ++++ FIGURE 37 - PARTIAL LOAD UNSYMM. TORSION COMPARISON ++++ ");
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value, 'UnsymmetricalLoadsComparisonSeventyPerCentLoad.pdf', 'ContentType', 'vector')
        exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Unsymmetrical_loads_comparison.value, 'UnsymmetricalLoadsComparisonSeventyPerCentLoad.png', 'ContentType', 'vector')
        % Saving figures inside correct folder
        fprintf('Saving UnsymmetricalLoadsComparisonSeventyPerCentLoad.pdf in: ');
        fprintf('\n'); 
        fprintf('%s\n', SaveFolder);
        % Moving file inside correct folder
        movefile UnsymmetricalLoadsComparisonSeventyPerCentLoad.pdf Output
        movefile UnsymmetricalLoadsComparisonSeventyPerCentLoad.png Output
        
        %% CRITICAL CONDITION FOR AILERON 
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.value = NaN;
        Condition1 = F_Aileron_C / F_aileron_A;
        Condition2 = F_Aileron_D / F_aileron_A;
        if Condition1 > Condition2 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.value = Condition1; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.expression = "FC/FA";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.description = 'Critical condition for aileron at Point C';
        elseif Condition2 > Condition1 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.value = Condition2; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.expression = "FD/FA";
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.description = 'Critical condition for aileron at Point D'; 
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.unit = "Non dimensional";
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_condition.Attributes.cs   = " 455 ";
        
        %% CRITICAL CONDITION FOR TORSION DUE TO AILERON DEFLECTION
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.value = NaN;
        if (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Global_torsion_full_load.value)) && (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value))
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.Attributes.description = 'Critical condition for torsion at Point D';
        elseif (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value)) && (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value))
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Global_torsion_full_load.value; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.Attributes.description = 'Critical condition for torsion at Point A';
        elseif (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Unsymmetrical_loads.Global_torsion_full_load.value)) && (abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value) > abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Unsymmetrical_loads.Global_torsion_full_load.value))
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Unsymmetrical_loads.Global_torsion_full_load.value; 
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.Attributes.description = 'Critical condition for torsion at Point C';
        end
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Aileron_critical_torsion.Attributes.unit = "N*m";          
        
end

