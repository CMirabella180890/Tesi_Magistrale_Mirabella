%% VERTICAL TAIL SURFACES 
%  Digitalization of the airworthiness rules for the vertical tail loads
%  calculation and structural sizing. 

disp(" ");
disp(" ++++ STARTING VERTICAL TAIL LOADS CALCULATION ++++ ");
disp(" ");
disp(" ---------------------------------------------------------------- ");
disp(" ++++ CS - VLA 441 - VERTICAL TAIL LOADS - MANOEUVRING LOADS ++++ ");
disp(" ---------------------------------------------------------------- ");

switch (Aircraft.Geometry.Vertical.empennage_flag.value)
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    % SINGLE FIN
    case 'Single fin'
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %% CS - VLA 441 MANOEUVRING LOADS 
        %
        %  (a) At speeds up to VA, the vertical tail surfaces must be designed to
        %      withistand the following conditions.In computing the tail loads, the
        %      yawing velocity may be assumed to be zero
        %      (1) whit the aeroplane in unaccelerated flight at zero yaw, it is
        %          assumed that the rudder control is suddenly displaced to the
        %          maximum deflection, as limited by the control stops or by limit
        %          pilot forces; 
        %      (2) with the rudder deflected as specified in sub-paragraph (a)(1)
        %          of this paragraph, it is assumed that the aeroplane yaws to the
        %          resulting sideslip angle. In lieu of a rational analysis, an
        %          overswing angle equal to 1.3 times the static sideslip angle of
        %          sub-paragraph (a)(3) of this paragraph may be assumed; 
        %      (3) A yaw angle of 15.0 degrees with the rudder control maintained
        %          in the neutral position, except as limited by pilot strength. 
        %
        %  (b) The average loading of Appendix B, B11 and figure B1 of Appendix B 
        %      and the distribution of figures B6, B7, B8 of Appendix B may be used
        %      instead of requirements of sub-paragraphs (a)(2), (a)(1) and (a)(3)
        %      of this paragraph, respectively.
        %
        %  (c) The yaw angles specified in sub-paragraph (a)(3) of this paragraph
        %      may be reduced if the yaw angle chosen for a particular speed 
        %      cannot be exceeded in 
        %      (1) steady sideslip conditions; and 
        %      (2) uncoordinated rolls from steep banks.

        %% AMC VLA 441 MANOEUVRING LOADS 
        %  -----------------------------
        %  INTERPRETATIVE MATERIAL AND ACCEPTABLE MEANS OF COMPLIANCE 
        % 
        %  For aeroplanes where the horizontal tail is supported by the vertical
        %  tail, the tail surfaces and their supporting structure, including the
        %  rear portion of the fuselage, should be designed to whitstand the
        %  prescribed loading on the vertica ltail and the roll-moments induced by
        %  the horizontal tail acting in the same direction.

        %  For T - Tails, in the abscence of a more rational analysis, the rolling
        %  moment induced by deflection of the vertical rudder may be computed as
        %  follows: 
        %                            rho0
        %  M_rudder = (0.3) * S_ht * ---- * beta * V^2 * b_ht 
        %                              2
        %  where 
        %
        %  M_rudder = induced roll - moment at horizontal tail (N * m) 
        %  b_ht     = span of horizonta tail (m) 
        %  V        = airspeed (m/s)
        %  S_ht     = area of horizontal tail (m^2)
        %  rho0     = air densidity at sea level (kg/m^3)
        %  beta     = angle of zero - lift line due to rudder deflection; this
        %             angle can be defined as follow         
        %                         dL
        %                       ------- * eta * f_n
        %                        d eta 
        %             where 
        %               dL
        %             ------- = change of zero - lift angle of eta * f_n = 1
        %              d eta         
        %       
        %             eta     = rudder deflection 
        %             f_n     = effectivity factor in accordance with angle of
        %                       rudder deflection

        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.global_requirement.Attributes.cs  = " 441 ";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.global_requirement.Attributes.amc = " AMC 441 "; 

        % APPLICABLE YAW ANGLE 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.value = 15.0;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.Attributes.unit = "deg";

        % APPLICABLE LOAD FACTOR 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.load_factor_n.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.load_factor_n.Attributes.unit = "g's";
        % WING LOADING IN PASCALS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.value = Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.load_factor_n.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.Attributes.unit = "Pa"; 
        % AREA OF THE VERTICAL TAIL 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value = Aircraft.Geometry.Vertical.S.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.Attributes.unit = "m^2";
        % VERTICAL TAIL SURFACE RATIO: S_vt / S_wing
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value/Aircraft.Geometry.Wing.S.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.Attributes.unit = " Non dimensional";
        % PRODUCT OF WING LOADING TIMES THE VERTICAL TAIL SURFACES 
        %
        % P = WS * S * (1/g) --> [P] = [kg]
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Reference_mass.value = (1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Reference_mass.Attributes.unit = "kg";
        % VERTICAL TAIL SPAN b_vt
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.b_vt.value = Aircraft.Geometry.Vertical.b.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.b_vt.Attributes.unit = "m";
        % CALCULATING THE RATIO OF WING LOADING IN KILOGRAMS PER UNIT SPAN 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Spanwise_Load_ratio.value = (1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.b_vt.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Spanwise_Load_ratio.Attributes.unit = "kg/m";
        % VERTICAL TAIL ROOT CHORD 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.value = Aircraft.Geometry.Vertical.croot.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.Attributes.unit = "m";
        % QUARTER CHORD DISTANCE: C_ROOT_VT / 4
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_quarter.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.value/4;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_quarter.Attributes.unit = "m";
        % CHORDWISE KILOGRAMS OF LOAD DISTRIBUTION 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Chordwise_Load_ratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Reference_mass.value/Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Chordwise_Load_ratio.Attributes.unit = "kg/m";

        %% CS - VLA 441 CASE (a)(1)
        % A VECTOR CALLED DR
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.dr.value = [0.0 20.0]';
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.dr.Attributes.unit = "Unknown";

        % LATERAL FORCE COEFFICIENT AT TWO BETA - FROM CFD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_vector.value = [0.00003615686 0.01291153]';
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_vector.Attributes.unit = "Non dimensional";

        % LATERAL FORCE COEFFICIENT WHEN BETA = 0 AND DR = 0
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_0.value = 0.000036156860;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_0.Attributes.unit = "Non dimensional";

        % LATERAL FORCE GRADIENT C_Y / dr
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CYdr.value = 0.000644;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CYdr.Attributes.unit = "1/deg";

        % MAXIMUM RUDDER DEFLECTION 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.value = 30.0;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.Attributes.unit = "deg";

        % LATERAL FORCE COEFFICIENT WHEN THE RUDDER IS AT MAXIMUM DEFLECTION:
        %              d CY
        % CY = CY_0 + ------- * max(delta_rudder)
        %               d r
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_0.value + Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CYdr.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY.Attributes.unit = "Non dimensional";

        % CALCULATION OF: CY / S_ratio
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY.value/Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.Attributes.unit = "Non dimensional";

        % CALCULATION OF LATERAL FORCE

        % DYNAMIC PRESSURE AT VA
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.Attributes.unit = "Pa";

        % LATERAL FORCE 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.Attributes.unit = "N";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.Attributes.cs = " 441(a)(1) ";

        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_decanewton.value = (1e-1)*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_decanewton.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_decanewton.Attributes.cs = " 441(a)(1) ";

        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_kilograms.value = (1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_kilograms.Attributes.unit = "kg";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_kilograms.Attributes.cs = " 441(a)(1) ";

%         % THE TOTAL LATERAL FORCE MUST BE DIVIDED BY TWO, TO HAVE THE AIRLOADS
%         % ACTING ON A SINGLE VERTICAL FIN
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.value*0.5;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin.Attributes.unit = "N"; 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin.Attributes.cs = " 441(a)(1) ";
% 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.value = (1e-1)*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.value*0.5;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.Attributes.unit = "daN"; 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.Attributes.cs = " 441(a)(1) ";
% 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_kilograms.value = (1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.value*0.5;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_kilograms.Attributes.unit = "kg"; 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_kilograms.Attributes.cs = " 441(a)(1) ";

        % DISPLAYING RESULTS 
        disp( " ")
        disp(" ++++ CS - VLA 441 - CASE (a)(1) ++++ ");
        disp(" ------------------------------------ ");

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_decanewton.value];
        % Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.value];
        disp(" ++++++++++ Vertical tail loads - CASE (a)(1) [daN] ++++++++++ ")
        format = ' %6.6f\n';
        label  = ' Y\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

        %% CS - VLA 441 CASE (a)(2)

        % DYNAMIC PRESSURE AT VA
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.Attributes.unit = "Pa";

        % CALCULATION OF (1.3)*(YAW_ANGLE)
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.yaw_angle.value = (1.3)*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.yaw_angle.Attributes.unit = "deg";

        % LATERAL FORCE COEFFICIENT 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_VTP.value = 0.0245;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_VTP.Attributes.unit = "Non dimensional";

        % LATERAL FORCE COEFFICIENT 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_VTP.value/Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.Attributes.unit = "Non dimensional";

        % LATERAL FORCE - NEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.Attributes.unit = "N";

        % LATERAL FORCE - DECANEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.value = (1e-1)*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.Attributes.unit = "daN";

        % LATERAL FORCE - KILOGRAMS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.value = (1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.Attributes.unit = "N";

%         % -------------------------------------------------------------------------
%         % LATERAL FORCE ON A SINGLE FIN - NEWTON
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.value = (0.5)*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.Attributes.unit = "N";
         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.Attributes.cs = " 441(a)(2) ";
% 
%         % LATERAL FORCE ON A SINGLE FIN - DECANEWTON
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.value = (0.5)*(1e-1)*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.Attributes.unit = "daN";
         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.Attributes.cs = " 441(a)(2) ";
% 
%         % LATERAL FORCE ON A SINGLE FIN - KILOGRAMS
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.value = (0.5)*(1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.Attributes.unit = "N";
         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.Attributes.cs = " 441(a)(2) ";

        % DISPLAYING RESULTS 
        disp( " ")
        disp(" ++++ CS - VLA 441 - CASE (a)(2) ++++ ");
        disp(" ------------------------------------ ");

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.value];
        disp(" ++++++++++ Vertical tail loads - CASE (a)(2) [daN] ++++++++++ ")
        format = ' %6.6f\n';
        label  = ' Y\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

        %% CS - VLA 441 CASE (a)(3)

        % DYNAMIC PRESSURE AT VA
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.Attributes.unit = "Pa";

        % CALCULATION OF (1.3)*(YAW_ANGLE)
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.yaw_angle.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.yaw_angle.Attributes.unit = "deg";

        % LATERAL FORCE COEFFICIENT 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_VTP.value = 0.0233;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_VTP.Attributes.unit = "Non dimensional";

        % LATERAL FORCE COEFFICIENT 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_VTP.value/Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.Attributes.unit = "Non dimensional";

        % LATERAL FORCE - NEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.Attributes.unit = "N";

        % LATERAL FORCE - DECANEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value = (1e-1)*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.Attributes.unit = "daN";

        % LATERAL FORCE - KILOGRAMS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.value = (1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.Attributes.unit = "N";

%         % -------------------------------------------------------------------------
%         % LATERAL FORCE ON A SINGLE FIN - NEWTON
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.value = (0.5)*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.Attributes.unit = "N";
         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.Attributes.cs = " 441(a)(3) ";
% 
%         % LATERAL FORCE ON A SINGLE FIN - DECANEWTON
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value = (0.5)*(1e-1)*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.Attributes.unit = "daN";
         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.Attributes.cs = " 441(a)(3) ";
% 
%         % LATERAL FORCE ON A SINGLE FIN - KILOGRAMS
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.value = (0.5)*(1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.Attributes.unit = "N";
         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.Attributes.cs = " 441(a)(3) ";

        % DISPLAYING RESULTS 
        disp( " ")
        disp(" ++++ CS - VLA 441 - CASE (a)(3) ++++ ");
        disp(" ------------------------------------ ");

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value];
        disp(" ++++++++++ Vertical tail loads - CASE (a)(3) [daN] ++++++++++ ")
        format = ' %f\n';
        label  = ' Y\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

        %% DEFINING CRITICAL LOADS FOR THE VERTICAL FINS
        tl_0 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_decanewton.value;
        % tl_0 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.value;
        tl_1 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.value;
        tl_2 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value;
        if (abs(tl_0)>abs(tl_1)) && (abs(tl_0)>abs(tl_2)) 
            disp(" ++++++++++ Vertical tail loads - MOST CRITICAL LOADS [daN] ++++++++++ ");
            disp(" CRITICAL LOADS PER ---> CS - VLA 441 (a)(1) "); 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.value = tl_0;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.unit = "daN"; 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.flag = " CS - VLA 441 (a)(1) ";
        elseif (abs(tl_1)>abs(tl_0)) && (abs(tl_1)>abs(tl_0))
            disp(" ++++++++++ Vertical tail loads - MOST CRITICAL LOADS [daN] ++++++++++ ");
            disp(" CRITICAL LOADS PER ---> CS - VLA 441 (a)(2) "); 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.value = tl_1;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.unit = "daN"; 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.flag = " CS - VLA 441 (a)(2) ";
        elseif (abs(tl_2)>abs(tl_1)) && (abs(tl_2)>abs(tl_0))    
            disp(" ++++++++++ Vertical tail loads - MOST CRITICAL LOADS [daN] ++++++++++ ");
            disp(" CRITICAL LOADS PER ---> CS - VLA 441 (a)(3) "); 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.value = tl_2;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.unit = "daN"; 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.flag = " CS - VLA 441 (a)(3) ";
        end

        %% CS - VLA 443 GUST LOADS 
        %  
        %  (a) Vertical tail surfaces must be designed to whithstand, in
        %      unaccelerated flight at speed V = VC, lateral gusts of the values
        %      prescribed for VC in CS - VLA 333 (c).
        %
        %  (b) In the absence of a more rational analysis, the gust load must be
        %      computed as follows 
        %
        %             K_gt * U_de * V * a_vt * S_vt
        %      L_vt = -----------------------------
        %                          16.3
        %
        %      where 
        %
        %      U_de  = derived gust velocities (m/s);
        % 
        %      L_vt  = vertical tail loads (daN);
        %     
        %               (0.88) * mu_gt 
        %      K_gt  = --------------- = gust alleviation factor;
        %                5.3 + mu_gt
        %
        %                          2 * M                   K^2
        %      mu_gt = ------------------------------- * ------- = lat. mass ratio;
        %              rho * c_bar_t * g * a_vt * S_vt   (l_t^2)
        %              
        %              where 
        %                
        %              M    = aeroplane mass (kg); 
        %              rho  = air density (kg/m^3); 
        %              S_vt = area of vertical tail (m^2);
        %              l_t  = distance from aeroplane c.g. to lift centre of
        %                     vertical surface (m); 
        %              a_vt = lift curve slope of vertical tail (1/rad); 
        %              V    = aeroplane equivalent speed (m/s); 
        %              K    = radius of gyration in yaw (m);
        %              g    = acceleration due to gravity (m/s^2).
        %
        %  (c) The average loading in figure B5 and the distribution in figure B8
        %      of Appendix B may be used. 

        %% AMC VLA 443 GUST LOADS 
        %  -----------------------------
        %  INTERPRETATIVE MATERIAL AND ACCEPTABLE MEANS OF COMPLIANCE 
        %
        %  For aeroplanes where the horizontal tail is supported by the vertical
        %  tail, the tail surfaces and their supporting structure including the
        %  rear portion of the fuselage should be designed to whithstand the
        %  prescribed loadings on the vertical tail and the roll - moments induced
        %  by the horizontal tail acting in the same direction. 
        %
        %  For T - Tails in the abscence of a more rational analysis, the rolling
        %  moment induced by gust load may be computed as follow 
        % 
        %                             rho0
        %  M_rudder = (0.3) * S_ht * ------ * V * U * b_ht * K_gf
        %                               2
        %    
        %  where 
        %
        %  M_rudder = induced roll - moment at horizontal tail (N * m);
        %  K_gf     = gust factor = 1.2; 
        %  b_ht     = span of horizontal tail (m); 
        %  S_ht     = area of horizontal tail (m^2); 
        %  rho0     = density of air at sea level (kg/m^3);
        %  V        = airspeed of flight (m/s);
        %  U        = gust speed (m/s).

        %% CS - VLA 445 OUTBOARD FINS 
        %
        %  (a) If outboard fins are on the horizonta tail surface, the tail
        %      surfaces must be designed for the maximum horizontal surface load in
        %      combination with the corresponding loads induced on the vertical
        %      surfaces by endplates effects. These induced effects need not be
        %      combined with other vertical surface loads. 
        %
        %  (b) If outboard fins extend above and below the horizontal surface, the 
        %      critical vertical surface loading, the load per unit area as
        %      determined under CS - VLA 441 and 443, must be applied to 
        %      (1) the part of the vertical surfaces above the horizontal surface
        %          with 80% of that loading applied to the part below the 
        %          horizontal surface; and 
        %      (2) the part of the vertical surfaces below the horizontal surface
        %          with 80% of that loading applied to the part above the 
        %          horizontal surface. 
        %
        %  (c) The endplate effects of outboard fins must be taken into account in
        %      applying the yawing conditions of CS - VLA 441 and 443 to the
        %      vertical surfaces in sub - paragraph (b) of this paragraph.

        %% GUST LOADS CALCULATION 

        % GUST LOAD AT VC
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVC.value = 15.24;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVC.Attributes.unit = "m/s";

        % GUST LOAD AT VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVD.value = 7.62;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVD.Attributes.unit = "m/s";

        % VERTICAL TAIL LIFT SLOPE COEFFICIENT
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value = Aircraft.Certification.Aerodynamic_data.Vertical.a_vt.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.Attributes.unit = "1/rad";

        % VERTICAL TAIL SURFACES
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.Attributes.unit = "m^2";

        % ATMOSPHERE PROPERTIES LOCALLY INSTANTIATED 
        h = 3000; % [m]
        [T, acc, p, rho] = atmosisa(h);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.value = rho;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.Attributes.unit = "kg/m^3";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.Attributes.height = "3000 m";

        % RADIUS OF GYRATION 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.value = 0.3; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.Attributes.unit = "m"; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.Attributes.quantity = "Radius of gyration";

        % AIRSPEED - LOAD AT VC
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VC.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VC.Attributes.unit = "m/s";

        % AIRSPEED - LOAD AT VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VD.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VD.Attributes.unit = "m/s";

        % ++++++++++++++++
        % GUST LOADS AT VC 
        % ++++++++++++++++

        % LATERAL MASS RATIO - VC
        mu_g = @(M, rho, mac_vt, g, a_vt, S_vt, K, l_vt) ((2 * M)/(rho * mac_vt * g * a_vt * S_vt))*(K/l_vt)^2;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.mu_g.value = mu_g(Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.value, ...
                                                                                                Aircraft.Geometry.Vertical.MAC.value, ... 
                                                                                                Aircraft.Constants.g.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.value, ...
                                                                                                Aircraft.Geometry.Vertical.l_vt.value);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.mu_g.Attributes.cs = " 443 ";
        % GUST ALLEVIATION FACTOR - VC
        k_g = @(mu_g) ((0.88)*mu_g)/(5.3 + mu_g);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.K_g.value = k_g(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.mu_g.value);

        % GUST LOADS - VC
        Gust_loads = @(K_gt, U_de, V, a_vt, S_vt) (K_gt * U_de * V * a_vt * S_vt)/(16.3);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.value = Gust_loads(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.K_g.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVC.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VC.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.Attributes.unit = "daN";    
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.Attributes.cs = " 445 ";

        % ++++++++++++++++
        % GUST LOADS AT VD 
        % ++++++++++++++++

        % LATERAL MASS RATIO - VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.mu_g.value = mu_g(Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.value, ...
                                                                                                Aircraft.Geometry.Vertical.MAC.value, ... 
                                                                                                Aircraft.Constants.g.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.value, ...
                                                                                                Aircraft.Geometry.Vertical.l_vt.value);
        % GUST ALLEVIATION FACTOR - VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.K_g.value = k_g(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.mu_g.value);

        % GUST LOADS - VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.value = Gust_loads(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.K_g.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVD.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VD.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.Attributes.unit = "daN"; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.Attributes.cs = " 445 ";

        % DISPLAYING RESULTS 
        disp( " ")
        disp(" ++++ CS - VLA 445 ++++ ");
        disp(" ---------------------- ");

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.value, ...
                     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.value];
        disp(" ++++++++++ Vertical tail - GUST LOADS [daN] ++++++++++ ")
        format = ' %6.6f         %6.6f\n';
        label  = ' GUST - VC        GUST - VD\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

        % CRITICAL GUST LOADS 
        gl_0 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.value;
        gl_1 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.value;

        if (abs(gl_0)>abs(gl_1))
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.value = gl_0;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.flag = " at VC";
        elseif (abs(gl_1)>abs(gl_0))
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.value = gl_1;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.flag = " at VD";
        end

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.value];
        disp(" ++++++++++ Vertical tail - CRITICAL GUST LOADS [daN] ++++++++++ ")
        disp( Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.flag );
        format = ' %6.6f\n';
        label  = ' CRITICAL GUST LOADS\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    % DOUBLE FIN
    case 'Double fin'
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
        %% CS - VLA 441 MANOEUVRING LOADS 
        %
        %  (a) At speeds up to VA, the vertical tail surfaces must be designed to
        %      withistand the following conditions.In computing the tail loads, the
        %      yawing velocity may be assumed to be zero
        %      (1) whit the aeroplane in unaccelerated flight at zero yaw, it is
        %          assumed that the rudder control is suddenly displaced to the
        %          maximum deflection, as limited by the control stops or by limit
        %          pilot forces; 
        %      (2) with the rudder deflected as specified in sub-paragraph (a)(1)
        %          of this paragraph, it is assumed that the aeroplane yaws to the
        %          resulting sideslip angle. In lieu of a rational analysis, an
        %          overswing angle equal to 1.3 times the static sideslip angle of
        %          sub-paragraph (a)(3) of this paragraph may be assumed; 
        %      (3) A yaw angle of 15.0 degrees with the rudder control maintained
        %          in the neutral position, except as limited by pilot strength. 
        %
        %  (b) The average loading of Appendix B, B11 and figure B1 of Appendix B 
        %      and the distribution of figures B6, B7, B8 of Appendix B may be used
        %      instead of requirements of sub-paragraphs (a)(2), (a)(1) and (a)(3)
        %      of this paragraph, respectively.
        %
        %  (c) The yaw angles specified in sub-paragraph (a)(3) of this paragraph
        %      may be reduced if the yaw angle chosen for a particular speed 
        %      cannot be exceeded in 
        %      (1) steady sideslip conditions; and 
        %      (2) uncoordinated rolls from steep banks.

        %% AMC VLA 441 MANOEUVRING LOADS 
        %  -----------------------------
        %  INTERPRETATIVE MATERIAL AND ACCEPTABLE MEANS OF COMPLIANCE 
        % 
        %  For aeroplanes where the horizontal tail is supported by the vertical
        %  tail, the tail surfaces and their supporting structure, including the
        %  rear portion of the fuselage, should be designed to whitstand the
        %  prescribed loading on the vertica ltail and the roll-moments induced by
        %  the horizontal tail acting in the same direction.

        %  For T - Tails, in the abscence of a more rational analysis, the rolling
        %  moment induced by deflection of the vertical rudder may be computed as
        %  follows: 
        %                            rho0
        %  M_rudder = (0.3) * S_ht * ---- * beta * V^2 * b_ht 
        %                              2
        %  where 
        %
        %  M_rudder = induced roll - moment at horizontal tail (N * m) 
        %  b_ht     = span of horizonta tail (m) 
        %  V        = airspeed (m/s)
        %  S_ht     = area of horizontal tail (m^2)
        %  rho0     = air densidity at sea level (kg/m^3)
        %  beta     = angle of zero - lift line due to rudder deflection; this
        %             angle can be defined as follow         
        %                         dL
        %                       ------- * eta * f_n
        %                        d eta 
        %             where 
        %               dL
        %             ------- = change of zero - lift angle of eta * f_n = 1
        %              d eta         
        %       
        %             eta     = rudder deflection 
        %             f_n     = effectivity factor in accordance with angle of
        %                       rudder deflection

        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.global_requirement.Attributes.cs  = " 441 ";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.global_requirement.Attributes.amc = " AMC 441 "; 

%         % APPLICABLE YAW ANGLE 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.value = 15.0;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.Attributes.unit = "deg";

        % APPLICABLE LOAD FACTOR 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.load_factor_n.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.load_factor_n.Attributes.unit = "g's";
        % WING LOADING IN PASCALS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.value = Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.load_factor_n.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.Attributes.unit = "Pa"; 
        % AREA OF THE VERTICAL TAIL 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value = Aircraft.Geometry.Vertical.S.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.Attributes.unit = "m^2";
        % VERTICAL TAIL SURFACE RATIO: S_vt / S_wing
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value/Aircraft.Geometry.Wing.S.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.Attributes.unit = " Non dimensional";
        % PRODUCT OF WING LOADING TIMES THE VERTICAL TAIL SURFACES 
        %
        % P = WS * S * (1/g) --> [P] = [kg]
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Reference_mass.value = (1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Reference_mass.Attributes.unit = "kg";
        % VERTICAL TAIL SPAN b_vt
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.b_vt.value = Aircraft.Geometry.Vertical.b.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.b_vt.Attributes.unit = "m";
        % CALCULATING THE RATIO OF WING LOADING IN KILOGRAMS PER UNIT SPAN 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Spanwise_Load_ratio.value = (1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.b_vt.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Spanwise_Load_ratio.Attributes.unit = "kg/m";
        % VERTICAL TAIL ROOT CHORD 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.value = Aircraft.Geometry.Vertical.croot.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.Attributes.unit = "m";
        % QUARTER CHORD DISTANCE: C_ROOT_VT / 4
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_quarter.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.value/4;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_quarter.Attributes.unit = "m";
        % CHORDWISE KILOGRAMS OF LOAD DISTRIBUTION 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Chordwise_Load_ratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Reference_mass.value/Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Chordwise_Load_ratio.Attributes.unit = "kg/m";

        %% CS - VLA 441 CASE (a)(1)
        % A VECTOR CALLED DR
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.dr.value = [0.0 20.0]';
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.dr.Attributes.unit = "Unknown";

%         % LATERAL FORCE COEFFICIENT AT TWO BETA - FROM CFD
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_vector.value = [0.00003615686 0.01291153]';
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_vector.Attributes.unit = "Non dimensional";

%         % LATERAL FORCE COEFFICIENT WHEN BETA = 0 AND DR = 0
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_0.value = 0.000036156860;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_0.Attributes.unit = "Non dimensional";

%         % LATERAL FORCE GRADIENT C_Y / dr
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CYdr.value = 0.000644;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CYdr.Attributes.unit = "1/deg";

        % MAXIMUM RUDDER DEFLECTION 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.value = 30.0;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.Attributes.unit = "deg";

        % LATERAL FORCE COEFFICIENT WHEN THE RUDDER IS AT MAXIMUM DEFLECTION:
        %              d CY
        % CY = CY_0 + ------- * max(delta_rudder)
        %               d r
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_0.value + Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CYdr.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY.Attributes.unit = "Non dimensional";

        % CALCULATION OF: CY / S_ratio
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY.value/Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.Attributes.unit = "Non dimensional";

        % CALCULATION OF LATERAL FORCE

        % DYNAMIC PRESSURE AT VA
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.Attributes.unit = "Pa";

        % LATERAL FORCE 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.value = 2*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.Attributes.unit = "N";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.Attributes.cs = " 441(a)(1) ";

        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_decanewton.value = (1e-1)*2*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_decanewton.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_decanewton.Attributes.cs = " 441(a)(1) ";

        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_kilograms.value = (1/(Aircraft.Constants.g.value))*2*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_kilograms.Attributes.unit = "kg";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_kilograms.Attributes.cs = " 441(a)(1) ";

        % THE TOTAL LATERAL FORCE MUST BE DIVIDED BY TWO, TO HAVE THE AIRLOADS
        % ACTING ON A SINGLE VERTICAL FIN
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.value*0.5;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin.Attributes.unit = "N"; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin.Attributes.cs = " 441(a)(1) ";

        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.value = (1e-1)*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.value*0.5;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.Attributes.unit = "daN"; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.Attributes.cs = " 441(a)(1) ";

        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_kilograms.value = (1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.value*0.5;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_kilograms.Attributes.unit = "kg"; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_kilograms.Attributes.cs = " 441(a)(1) ";

        % DISPLAYING RESULTS 
        disp( " ")
        disp(" ++++ CS - VLA 441 - CASE (a)(1) ++++ ");
        disp(" ------------------------------------ ");

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.value];
        disp(" ++++++++++ Vertical tail loads - CASE (a)(1) [daN] ++++++++++ ")
        format = ' %6.6f\n';
        label  = ' Y\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

        %% CS - VLA 441 CASE (a)(2)

        % DYNAMIC PRESSURE AT VA
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.Attributes.unit = "Pa";

        % CALCULATION OF (1.3)*(YAW_ANGLE)
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.yaw_angle.value = (1.3)*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.yaw_angle.Attributes.unit = "deg";

%         % LATERAL FORCE COEFFICIENT 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_VTP.value = 0.0245;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_VTP.Attributes.unit = "Non dimensional";

        % LATERAL FORCE COEFFICIENT 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_VTP.value/Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.Attributes.unit = "Non dimensional";

        % LATERAL FORCE - NEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.value = 2*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.Attributes.unit = "N";

        % LATERAL FORCE - DECANEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.value = (1e-1)*2*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.Attributes.unit = "daN";

        % LATERAL FORCE - KILOGRAMS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.value = (1/(Aircraft.Constants.g.value))*2*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.Attributes.unit = "N";

        % -------------------------------------------------------------------------
        % LATERAL FORCE ON A SINGLE FIN - NEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.value = (0.5)*2*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.Attributes.unit = "N";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.Attributes.cs = " 441(a)(2) ";

        % LATERAL FORCE ON A SINGLE FIN - DECANEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.value = (0.5)*(1e-1)*2*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.Attributes.cs = " 441(a)(2) ";

        % LATERAL FORCE ON A SINGLE FIN - KILOGRAMS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.value = (0.5)*(1/(Aircraft.Constants.g.value))*2*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.Attributes.unit = "N";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.Attributes.cs = " 441(a)(2) ";

        % DISPLAYING RESULTS 
        disp( " ")
        disp(" ++++ CS - VLA 441 - CASE (a)(2) ++++ ");
        disp(" ------------------------------------ ");

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.value];
        disp(" ++++++++++ Vertical tail loads - CASE (a)(2) [daN] ++++++++++ ")
        format = ' %6.6f\n';
        label  = ' Y\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

        %% CS - VLA 441 CASE (a)(3)

        % DYNAMIC PRESSURE AT VA
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.Attributes.unit = "Pa";

        % CALCULATION OF (1.3)*(YAW_ANGLE)
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.yaw_angle.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.yaw_angle.Attributes.unit = "deg";

%         % LATERAL FORCE COEFFICIENT 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_VTP.value = 0.0233;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_VTP.Attributes.unit = "Non dimensional";

        % LATERAL FORCE COEFFICIENT 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_VTP.value/Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.Attributes.unit = "Non dimensional";

        % LATERAL FORCE - NEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.value = 2*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.Attributes.unit = "N";

        % LATERAL FORCE - DECANEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value = (1e-1)*2*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.Attributes.unit = "daN";

        % LATERAL FORCE - KILOGRAMS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.value = (1/(Aircraft.Constants.g.value))*2*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.Attributes.unit = "N";

        % -------------------------------------------------------------------------
        % LATERAL FORCE ON A SINGLE FIN - NEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.value = (0.5)*2*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.Attributes.unit = "N";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.Attributes.cs = " 441(a)(3) ";

        % LATERAL FORCE ON A SINGLE FIN - DECANEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value = (0.5)*(1e-1)*2*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.Attributes.cs = " 441(a)(3) ";

        % LATERAL FORCE ON A SINGLE FIN - KILOGRAMS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.value = (0.5)*(1/(Aircraft.Constants.g.value))*2*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.Attributes.unit = "N";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.Attributes.cs = " 441(a)(3) ";

        % DISPLAYING RESULTS 
        disp( " ")
        disp(" ++++ CS - VLA 441 - CASE (a)(3) ++++ ");
        disp(" ------------------------------------ ");

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value];
        disp(" ++++++++++ Vertical tail loads - CASE (a)(3) [daN] ++++++++++ ")
        format = ' %f\n';
        label  = ' Y\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

        %% DEFINING CRITICAL LOADS FOR THE VERTICAL FINS
        tl_0 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.value;
        tl_1 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.value;
        tl_2 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value;
        if (abs(tl_0)>abs(tl_1)) && (abs(tl_0)>abs(tl_2)) 
            disp(" ++++++++++ Vertical tail loads - MOST CRITICAL LOADS [daN] ++++++++++ ");
            disp(" CRITICAL LOADS PER ---> CS - VLA 441 (a)(1) "); 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.value = tl_0;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.unit = "daN"; 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.flag = " CS - VLA 441 (a)(1) ";
        elseif (abs(tl_1)>abs(tl_0)) && (abs(tl_1)>abs(tl_0))
            disp(" ++++++++++ Vertical tail loads - MOST CRITICAL LOADS [daN] ++++++++++ ");
            disp(" CRITICAL LOADS PER ---> CS - VLA 441 (a)(2) "); 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.value = tl_1;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.unit = "daN"; 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.flag = " CS - VLA 441 (a)(2) ";
        elseif (abs(tl_2)>abs(tl_1)) && (abs(tl_2)>abs(tl_0))    
            disp(" ++++++++++ Vertical tail loads - MOST CRITICAL LOADS [daN] ++++++++++ ");
            disp(" CRITICAL LOADS PER ---> CS - VLA 441 (a)(3) "); 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.value = tl_2;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.unit = "daN"; 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.flag = " CS - VLA 441 (a)(3) ";
        end

        %% CS - VLA 443 GUST LOADS 
        %  
        %  (a) Vertical tail surfaces must be designed to whithstand, in
        %      unaccelerated flight at speed V = VC, lateral gusts of the values
        %      prescribed for VC in CS - VLA 333 (c).
        %
        %  (b) In the absence of a more rational analysis, the gust load must be
        %      computed as follows 
        %
        %             K_gt * U_de * V * a_vt * S_vt
        %      L_vt = -----------------------------
        %                          16.3
        %
        %      where 
        %
        %      U_de  = derived gust velocities (m/s);
        % 
        %      L_vt  = vertical tail loads (daN);
        %     
        %               (0.88) * mu_gt 
        %      K_gt  = --------------- = gust alleviation factor;
        %                5.3 + mu_gt
        %
        %                          2 * M                   K^2
        %      mu_gt = ------------------------------- * ------- = lat. mass ratio;
        %              rho * c_bar_t * g * a_vt * S_vt   (l_t^2)
        %              
        %              where 
        %                
        %              M    = aeroplane mass (kg); 
        %              rho  = air density (kg/m^3); 
        %              S_vt = area of vertical tail (m^2);
        %              l_t  = distance from aeroplane c.g. to lift centre of
        %                     vertical surface (m); 
        %              a_vt = lift curve slope of vertical tail (1/rad); 
        %              V    = aeroplane equivalent speed (m/s); 
        %              K    = radius of gyration in yaw (m);
        %              g    = acceleration due to gravity (m/s^2).
        %
        %  (c) The average loading in figure B5 and the distribution in figure B8
        %      of Appendix B may be used. 

        %% AMC VLA 443 GUST LOADS 
        %  -----------------------------
        %  INTERPRETATIVE MATERIAL AND ACCEPTABLE MEANS OF COMPLIANCE 
        %
        %  For aeroplanes where the horizontal tail is supported by the vertical
        %  tail, the tail surfaces and their supporting structure including the
        %  rear portion of the fuselage should be designed to whithstand the
        %  prescribed loadings on the vertical tail and the roll - moments induced
        %  by the horizontal tail acting in the same direction. 
        %
        %  For T - Tails in the abscence of a more rational analysis, the rolling
        %  moment induced by gust load may be computed as follow 
        % 
        %                             rho0
        %  M_rudder = (0.3) * S_ht * ----- * V * U * b_ht * K_gf
        %                              2
        %    
        %  where 
        %
        %  M_rudder = induced roll - moment at horizontal tail (N * m);
        %  K_gf     = gust factor = 1.2; 
        %  b_ht     = span of horizontal tail (m); 
        %  S_ht     = area of horizontal tail (m^2); 
        %  rho0     = density of air at sea level (kg/m^3);
        %  V        = airspeed of flight (m/s);
        %  U        = gust speed (m/s).

        %% CS - VLA 445 OUTBOARD FINS 
        %
        %  (a) If outboard fins are on the horizonta tail surface, the tail
        %      surfaces must be designed for the maximum horizontal surface load in
        %      combination with the corresponding loads induced on the vertical
        %      surfaces by endplates effects. These induced effects need not be
        %      combined with other vertical surface loads. 
        %
        %  (b) If outboard fins extend above and below the horizontal surface, the 
        %      critical vertical surface loading, the load per unit area as
        %      determined under CS - VLA 441 and 443, must be applied to 
        %      (1) the part of the vertical surfaces above the horizontal surface
        %          with 80% of that loading applied to the part below the 
        %          horizontal surface; and 
        %      (2) the part of the vertical surfaces below the horizontal surface
        %          with 80% of that loading applied to the part above the 
        %          horizontal surface. 
        %
        %  (c) The endplate effects of outboard fins must be taken into account in
        %      applying the yawing conditions of CS - VLA 441 and 443 to the
        %      vertical surfaces in sub - paragraph (b) of this paragraph.

        %% GUST LOADS CALCULATION 

        % GUST LOAD AT VC
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVC.value = 15.24;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVC.Attributes.unit = "m/s";

        % GUST LOAD AT VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVD.value = 7.62;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVD.Attributes.unit = "m/s";

        % VERTICAL TAIL LIFT SLOPE COEFFICIENT
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value = Aircraft.Certification.Aerodynamic_data.Vertical.a_vt.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.Attributes.unit = "1/rad";

        % VERTICAL TAIL SURFACES
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.Attributes.unit = "m^2";

        % ATMOSPHERE PROPERTIES LOCALLY INSTANTIATED 
        % h = 3000; % [m]
        h = Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.value;
        [T, acc, p, rho] = atmosisa(h);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.value = rho;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.Attributes.unit = "kg/m^3";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.Attributes.height = "3500 m";

%         % RADIUS OF GYRATION 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.value = 0.3; 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.Attributes.unit = "m"; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.Attributes.quantity = "Radius of gyration";

        % AIRSPEED - LOAD AT VC
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VC.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VC.Attributes.unit = "m/s";

        % AIRSPEED - LOAD AT VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VD.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VD.Attributes.unit = "m/s";

        % ++++++++++++++++
        % GUST LOADS AT VC 
        % ++++++++++++++++

        % LATERAL MASS RATIO - VC
        mu_g = @(M, rho, mac_vt, g, a_vt, S_vt, K, l_vt) ((2 * M)/(rho * mac_vt * g * a_vt * S_vt))*(K/l_vt)^2;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.mu_g.value = mu_g(Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.value, ...
                                                                                                Aircraft.Geometry.Vertical.MAC.value, ... 
                                                                                                Aircraft.Constants.g.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.value, ...
                                                                                                Aircraft.Geometry.Vertical.l_vt.value);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.mu_g.Attributes.cs = " 443 ";
        % GUST ALLEVIATION FACTOR - VC
        k_g = @(mu_g) ((0.88)*mu_g)/(5.3 + mu_g);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.K_g.value = k_g(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.mu_g.value);

        % GUST LOADS - VC
        Gust_loads = @(K_gt, U_de, V, a_vt, S_vt) (K_gt * U_de * V * a_vt * S_vt)/(16.3);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.value = Gust_loads(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.K_g.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVC.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VC.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.Attributes.unit = "daN";    
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.Attributes.cs = " 445 ";

        % ++++++++++++++++
        % GUST LOADS AT VD 
        % ++++++++++++++++

        % LATERAL MASS RATIO - VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.mu_g.value = mu_g(Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.value, ...
                                                                                                Aircraft.Geometry.Vertical.MAC.value, ... 
                                                                                                Aircraft.Constants.g.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.value, ...
                                                                                                Aircraft.Geometry.Vertical.l_vt.value);
        % GUST ALLEVIATION FACTOR - VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.K_g.value = k_g(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.mu_g.value);

        % GUST LOADS - VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.value = Gust_loads(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.K_g.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVD.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VD.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.Attributes.unit = "daN"; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.Attributes.cs = " 445 ";

        % DISPLAYING RESULTS 
        disp( " ")
        disp(" ++++ CS - VLA 445 ++++ ");
        disp(" ---------------------- ");

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.value, ...
                     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.value];
        disp(" ++++++++++ Vertical tail - GUST LOADS [daN] ++++++++++ ")
        format = ' %6.6f         %6.6f\n';
        label  = ' GUST - VC        GUST - VD\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

        % CRITICAL GUST LOADS 
        gl_0 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.value;
        gl_1 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.value;

        if (abs(gl_0)>abs(gl_1))
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.value = gl_0;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.flag = " at VC";
        elseif (abs(gl_1)>abs(gl_0))
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.value = gl_1;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.flag = " at VD";
        end

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.value];
        disp(" ++++++++++ Vertical tail - CRITICAL GUST LOADS [daN] ++++++++++ ")
        disp( Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.flag );
        format = ' %6.6f\n';
        label  = ' CRITICAL GUST LOADS\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % MULTIPLE FINS
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    case 'Multiple fin'
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
        %% CS - VLA 441 MANOEUVRING LOADS 
        %
        %  (a) At speeds up to VA, the vertical tail surfaces must be designed to
        %      withistand the following conditions.In computing the tail loads, the
        %      yawing velocity may be assumed to be zero
        %      (1) whit the aeroplane in unaccelerated flight at zero yaw, it is
        %          assumed that the rudder control is suddenly displaced to the
        %          maximum deflection, as limited by the control stops or by limit
        %          pilot forces; 
        %      (2) with the rudder deflected as specified in sub-paragraph (a)(1)
        %          of this paragraph, it is assumed that the aeroplane yaws to the
        %          resulting sideslip angle. In lieu of a rational analysis, an
        %          overswing angle equal to 1.3 times the static sideslip angle of
        %          sub-paragraph (a)(3) of this paragraph may be assumed; 
        %      (3) A yaw angle of 15.0 degrees with the rudder control maintained
        %          in the neutral position, except as limited by pilot strength. 
        %
        %  (b) The average loading of Appendix B, B11 and figure B1 of Appendix B 
        %      and the distribution of figures B6, B7, B8 of Appendix B may be used
        %      instead of requirements of sub-paragraphs (a)(2), (a)(1) and (a)(3)
        %      of this paragraph, respectively.
        %
        %  (c) The yaw angles specified in sub-paragraph (a)(3) of this paragraph
        %      may be reduced if the yaw angle chosen for a particular speed 
        %      cannot be exceeded in 
        %      (1) steady sideslip conditions; and 
        %      (2) uncoordinated rolls from steep banks.

        %% AMC VLA 441 MANOEUVRING LOADS 
        %  -----------------------------
        %  INTERPRETATIVE MATERIAL AND ACCEPTABLE MEANS OF COMPLIANCE 
        % 
        %  For aeroplanes where the horizontal tail is supported by the vertical
        %  tail, the tail surfaces and their supporting structure, including the
        %  rear portion of the fuselage, should be designed to whitstand the
        %  prescribed loading on the vertica ltail and the roll-moments induced by
        %  the horizontal tail acting in the same direction.

        %  For T - Tails, in the abscence of a more rational analysis, the rolling
        %  moment induced by deflection of the vertical rudder may be computed as
        %  follows: 
        %                            rho0
        %  M_rudder = (0.3) * S_ht * ---- * beta * V^2 * b_ht 
        %                              2
        %  where 
        %
        %  M_rudder = induced roll - moment at horizontal tail (N * m) 
        %  b_ht     = span of horizonta tail (m) 
        %  V        = airspeed (m/s)
        %  S_ht     = area of horizontal tail (m^2)
        %  rho0     = air densidity at sea level (kg/m^3)
        %  beta     = angle of zero - lift line due to rudder deflection; this
        %             angle can be defined as follow         
        %                         dL
        %                       ------- * eta * f_n
        %                        d eta 
        %             where 
        %               dL
        %             ------- = change of zero - lift angle of eta * f_n = 1
        %              d eta         
        %       
        %             eta     = rudder deflection 
        %             f_n     = effectivity factor in accordance with angle of
        %                       rudder deflection

        if Aircraft.Geometry.Vertical.empennage_flag.Attributes.number_of_fin ~ NaN 
            n_fin = Aircraft.Geometry.Vertical.empennage_flag.Attributes.number_of_fin;
        end
        
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.global_requirement.Attributes.cs  = " 441 ";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.global_requirement.Attributes.amc = " AMC 441 "; 

%         % APPLICABLE YAW ANGLE 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.value = 15.0;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.Attributes.unit = "deg";

        % APPLICABLE LOAD FACTOR 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.load_factor_n.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.load_factor_n.Attributes.unit = "g's";
        % WING LOADING IN PASCALS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.value = Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.load_factor_n.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.Attributes.unit = "Pa"; 
        % AREA OF THE VERTICAL TAIL 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value = Aircraft.Geometry.Vertical.S.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.Attributes.unit = "m^2";
        % VERTICAL TAIL SURFACE RATIO: S_vt / S_wing
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value/Aircraft.Geometry.Wing.S.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.Attributes.unit = " Non dimensional";
        % PRODUCT OF WING LOADING TIMES THE VERTICAL TAIL SURFACES 
        %
        % P = WS * S * (1/g) --> [P] = [kg]
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Reference_mass.value = (1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Reference_mass.Attributes.unit = "kg";
        % VERTICAL TAIL SPAN b_vt
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.b_vt.value = Aircraft.Geometry.Vertical.b.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.b_vt.Attributes.unit = "m";
        % CALCULATING THE RATIO OF WING LOADING IN KILOGRAMS PER UNIT SPAN 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Spanwise_Load_ratio.value = (1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.b_vt.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Spanwise_Load_ratio.Attributes.unit = "kg/m";
        % VERTICAL TAIL ROOT CHORD 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.value = Aircraft.Geometry.Vertical.croot.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.Attributes.unit = "m";
        % QUARTER CHORD DISTANCE: C_ROOT_VT / 4
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_quarter.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.value/4;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_quarter.Attributes.unit = "m";
        % CHORDWISE KILOGRAMS OF LOAD DISTRIBUTION 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Chordwise_Load_ratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Reference_mass.value/Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Chordwise_Load_ratio.Attributes.unit = "kg/m";

        %% CS - VLA 441 CASE (a)(1)
        % A VECTOR CALLED DR
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.dr.value = [0.0 20.0]';
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.dr.Attributes.unit = "Unknown";

%         % LATERAL FORCE COEFFICIENT AT TWO BETA - FROM CFD
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_vector.value = [0.00003615686 0.01291153]';
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_vector.Attributes.unit = "Non dimensional";

%         % LATERAL FORCE COEFFICIENT WHEN BETA = 0 AND DR = 0
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_0.value = 0.000036156860;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_0.Attributes.unit = "Non dimensional";

%         % LATERAL FORCE GRADIENT C_Y / dr
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CYdr.value = 0.000644;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CYdr.Attributes.unit = "1/deg";

        % MAXIMUM RUDDER DEFLECTION 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.value = 30.0;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.Attributes.unit = "deg";

        % LATERAL FORCE COEFFICIENT WHEN THE RUDDER IS AT MAXIMUM DEFLECTION:
        %              d CY
        % CY = CY_0 + ------- * max(delta_rudder)
        %               d r
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_0.value + Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CYdr.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY.Attributes.unit = "Non dimensional";

        % CALCULATION OF: CY / S_ratio
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY.value/Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.Attributes.unit = "Non dimensional";

        % CALCULATION OF LATERAL FORCE

        % DYNAMIC PRESSURE AT VA
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.Attributes.unit = "Pa";

        % LATERAL FORCE 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.value = n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.Attributes.unit = "N";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.Attributes.cs = " 441(a)(1) ";

        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_decanewton.value = (1e-1)*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_decanewton.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_decanewton.Attributes.cs = " 441(a)(1) ";

        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_kilograms.value = (1/(Aircraft.Constants.g.value))*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_kilograms.Attributes.unit = "kg";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_kilograms.Attributes.cs = " 441(a)(1) ";

        % THE TOTAL LATERAL FORCE MUST BE DIVIDED BY TWO, TO HAVE THE AIRLOADS
        % ACTING ON A SINGLE VERTICAL FIN
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.value / n_fin;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin.Attributes.unit = "N"; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin.Attributes.cs = " 441(a)(1) ";

        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.value = (1e-1)*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.value / n_fin;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.Attributes.unit = "daN"; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.Attributes.cs = " 441(a)(1) ";

        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_kilograms.value = (1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.value / n_fin;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_kilograms.Attributes.unit = "kg"; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_kilograms.Attributes.cs = " 441(a)(1) ";

        % DISPLAYING RESULTS 
        disp( " ")
        disp(" ++++ CS - VLA 441 - CASE (a)(1) ++++ ");
        disp(" ------------------------------------ ");

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.value];
        disp(" ++++++++++ Vertical tail loads - CASE (a)(1) [daN] ++++++++++ ")
        format = ' %6.6f\n';
        label  = ' Y\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

        %% CS - VLA 441 CASE (a)(2)

        % DYNAMIC PRESSURE AT VA
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.Attributes.unit = "Pa";

        % CALCULATION OF (1.3)*(YAW_ANGLE)
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.yaw_angle.value = (1.3)*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.yaw_angle.Attributes.unit = "deg";

%         % LATERAL FORCE COEFFICIENT 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_VTP.value = 0.0245;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_VTP.Attributes.unit = "Non dimensional";

        % LATERAL FORCE COEFFICIENT 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_VTP.value/Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.Attributes.unit = "Non dimensional";

        % LATERAL FORCE - NEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.value = n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.Attributes.unit = "N";

        % LATERAL FORCE - DECANEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.value = (1e-1)*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.Attributes.unit = "daN";

        % LATERAL FORCE - KILOGRAMS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.value = (1/(Aircraft.Constants.g.value))*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.Attributes.unit = "N";

        % -------------------------------------------------------------------------
        % LATERAL FORCE ON A SINGLE FIN - NEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.value = (0.5)*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.Attributes.unit = "N";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.Attributes.cs = " 441(a)(2) ";

        % LATERAL FORCE ON A SINGLE FIN - DECANEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.value = (0.5)*(1e-1)*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.Attributes.cs = " 441(a)(2) ";

        % LATERAL FORCE ON A SINGLE FIN - KILOGRAMS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.value = (0.5)*(1/(Aircraft.Constants.g.value))*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.Attributes.unit = "N";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.Attributes.cs = " 441(a)(2) ";

        % DISPLAYING RESULTS 
        disp( " ")
        disp(" ++++ CS - VLA 441 - CASE (a)(2) ++++ ");
        disp(" ------------------------------------ ");

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.value];
        disp(" ++++++++++ Vertical tail loads - CASE (a)(2) [daN] ++++++++++ ")
        format = ' %6.6f\n';
        label  = ' Y\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

        %% CS - VLA 441 CASE (a)(3)

        % DYNAMIC PRESSURE AT VA
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.Attributes.unit = "Pa";

        % CALCULATION OF (1.3)*(YAW_ANGLE)
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.yaw_angle.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.yaw_angle.Attributes.unit = "deg";

%         % LATERAL FORCE COEFFICIENT 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_VTP.value = 0.0233;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_VTP.Attributes.unit = "Non dimensional";

        % LATERAL FORCE COEFFICIENT 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_VTP.value/Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.Attributes.unit = "Non dimensional";

        % LATERAL FORCE - NEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.value = n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.Attributes.unit = "N";

        % LATERAL FORCE - DECANEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value = (1e-1)*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.Attributes.unit = "daN";

        % LATERAL FORCE - KILOGRAMS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.value = (1/(Aircraft.Constants.g.value))*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.Attributes.unit = "N";

        % -------------------------------------------------------------------------
        % LATERAL FORCE ON A SINGLE FIN - NEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.value = (0.5)*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.Attributes.unit = "N";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.Attributes.cs = " 441(a)(3) ";

        % LATERAL FORCE ON A SINGLE FIN - DECANEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value = (0.5)*(1e-1)*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.Attributes.cs = " 441(a)(3) ";

        % LATERAL FORCE ON A SINGLE FIN - KILOGRAMS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.value = (0.5)*(1/(Aircraft.Constants.g.value))*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.Attributes.unit = "N";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.Attributes.cs = " 441(a)(3) ";

        % DISPLAYING RESULTS 
        disp( " ")
        disp(" ++++ CS - VLA 441 - CASE (a)(3) ++++ ");
        disp(" ------------------------------------ ");

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value];
        disp(" ++++++++++ Vertical tail loads - CASE (a)(3) [daN] ++++++++++ ")
        format = ' %f\n';
        label  = ' Y\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

        %% DEFINING CRITICAL LOADS FOR THE VERTICAL FINS
        tl_0 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.value;
        tl_1 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.value;
        tl_2 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value;
        if (abs(tl_0)>abs(tl_1)) && (abs(tl_0)>abs(tl_2)) 
            disp(" ++++++++++ Vertical tail loads - MOST CRITICAL LOADS [daN] ++++++++++ ");
            disp(" CRITICAL LOADS PER ---> CS - VLA 441 (a)(1) "); 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.value = tl_0;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.unit = "daN"; 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.flag = " CS - VLA 441 (a)(1) ";
        elseif (abs(tl_1)>abs(tl_0)) && (abs(tl_1)>abs(tl_0))
            disp(" ++++++++++ Vertical tail loads - MOST CRITICAL LOADS [daN] ++++++++++ ");
            disp(" CRITICAL LOADS PER ---> CS - VLA 441 (a)(2) "); 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.value = tl_1;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.unit = "daN"; 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.flag = " CS - VLA 441 (a)(2) ";
        elseif (abs(tl_2)>abs(tl_1)) && (abs(tl_2)>abs(tl_0))    
            disp(" ++++++++++ Vertical tail loads - MOST CRITICAL LOADS [daN] ++++++++++ ");
            disp(" CRITICAL LOADS PER ---> CS - VLA 441 (a)(3) "); 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.value = tl_2;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.unit = "daN"; 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.flag = " CS - VLA 441 (a)(3) ";
        end

        %% CS - VLA 443 GUST LOADS 
        %  
        %  (a) Vertical tail surfaces must be designed to whithstand, in
        %      unaccelerated flight at speed V = VC, lateral gusts of the values
        %      prescribed for VC in CS - VLA 333 (c).
        %
        %  (b) In the absence of a more rational analysis, the gust load must be
        %      computed as follows 
        %
        %             K_gt * U_de * V * a_vt * S_vt
        %      L_vt = -----------------------------
        %                          16.3
        %
        %      where 
        %
        %      U_de  = derived gust velocities (m/s);
        % 
        %      L_vt  = vertical tail loads (daN);
        %     
        %               (0.88) * mu_gt 
        %      K_gt  = --------------- = gust alleviation factor;
        %                5.3 + mu_gt
        %
        %                          2 * M                   K^2
        %      mu_gt = ------------------------------- * ------- = lat. mass ratio;
        %              rho * c_bar_t * g * a_vt * S_vt   (l_t^2)
        %              
        %              where 
        %                
        %              M    = aeroplane mass (kg); 
        %              rho  = air density (kg/m^3); 
        %              S_vt = area of vertical tail (m^2);
        %              l_t  = distance from aeroplane c.g. to lift centre of
        %                     vertical surface (m); 
        %              a_vt = lift curve slope of vertical tail (1/rad); 
        %              V    = aeroplane equivalent speed (m/s); 
        %              K    = radius of gyration in yaw (m);
        %              g    = acceleration due to gravity (m/s^2).
        %
        %  (c) The average loading in figure B5 and the distribution in figure B8
        %      of Appendix B may be used. 

        %% AMC VLA 443 GUST LOADS 
        %  -----------------------------
        %  INTERPRETATIVE MATERIAL AND ACCEPTABLE MEANS OF COMPLIANCE 
        %
        %  For aeroplanes where the horizontal tail is supported by the vertical
        %  tail, the tail surfaces and their supporting structure including the
        %  rear portion of the fuselage should be designed to whithstand the
        %  prescribed loadings on the vertical tail and the roll - moments induced
        %  by the horizontal tail acting in the same direction. 
        %
        %  For T - Tails in the abscence of a more rational analysis, the rolling
        %  moment induced by gust load may be computed as follow 
        % 
        %                             rho0
        %  M_rudder = (0.3) * S_ht * ----- * V * U * b_ht * K_gf
        %                              2
        %    
        %  where 
        %
        %  M_rudder = induced roll - moment at horizontal tail (N * m);
        %  K_gf     = gust factor = 1.2; 
        %  b_ht     = span of horizontal tail (m); 
        %  S_ht     = area of horizontal tail (m^2); 
        %  rho0     = density of air at sea level (kg/m^3);
        %  V        = airspeed of flight (m/s);
        %  U        = gust speed (m/s).

        %% CS - VLA 445 OUTBOARD FINS 
        %
        %  (a) If outboard fins are on the horizonta tail surface, the tail
        %      surfaces must be designed for the maximum horizontal surface load in
        %      combination with the corresponding loads induced on the vertical
        %      surfaces by endplates effects. These induced effects need not be
        %      combined with other vertical surface loads. 
        %
        %  (b) If outboard fins extend above and below the horizontal surface, the 
        %      critical vertical surface loading, the load per unit area as
        %      determined under CS - VLA 441 and 443, must be applied to 
        %      (1) the part of the vertical surfaces above the horizontal surface
        %          with 80% of that loading applied to the part below the 
        %          horizontal surface; and 
        %      (2) the part of the vertical surfaces below the horizontal surface
        %          with 80% of that loading applied to the part above the 
        %          horizontal surface. 
        %
        %  (c) The endplate effects of outboard fins must be taken into account in
        %      applying the yawing conditions of CS - VLA 441 and 443 to the
        %      vertical surfaces in sub - paragraph (b) of this paragraph.

        %% GUST LOADS CALCULATION 

        % GUST LOAD AT VC
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVC.value = 15.24;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVC.Attributes.unit = "m/s";

        % GUST LOAD AT VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVD.value = 7.62;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVD.Attributes.unit = "m/s";

        % VERTICAL TAIL LIFT SLOPE COEFFICIENT
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value = Aircraft.Certification.Aerodynamic_data.Vertical.a_vt.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.Attributes.unit = "1/rad";

        % VERTICAL TAIL SURFACES
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.Attributes.unit = "m^2";

        % ATMOSPHERE PROPERTIES LOCALLY INSTANTIATED 
        % h = 3000; % [m]
        h = Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.value;
        [T, acc, p, rho] = atmosisa(h);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.value = rho;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.Attributes.unit = "kg/m^3";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.Attributes.height = "3500 m";

%         % RADIUS OF GYRATION 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.value = 0.3; 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.Attributes.unit = "m"; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.Attributes.quantity = "Radius of gyration";

        % AIRSPEED - LOAD AT VC
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VC.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VC.Attributes.unit = "m/s";

        % AIRSPEED - LOAD AT VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VD.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VD.Attributes.unit = "m/s";

        % ++++++++++++++++
        % GUST LOADS AT VC 
        % ++++++++++++++++

        % LATERAL MASS RATIO - VC
        mu_g = @(M, rho, mac_vt, g, a_vt, S_vt, K, l_vt) ((2 * M)/(rho * mac_vt * g * a_vt * S_vt))*(K/l_vt)^2;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.mu_g.value = mu_g(Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.value, ...
                                                                                                Aircraft.Geometry.Vertical.MAC.value, ... 
                                                                                                Aircraft.Constants.g.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.value, ...
                                                                                                Aircraft.Geometry.Vertical.l_vt.value);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.mu_g.Attributes.cs = " 443 ";
        % GUST ALLEVIATION FACTOR - VC
        k_g = @(mu_g) ((0.88)*mu_g)/(5.3 + mu_g);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.K_g.value = k_g(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.mu_g.value);

        % GUST LOADS - VC
        Gust_loads = @(K_gt, U_de, V, a_vt, S_vt) (K_gt * U_de * V * a_vt * S_vt)/(16.3);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.value = Gust_loads(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.K_g.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVC.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VC.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.Attributes.unit = "daN";    
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.Attributes.cs = " 445 ";

        % ++++++++++++++++
        % GUST LOADS AT VD 
        % ++++++++++++++++

        % LATERAL MASS RATIO - VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.mu_g.value = mu_g(Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.value, ...
                                                                                                Aircraft.Geometry.Vertical.MAC.value, ... 
                                                                                                Aircraft.Constants.g.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.value, ...
                                                                                                Aircraft.Geometry.Vertical.l_vt.value);
        % GUST ALLEVIATION FACTOR - VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.K_g.value = k_g(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.mu_g.value);

        % GUST LOADS - VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.value = Gust_loads(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.K_g.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVD.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VD.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.Attributes.unit = "daN"; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.Attributes.cs = " 445 ";

        % DISPLAYING RESULTS 
        disp( " ")
        disp(" ++++ CS - VLA 445 ++++ ");
        disp(" ---------------------- ");

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.value, ...
                     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.value];
        disp(" ++++++++++ Vertical tail - GUST LOADS [daN] ++++++++++ ")
        format = ' %6.6f         %6.6f\n';
        label  = ' GUST - VC        GUST - VD\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

        % CRITICAL GUST LOADS 
        gl_0 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.value;
        gl_1 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.value;

        if (abs(gl_0)>abs(gl_1))
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.value = gl_0;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.flag = " at VC";
        elseif (abs(gl_1)>abs(gl_0))
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.value = gl_1;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.flag = " at VD";
        end

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.value];
        disp(" ++++++++++ Vertical tail - CRITICAL GUST LOADS [daN] ++++++++++ ")
        disp( Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.flag );
        format = ' %6.6f\n';
        label  = ' CRITICAL GUST LOADS\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")      

    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % T TAIL 
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    case 'T-tail'
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
        %% CS - VLA 441 MANOEUVRING LOADS 
        %
        %  (a) At speeds up to VA, the vertical tail surfaces must be designed to
        %      withistand the following conditions.In computing the tail loads, the
        %      yawing velocity may be assumed to be zero
        %      (1) whit the aeroplane in unaccelerated flight at zero yaw, it is
        %          assumed that the rudder control is suddenly displaced to the
        %          maximum deflection, as limited by the control stops or by limit
        %          pilot forces; 
        %      (2) with the rudder deflected as specified in sub-paragraph (a)(1)
        %          of this paragraph, it is assumed that the aeroplane yaws to the
        %          resulting sideslip angle. In lieu of a rational analysis, an
        %          overswing angle equal to 1.3 times the static sideslip angle of
        %          sub-paragraph (a)(3) of this paragraph may be assumed; 
        %      (3) A yaw angle of 15.0 degrees with the rudder control maintained
        %          in the neutral position, except as limited by pilot strength. 
        %
        %  (b) The average loading of Appendix B, B11 and figure B1 of Appendix B 
        %      and the distribution of figures B6, B7, B8 of Appendix B may be used
        %      instead of requirements of sub-paragraphs (a)(2), (a)(1) and (a)(3)
        %      of this paragraph, respectively.
        %
        %  (c) The yaw angles specified in sub-paragraph (a)(3) of this paragraph
        %      may be reduced if the yaw angle chosen for a particular speed 
        %      cannot be exceeded in 
        %      (1) steady sideslip conditions; and 
        %      (2) uncoordinated rolls from steep banks.

        %% AMC VLA 441 MANOEUVRING LOADS 
        %  -----------------------------
        %  INTERPRETATIVE MATERIAL AND ACCEPTABLE MEANS OF COMPLIANCE 
        % 
        %  For aeroplanes where the horizontal tail is supported by the vertical
        %  tail, the tail surfaces and their supporting structure, including the
        %  rear portion of the fuselage, should be designed to whitstand the
        %  prescribed loading on the vertica ltail and the roll-moments induced by
        %  the horizontal tail acting in the same direction.

        %  For T - Tails, in the abscence of a more rational analysis, the rolling
        %  moment induced by deflection of the vertical rudder may be computed as
        %  follows: 
        %                            rho0
        %  M_rudder = (0.3) * S_ht * ---- * beta * V^2 * b_ht 
        %                              2
        %  where 
        %
        %  M_rudder = induced roll - moment at horizontal tail (N * m) 
        %  b_ht     = span of horizonta tail (m) 
        %  V        = airspeed (m/s)
        %  S_ht     = area of horizontal tail (m^2)
        %  rho0     = air densidity at sea level (kg/m^3)
        %  beta     = angle of zero - lift line due to rudder deflection; this
        %             angle can be defined as follow         
        %                         dL
        %                       ------- * eta * f_n
        %                        d eta 
        %             where 
        %               dL
        %             ------- = change of zero - lift angle of eta * f_n = 1
        %              d eta         
        %       
        %             eta     = rudder deflection 
        %             f_n     = effectivity factor in accordance with angle of
        %                       rudder deflection
        
        % MANOEUVRING AIR SPEED
        VA     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
        VC     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
        VD     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
        S_ht   = Aircraft.Geometry.Horizontal.S.value;
        b_ht   = Aircraft.Geometry.Horizontal.b.value;
        rho0   = Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value;
        dLdeta = 1.0;
        f_eta  = 0.7;

        if Aircraft.Geometry.Vertical.empennage_flag.Attributes.number_of_fin ~ NaN 
            n_fin = Aircraft.Geometry.Vertical.empennage_flag.Attributes.number_of_fin;
        end
        
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.global_requirement.Attributes.cs  = " 441 ";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.global_requirement.Attributes.amc = " AMC 441 "; 

%         % APPLICABLE YAW ANGLE 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.value = 15.0;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.Attributes.unit = "deg";

        % APPLICABLE LOAD FACTOR 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.load_factor_n.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.load_factor_n.Attributes.unit = "g's";
        % WING LOADING IN PASCALS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.value = Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.load_factor_n.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.Attributes.unit = "Pa"; 
        % AREA OF THE VERTICAL TAIL 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value = Aircraft.Geometry.Vertical.S.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.Attributes.unit = "m^2";
        % VERTICAL TAIL SURFACE RATIO: S_vt / S_wing
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value/Aircraft.Geometry.Wing.S.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.Attributes.unit = " Non dimensional";
        % PRODUCT OF WING LOADING TIMES THE VERTICAL TAIL SURFACES 
        %
        % P = WS * S * (1/g) --> [P] = [kg]
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Reference_mass.value = (1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Reference_mass.Attributes.unit = "kg";
        % VERTICAL TAIL SPAN b_vt
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.b_vt.value = Aircraft.Geometry.Vertical.b.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.b_vt.Attributes.unit = "m";
        % CALCULATING THE RATIO OF WING LOADING IN KILOGRAMS PER UNIT SPAN 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Spanwise_Load_ratio.value = (1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.b_vt.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Spanwise_Load_ratio.Attributes.unit = "kg/m";
        % VERTICAL TAIL ROOT CHORD 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.value = Aircraft.Geometry.Vertical.croot.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.Attributes.unit = "m";
        % QUARTER CHORD DISTANCE: C_ROOT_VT / 4
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_quarter.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.value/4;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_quarter.Attributes.unit = "m";
        % CHORDWISE KILOGRAMS OF LOAD DISTRIBUTION 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Chordwise_Load_ratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Reference_mass.value/Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Chordwise_Load_ratio.Attributes.unit = "kg/m";

        %% CS - VLA 441 CASE (a)(1)
        % A VECTOR CALLED DR
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.dr.value = [0.0 20.0]';
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.dr.Attributes.unit = "Unknown";

%         % LATERAL FORCE COEFFICIENT AT TWO BETA - FROM CFD
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_vector.value = [0.00003615686 0.01291153]';
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_vector.Attributes.unit = "Non dimensional";

%         % LATERAL FORCE COEFFICIENT WHEN BETA = 0 AND DR = 0
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_0.value = 0.000036156860;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_0.Attributes.unit = "Non dimensional";

%         % LATERAL FORCE GRADIENT C_Y / dr
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CYdr.value = 0.000644;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CYdr.Attributes.unit = "1/deg";

        % MAXIMUM RUDDER DEFLECTION 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.value = 30.0;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.Attributes.unit = "deg";

        % LATERAL FORCE COEFFICIENT WHEN THE RUDDER IS AT MAXIMUM DEFLECTION:
        %              d CY
        % CY = CY_0 + ------- * max(delta_rudder)
        %               d r
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_0.value + Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CYdr.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY.Attributes.unit = "Non dimensional";

        % CALCULATION OF: CY / S_ratio
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY.value/Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.Attributes.unit = "Non dimensional";

        % CALCULATION OF LATERAL FORCE

        % DYNAMIC PRESSURE AT VA
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.Attributes.unit = "Pa";

        % LATERAL FORCE 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.value = n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.Attributes.unit = "N";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.Attributes.cs = " 441(a)(1) ";

        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_decanewton.value = (1e-1)*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_decanewton.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_decanewton.Attributes.cs = " 441(a)(1) ";

        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_kilograms.value = (1/(Aircraft.Constants.g.value))*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_kilograms.Attributes.unit = "kg";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force_kilograms.Attributes.cs = " 441(a)(1) ";

        % THE TOTAL LATERAL FORCE MUST BE DIVIDED BY TWO, TO HAVE THE AIRLOADS
        % ACTING ON A SINGLE VERTICAL FIN
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.value / n_fin;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin.Attributes.unit = "N"; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin.Attributes.cs = " 441(a)(1) ";

        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.value = (1e-1)*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.value / n_fin;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.Attributes.unit = "daN"; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.Attributes.cs = " 441(a)(1) ";

        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_kilograms.value = (1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Lateral_force.value / n_fin;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_kilograms.Attributes.unit = "kg"; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_kilograms.Attributes.cs = " 441(a)(1) ";

        % DISPLAYING RESULTS 
        disp( " ")
        disp(" ++++ CS - VLA 441 - CASE (a)(1) ++++ ");
        disp(" ------------------------------------ ");

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.value];
        disp(" ++++++++++ Vertical tail loads - CASE (a)(1) [daN] ++++++++++ ")
        format = ' %6.6f\n';
        label  = ' Y\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

        %% CS - VLA 441 CASE (a)(2)

        % DYNAMIC PRESSURE AT VA
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.Attributes.unit = "Pa";

        % CALCULATION OF (1.3)*(YAW_ANGLE)
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.yaw_angle.value = (1.3)*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.yaw_angle.Attributes.unit = "deg";

%         % LATERAL FORCE COEFFICIENT 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_VTP.value = 0.0245;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_VTP.Attributes.unit = "Non dimensional";

        % LATERAL FORCE COEFFICIENT 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_VTP.value/Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.Attributes.unit = "Non dimensional";

        % LATERAL FORCE - NEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.value = n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.Attributes.unit = "N";

        % LATERAL FORCE - DECANEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.value = (1e-1)*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.Attributes.unit = "daN";

        % LATERAL FORCE - KILOGRAMS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.value = (1/(Aircraft.Constants.g.value))*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.Attributes.unit = "N";

        % -------------------------------------------------------------------------
        % LATERAL FORCE ON A SINGLE FIN - NEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.value = (0.5)*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.Attributes.unit = "N";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force.Attributes.cs = " 441(a)(2) ";

        % LATERAL FORCE ON A SINGLE FIN - DECANEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.value = (0.5)*(1e-1)*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.Attributes.cs = " 441(a)(2) ";

        % LATERAL FORCE ON A SINGLE FIN - KILOGRAMS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.value = (0.5)*(1/(Aircraft.Constants.g.value))*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.Attributes.unit = "N";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_kilograms.Attributes.cs = " 441(a)(2) ";

        % DISPLAYING RESULTS 
        disp( " ")
        disp(" ++++ CS - VLA 441 - CASE (a)(2) ++++ ");
        disp(" ------------------------------------ ");

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.value];
        disp(" ++++++++++ Vertical tail loads - CASE (a)(2) [daN] ++++++++++ ")
        format = ' %6.6f\n';
        label  = ' Y\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

        %% CS - VLA 441 CASE (a)(3)

        % DYNAMIC PRESSURE AT VA
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.Attributes.unit = "Pa";

        % CALCULATION OF (1.3)*(YAW_ANGLE)
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.yaw_angle.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.yaw_angle.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.yaw_angle.Attributes.unit = "deg";

%         % LATERAL FORCE COEFFICIENT 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_VTP.value = 0.0233;
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_VTP.Attributes.unit = "Non dimensional";

        % LATERAL FORCE COEFFICIENT 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_VTP.value/Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_ratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.Attributes.unit = "Non dimensional";

        % LATERAL FORCE - NEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.value = n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.Attributes.unit = "N";

        % LATERAL FORCE - DECANEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value = (1e-1)*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.Attributes.unit = "daN";

        % LATERAL FORCE - KILOGRAMS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.value = (1/(Aircraft.Constants.g.value))*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.Attributes.unit = "N";

        % -------------------------------------------------------------------------
        % LATERAL FORCE ON A SINGLE FIN - NEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.value = (0.5)*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.Attributes.unit = "N";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force.Attributes.cs = " 441(a)(3) ";

        % LATERAL FORCE ON A SINGLE FIN - DECANEWTON
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value = (0.5)*(1e-1)*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.Attributes.unit = "daN";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.Attributes.cs = " 441(a)(3) ";

        % LATERAL FORCE ON A SINGLE FIN - KILOGRAMS
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.value = (0.5)*(1/(Aircraft.Constants.g.value))*n_fin*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.qA.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.CY_dividedby_Sratio.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.Attributes.unit = "N";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_kilograms.Attributes.cs = " 441(a)(3) ";

        % DISPLAYING RESULTS 
        disp( " ")
        disp(" ++++ CS - VLA 441 - CASE (a)(3) ++++ ");
        disp(" ------------------------------------ ");

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value];
        disp(" ++++++++++ Vertical tail loads - CASE (a)(3) [daN] ++++++++++ ")
        format = ' %f\n';
        label  = ' Y\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

        %% DEFINING CRITICAL LOADS FOR THE VERTICAL FINS
        tl_0 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Y_on_single_fin_decanewton.value;
        tl_1 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_2.Lateral_force_decanewton.value;
        tl_2 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_3.Lateral_force_decanewton.value;
        if (abs(tl_0)>abs(tl_1)) && (abs(tl_0)>abs(tl_2)) 
            disp(" ++++++++++ Vertical tail loads - MOST CRITICAL LOADS [daN] ++++++++++ ");
            disp(" CRITICAL LOADS PER ---> CS - VLA 441 (a)(1) "); 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.value = tl_0;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.unit = "daN"; 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.flag = " CS - VLA 441 (a)(1) ";
        elseif (abs(tl_1)>abs(tl_0)) && (abs(tl_1)>abs(tl_0))
            disp(" ++++++++++ Vertical tail loads - MOST CRITICAL LOADS [daN] ++++++++++ ");
            disp(" CRITICAL LOADS PER ---> CS - VLA 441 (a)(2) "); 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.value = tl_1;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.unit = "daN"; 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.flag = " CS - VLA 441 (a)(2) ";
        elseif (abs(tl_2)>abs(tl_1)) && (abs(tl_2)>abs(tl_0))    
            disp(" ++++++++++ Vertical tail loads - MOST CRITICAL LOADS [daN] ++++++++++ ");
            disp(" CRITICAL LOADS PER ---> CS - VLA 441 (a)(3) "); 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.value = tl_2;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.unit = "daN"; 
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Attributes.flag = " CS - VLA 441 (a)(3) ";
        end
        
%% AMC 441 MANOEUVRING MOMENT 
        %                            rho0
        %  M_rudder = (0.3) * S_ht * ---- * beta * V^2 * b_ht 
        %                              2
        eta  = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Case_a_1.Maximum_delta_rudder.value;
        beta = dLdeta * f_eta * eta;  
        % M RUDDER 
        M_rudder_Manoeuvring = 0.3 * S_ht * rho0 * 0.5 * beta * VA^2 * b_ht; 
        disp(" ++++++++++ Vertical tail loads - MOMENT [daN] ++++++++++ ");
        disp(" CRITICAL MANOEUVRING MOMENT PER ---> CS - VLA AMC 441 "); 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Moment.Manoeuvring.value = M_rudder_Manoeuvring;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Moment.Manoeuvring.Attributes.unit = " N*m "; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Moment.Manoeuvring.Attributes.flag = " CS - VLA AMC 441 ";
        
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Critical_load.Moment.Manoeuvring.value];
        format = ' %f\n';
        label  = ' M_rudder\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

        %% CS - VLA 443 GUST LOADS 
        %  
        %  (a) Vertical tail surfaces must be designed to whithstand, in
        %      unaccelerated flight at speed V = VC, lateral gusts of the values
        %      prescribed for VC in CS - VLA 333 (c).
        %
        %  (b) In the absence of a more rational analysis, the gust load must be
        %      computed as follows 
        %
        %             K_gt * U_de * V * a_vt * S_vt
        %      L_vt = -----------------------------
        %                          16.3
        %
        %      where 
        %
        %      U_de  = derived gust velocities (m/s);
        % 
        %      L_vt  = vertical tail loads (daN);
        %     
        %               (0.88) * mu_gt 
        %      K_gt  = --------------- = gust alleviation factor;
        %                5.3 + mu_gt
        %
        %                          2 * M                   K^2
        %      mu_gt = ------------------------------- * ------- = lat. mass ratio;
        %              rho * c_bar_t * g * a_vt * S_vt   (l_t^2)
        %              
        %              where 
        %                
        %              M    = aeroplane mass (kg); 
        %              rho  = air density (kg/m^3); 
        %              S_vt = area of vertical tail (m^2);
        %              l_t  = distance from aeroplane c.g. to lift centre of
        %                     vertical surface (m); 
        %              a_vt = lift curve slope of vertical tail (1/rad); 
        %              V    = aeroplane equivalent speed (m/s); 
        %              K    = radius of gyration in yaw (m);
        %              g    = acceleration due to gravity (m/s^2).
        %
        %  (c) The average loading in figure B5 and the distribution in figure B8
        %      of Appendix B may be used. 

        %% AMC VLA 443 GUST LOADS 
        %  -----------------------------
        %  INTERPRETATIVE MATERIAL AND ACCEPTABLE MEANS OF COMPLIANCE 
        %
        %  For aeroplanes where the horizontal tail is supported by the vertical
        %  tail, the tail surfaces and their supporting structure including the
        %  rear portion of the fuselage should be designed to whithstand the
        %  prescribed loadings on the vertical tail and the roll - moments induced
        %  by the horizontal tail acting in the same direction. 
        %
        %  For T - Tails in the abscence of a more rational analysis, the rolling
        %  moment induced by gust load may be computed as follow 
        % 
        %                             rho0
        %  M_rudder = (0.3) * S_ht * ----- * V * U * b_ht * K_gf
        %                              2
        %    
        %  where 
        %
        %  M_rudder = induced roll - moment at horizontal tail (N * m);
        %  K_gf     = gust factor = 1.2; 
        %  b_ht     = span of horizontal tail (m); 
        %  S_ht     = area of horizontal tail (m^2); 
        %  rho0     = density of air at sea level (kg/m^3);
        %  V        = airspeed of flight (m/s);
        %  U        = gust speed (m/s).

        %% CS - VLA 445 OUTBOARD FINS 
        %
        %  (a) If outboard fins are on the horizonta tail surface, the tail
        %      surfaces must be designed for the maximum horizontal surface load in
        %      combination with the corresponding loads induced on the vertical
        %      surfaces by endplates effects. These induced effects need not be
        %      combined with other vertical surface loads. 
        %
        %  (b) If outboard fins extend above and below the horizontal surface, the 
        %      critical vertical surface loading, the load per unit area as
        %      determined under CS - VLA 441 and 443, must be applied to 
        %      (1) the part of the vertical surfaces above the horizontal surface
        %          with 80% of that loading applied to the part below the 
        %          horizontal surface; and 
        %      (2) the part of the vertical surfaces below the horizontal surface
        %          with 80% of that loading applied to the part above the 
        %          horizontal surface. 
        %
        %  (c) The endplate effects of outboard fins must be taken into account in
        %      applying the yawing conditions of CS - VLA 441 and 443 to the
        %      vertical surfaces in sub - paragraph (b) of this paragraph.

        %% GUST LOADS CALCULATION 

        % GUST LOAD AT VC
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVC.value = 15.24;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVC.Attributes.unit = "m/s";

        % GUST LOAD AT VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVD.value = 7.62;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVD.Attributes.unit = "m/s";

        % VERTICAL TAIL LIFT SLOPE COEFFICIENT
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value = Aircraft.Certification.Aerodynamic_data.Vertical.a_vt.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.Attributes.unit = "1/rad";

        % VERTICAL TAIL SURFACES
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.Attributes.unit = "m^2";

        % ATMOSPHERE PROPERTIES LOCALLY INSTANTIATED 
        % h = 3000; % [m]
        h = Aircraft.Certification.ISA_Condition.Theoretical_ceiling.Altitude.value;
        [T, acc, p, rho] = atmosisa(h);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.value = rho;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.Attributes.unit = "kg/m^3";
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.Attributes.height = "3500 m";

%         % RADIUS OF GYRATION 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.value = 0.3; 
%         Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.Attributes.unit = "m"; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.Attributes.quantity = "Radius of gyration";

        % AIRSPEED - LOAD AT VC
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VC.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VC.Attributes.unit = "m/s";

        % AIRSPEED - LOAD AT VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VD.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VD.Attributes.unit = "m/s";

        % ++++++++++++++++
        % GUST LOADS AT VC 
        % ++++++++++++++++

        % LATERAL MASS RATIO - VC
        mu_g = @(M, rho, mac_vt, g, a_vt, S_vt, K, l_vt) ((2 * M)/(rho * mac_vt * g * a_vt * S_vt))*(K/l_vt)^2;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.mu_g.value = mu_g(Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.value, ...
                                                                                                Aircraft.Geometry.Vertical.MAC.value, ... 
                                                                                                Aircraft.Constants.g.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.value, ...
                                                                                                Aircraft.Geometry.Vertical.l_vt.value);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.mu_g.Attributes.cs = " 443 ";
        % GUST ALLEVIATION FACTOR - VC
        k_g = @(mu_g) ((0.88)*mu_g)/(5.3 + mu_g);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.K_g.value = k_g(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.mu_g.value);

        % GUST LOADS - VC
        Gust_loads = @(K_gt, U_de, V, a_vt, S_vt) (K_gt * U_de * V * a_vt * S_vt)/(16.3);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.value = Gust_loads(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.K_g.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVC.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VC.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.Attributes.unit = "daN";    
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.Attributes.cs = " 445 ";

        % ++++++++++++++++
        % GUST LOADS AT VD 
        % ++++++++++++++++

        % LATERAL MASS RATIO - VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.mu_g.value = mu_g(Aircraft.Weight.I_Level.W_maxTakeOff.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.rho.value, ...
                                                                                                Aircraft.Geometry.Vertical.MAC.value, ... 
                                                                                                Aircraft.Constants.g.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value, ...
                                                                                                Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.K.value, ...
                                                                                                Aircraft.Geometry.Vertical.l_vt.value);
        % GUST ALLEVIATION FACTOR - VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.K_g.value = k_g(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.mu_g.value);

        % GUST LOADS - VD
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.value = Gust_loads(Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.K_g.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVD.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.VD.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.a_vt.value, ...
                                                                                                                   Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.S_vt.value);
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.Attributes.unit = "daN"; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.Attributes.cs = " 445 ";

        % DISPLAYING RESULTS 
        disp( " ")
        disp(" ++++ CS - VLA 445 ++++ ");
        disp(" ---------------------- ");

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.value, ...
                     Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.value];
        disp(" ++++++++++ Vertical tail - GUST LOADS [daN] ++++++++++ ")
        format = ' %6.6f         %6.6f\n';
        label  = ' GUST - VC        GUST - VD\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

        % CRITICAL GUST LOADS 
        gl_0 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vc.Gust_loads_VC.value;
        gl_1 = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.vd.Gust_loads_VD.value;

        if (abs(gl_0)>abs(gl_1))
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.value = gl_0;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.flag = " at VC";
        elseif (abs(gl_1)>abs(gl_0))
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.value = gl_1;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.flag = " at VD";
        end

        % Horizontal tail loads increments
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.value];
        disp(" ++++++++++ Vertical tail - CRITICAL GUST LOADS [daN] ++++++++++ ")
        disp( Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Attributes.flag );
        format = ' %6.6f\n';
        label  = ' CRITICAL GUST LOADS\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")       

%% AMC 443 GUST MOMENT ON THE VERTICAL TAIL         
        %  -----------------------------
        %                             rho0
        %  M_rudder = (0.3) * S_ht * ----- * V * U * b_ht * K_gf
        %                              2
        %    
        %  where 
        %
        %  M_rudder = induced roll - moment at horizontal tail (N * m);
        %  K_gf     = gust factor = 1.2; 
        %  b_ht     = span of horizontal tail (m); 
        %  S_ht     = area of horizontal tail (m^2); 
        %  rho0     = density of air at sea level (kg/m^3);
        %  V        = airspeed of flight (m/s);
        %  U        = gust speed (m/s).
        % GUST LOAD AT VC
        U_atVC = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVC.value;

        % GUST LOAD AT VD
        U_atVD = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.U_atVD.value;
        
        % VERTICAL TAIL GUST FACTOR
        K_gf_vert = 1.2; 
        
        % MOMENT AT VC 
        M_rudder_gust_VC = 0.3 * S_ht * rho0 * 0.5 * VC * U_atVC * b_ht * K_gf_vert;
        disp(" ++++++++++ Vertical tail loads - MOMENT AT VC [daN] ++++++++++ ");
        disp(" CRITICAL GUST MOMENT PER ---> CS - VLA AMC 441 "); 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Moment.Gust.VC.value = M_rudder_gust_VC;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Moment.Gust.VC.Attributes.unit = " N*m "; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Moment.Gust.VC.Attributes.flag = " CS - VLA AMC 443 at VC ";
        % MOMENT AT VD 
        M_rudder_gust_VD = 0.3 * S_ht * rho0 * 0.5 * VD * U_atVD * b_ht * K_gf_vert;
        disp(" ++++++++++ Vertical tail loads - MOMENT AT VD [daN] ++++++++++ ");
        disp(" CRITICAL GUST MOMENT PER ---> CS - VLA AMC 441 "); 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Moment.Gust.VD.value = M_rudder_gust_VD;
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Moment.Gust.VD.Attributes.unit = " N*m "; 
        Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Moment.Gust.VD.Attributes.flag = " CS - VLA AMC 443 at VD ";
        
        % ASSIGNING CRITICAL LOADS 
        % CRITICAL GUST MOMENT 
        mgl_0 = M_rudder_gust_VC;
        mgl_1 = M_rudder_gust_VD;

        if (abs(mgl_0)>abs(mgl_1))
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Moment.value = mgl_0;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Moment.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Moment.Attributes.flag = " at VC";
        elseif (abs(mgl_1)>abs(mgl_0))
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Moment.value = mgl_1;
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Moment.Attributes.unit = "daN";
            Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Moment.Attributes.flag = " at VD";
        end
        % DISPLAYING CRITICAL GUST MOMENT
        Increment = [Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Moment.value];
        disp(" ++++++++++ Vertical tail - CRITICAL GUST LOADS [daN] ++++++++++ ")
        disp( Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.GustLoads.Critical_gustloads.Moment.Attributes.flag );
        format = ' %6.6f\n';
        label  = ' CRITICAL GUST MOMENT\n';
        fprintf(label);
        fprintf(format, Increment.');
        disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")    

end