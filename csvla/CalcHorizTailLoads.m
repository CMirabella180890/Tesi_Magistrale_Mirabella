%% CalcHorizTailLoads
% Script to evaluate unsymmetrical load conditions associated with elevator
% deflection.
% =========================================================================
%   DESCRIPTION
%   It is useful to remember the following airworthiness rules: 
%   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%   CS-VLA 421 Balancing loads 
%
%   (a) A horizontal tail balancing load is a load necessary to maintain
%       equilibrium in any specified flight condition with no pitching
%       acceleration.
%
%   (b) Horizontal tail surfaces must be designed for the balancing loads
%       occuring at any point on the limit manoeuvring envelope and in the
%       flap conditions specified in CS - VLA 345. The distribution in
%       figure B6 of Appendix B may be used.
%       
%   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%   CS-VLA 423 Manoeuvring loads
%   Each horizontal tail surface must be designed for manoeuvring loads
%   imposed by one of the following conditions (a) + (b), or (c), or (d):
%
%   (a) A sudden deflection of the elevator control at VA, to 
%       (1) the maximum updward deflection; and 
%       (2) the maximum downward deflection, as limited by the control
%           stops or pilot effort, whichever is critical. 
%       The average loading of B11 of Appendix B and the distribution in
%       figure B7 of Appendix B may be used. 
%
%   (b) A sudden upward deflection of the elevator, at speeds abobe VA,
%       followed by a downward deflection of the elevator, resulting in the
%       following combinations of normal and angular acceleration:
%
%       |Condition | Normal acceleration (n) | Angular acceleration (rad/s^2)
%       |          |                         |
%       |Download  |                         |     (20.1)
%       |          | 1.0                     |  + ------- * n_m * (n_m - 1.5)
%       |          |                         |       V
%       |-------------------------------------------------------------------
%       |          |                         |
%       |Upload    |                         |     (20.1)
%       |          | n_m                     |  - ------- * n_m * (n_m - 1.5)
%       |          |                         |       V
%       |-------------------------------------------------------------------
%       
%       where
%       (1) n_m = positive limit manoeuvring factor used in the design of
%                 the aeroplane; and
%       (2) V   = initial speed in m/s.
%
%        The conditions in this paragraph involve loads corresponding to
%        the loads that may occur in "checked manoeuvre": 
%
%         CHECKED 
%        MANOEUVRE = A manoeuvre in which the pitching control is suddenly
%                    displaced in one direction and then suddenly moved in
%                    the opposite direction.
%
%        The deflections and timing avoiding exceeding the limit
%        manoeuvring loads factor. The total tail load for both down and up
%        load conditions is the sum of the balancing tail loads at V and
%        the specified value of the normal load factor n, plus the
%        manoeuvring load increment due to the specified value of the
%        normal load factor n, plus the manoeuvring increment due to the 
%        specified value of the normal load factor n, plus the manoeuvring
%        load increment due to the specified value of the angular 
%        acceleration. The manoeuvring load increment in figure B2 of
%        Appendix B and the distributions in figure B7 (for down loads) and
%        in figure B8 (for up loads) of Appendix B may be used. 
%
%   (c) A sudden deflection of the elevator, the following cases must be
%       considered:
%       
%       (i)   Speed VA, maximum upward deflection;
%       (ii)  Speed VA, maximum downward deflection;
%       (iii) Speed VD, one-third maximum upward deflection;
%       (iv)  Speed VD, one-third maximum downward deflection.
%
%       The followind assumptions must be made: 
%       
%       (A) The aeroplane is initially in level flight, and its attitude
%           and air speed do not change.
%       (B) The loads are balanced by inertia forces. 
%
%   (d) A sudden deflection of the elevator such as to cause the normal
%       acceleration to change from an initial value to a final value, the
%       following cases being considered (see figure 1):
%       
%      SPEED | INITIAL CONDITION | FINAL CONDITION | LOAD FACTOR INCREMENT
%      --------------------------------------------------------------------
%       VA   |        A1         |        A        |        n1 - 1.0
%            |        A          |        A1       |        1.0 - n1
%            |        A1         |        G        |        n4 - 1.0
%            |        G          |        A1       |        1.0 - n4       
%      --------------------------------------------------------------------
%       VD   |        D1         |        D        |        n2 - 1.0
%            |        D          |        D1       |        1.0 - n2
%            |        D1         |        E        |        n3 - 1.0 
%            |        E          |        D1       |        1.0 - n3
%      --------------------------------------------------------------------
%
%       For the purpose of this calculation the difference in air speed
%       between VA and the value corresponding to point G on the
%       manoeuvring envelope can be ignored. The following assumptions must
%       be made:
%
%       (1) the aeroplane is initially in level flight, and its attitude
%           and airspeed do not change;
%       (2) the loads are balanced by inertia forces;
%       (3) the aerodynamic tail load increment is given by:
%      
%                               / X_cg    S_ht     a_ht           d epsilon       1                                 \                                                  
%       DeltaP = DeltaN * Mg * <  ---- - ------ * ------ * [ 1 - ----------- ] - --- * 0.5 * rho0 * S_ht * a_ht * lt >
%                               \  lt      S        a              d alpha        M                                 /
%
%           where
%
%           DeltaP    = horizontal tail load increment, positive upwards 
%                       (N)
%           DeltaN    = load factor increment 
%           M         = mass of the aeroplane (kg) 
%           g         = acceleration due to gravity (m/s^2) 
%           X_cg      = longitudinal distance of aeroplane c.g. aft of the
%                       aerodynamic centre of the aeroplane less horizontal
%                       tail (m)
%           S_ht      = horizontal tail area (m^2) 
%           
%           d epsilon 
%           --------- = rate of change of downwash angle with angle of
%            d alpha    attack
% 
%           rho0      = density of air at sea-level (kg/m^3) 
%           lt        = tail arm (m)
%           S         = wing area (m^2)
%           a         = slope of wing lift curve per radian
%   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%   CS - VLA 425 Gust loads 
%   
%   (a) Each horizontal tail surface must be designed for loads resulting
%       from 
%       (1) gust velocities specified in CS - VLA 333 (c) with flaps
%           retracted; and 
%       (2) positive and negative gust of 7.62 m/s nominal intensity at VF
%           corresponding to the flight conditions specified in CS - VLA
%           345 (a) (2).
%                                                                                 
%   (b) The average loading in figure B3 and the distribution of figure B8 
%       may be used to determine the incremental gust loads for the
%       requirements of subparagraph (a) applied as both up and down
%       increments for subparagraph (c).
%
%   (c) When determining the total load on the horizontal tail for the
%       conditions specified in subparagraph (a) of this paragraph, the
%       initial balancing tail loads for steady unaccelerated flight at the
%       pertinent design speeds VF, VC and VD must first be determined. The
%       incremental tail load resulting from the gusts must be added to the
%       initial balancing tail load to obtain the total tail load. 
%
%   (d) In the absence of a more rational analysis, the incremental tail
%       load due to the gust, must be computed as follows:
%       
%                     Kg * U_de * V * a_ht * S_ht          d epsilon
%       Delta L_ht = ----------------------------- * [1 - -----------]
%                              (16) * (3)                   d alpha
%
%       where 
%       Delta L_ht      = incremental horizontal tail load (daN);
%       Kg              = gust alleviation factor defined in CS - VLA 341;
%       U_de            = derived gust velocity (m/s);
%       V               = aeroplane equivalent speed (m/s);
%       a_ht            = slope of horizontal tail lift curve per radian; 
%       S_ht            = area of horizontal tail (m^2); 
%       
%            d epsilon 
%       1 - ----------- = downwash factor.
%             d alpha
%
% =========================================================================
%% CS - VLA 423 - METHOD A - MANOEUVRING AIRSPEED VA 

% SWITCH TO SELECT CORRECTLY THE HORIZ. CONTROL EFFICIENCY TAU
if Aircraft.Certification.Aerodynamic_data.Horizontal.tau.Attributes.flag == "Conventional"
    Aircraft.Certification.Aerodynamic_data.Horizontal.tau.value = 0.50;
elseif Aircraft.Certification.Aerodynamic_data.Horizontal.tau.Attributes.flag == "Full movable"
    Aircraft.Certification.Aerodynamic_data.Horizontal.tau.value = 1.0;
elseif Aircraft.Certification.Aerodynamic_data.Horizontal.tau.Attributes.flag == "Custom"
    prmpt = "Enter Horiz. control efficiency --> tau: ";
    Aircraft.Certification.Aerodynamic_data.Horizontal.tau.value = input(prmpt);
end

% METHOD (a) OF THE CS VLA 423
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.l_tail.value = Aircraft.Geometry.Horizontal.l.value; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.l_tail.Attributes.unit = "meters";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_tail.value = Aircraft.Geometry.Horizontal.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_tail.Attributes.unit = "squared meters";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.qA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.qA.Attributes.unit = "Pa";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.a_tail_rad.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.a_tail_rad.Attributes.unit = "1/radians";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.a_tail_grad.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value*(pi/180);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.a_tail_grad.Attributes.unit = "1/degrees";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.CL_delta_elevator.value = Aircraft.Certification.Aerodynamic_data.Horizontal.CL_delta_elevator.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.CL_delta_elevator.Attributes.unit = "1/radians"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.CL_delta_elevator_grad.value = Aircraft.Certification.Aerodynamic_data.Horizontal.CL_delta_elevator.value*(pi/180);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.CL_delta_elevator_grad.Attributes.unit = "1/degrees";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_step.value = 0.006;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_step.Attributes.unit = "seconds";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_interval_num.value = 50;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_interval_num.Attributes.unit = "Pure number";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value = linspace(0, ...
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_interval_num.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_step.value, ...
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_interval_num.value)';
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.Attributes.unit = "seconds";

% AMC 23.423 ADVISED TOTAL DEFLECTION TIME INTERVAL 
% BE CAREFUL: AMC 23.423 Suggest the use of various total deflection time
% interval for different aircraft categories (Normal, Utility, Commuter,
% Acrobatic, ...). Pleas, take care of the definition of the time interval.
if (Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes.flag1 == "Aerobatic") && (Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes.flag2 == "Stick")
    Aircraft.Geometry.Horizontal.Movable.total_deflection_time.value = 0.1;
    Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes = "seconds";
elseif (Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes.flag1 == "Aerobatic") && (Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes.flag2 == "Wheel")
    Aircraft.Geometry.Horizontal.Movable.total_deflection_time.value = 0.2;
    Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes = "seconds";
elseif (Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes.flag1 == "Normal") | (Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes.flag1 == "Utility") | (Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes.flag1 == "Normal")
    if Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes.flag2 == "Stick"
        Aircraft.Geometry.Horizontal.Movable.total_deflection_time.value = 0.2;
        Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes = "seconds";
    elseif Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes.flag2 == "Wheel"
        Aircraft.Geometry.Horizontal.Movable.total_deflection_time.value = 0.3;
        Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes = "seconds";
    end
elseif Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes.flag1 == "Normal"
    prmpt = "Enter total defl. time interval --> t_total_defl_time: ";
    Aircraft.Geometry.Horizontal.Movable.total_deflection_time.value = input(prmpt);
    Aircraft.Geometry.Horizontal.Movable.total_deflection_time.Attributes = "seconds";
end

Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.total_time_of_deflection.value = 0.3;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.total_time_of_deflection.Attributes.unit = "sec";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.damping_factor.value = 0.3; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.damping_factor.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.IY.value = Aircraft.Weight.I_Level.IY.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.IY.Attributes.unit = "kg * m^2"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dvVA.value = 0.01*(180/pi); 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dvVA.Attributes.unit = "degree";

% RATIO L_TAIL / M.A.C. 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.L_ratio.value =  Aircraft.Geometry.Horizontal.l.value/Aircraft.Geometry.Wing.mac.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.L_ratio.Attributes.unit = "Non dimensional";

% RATIO S_TAIL / S_WING 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_ratio.Attributes.unit = "Non dimensional"; 

% HORIZONTAL TAIL VOLUME RATIO 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.Horizontal_Tail_Volume_Ratio.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.L_ratio.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_ratio.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.Horizontal_Tail_Volume_Ratio.Attributes.unit = "Non dimensional";

% SOLVING THE EQUATION OF MOTION 
% 
%     d^2 theta   (q * S_tail * a_tail * d)   /           delta_v     \
%     --------- = ------------------------- *| omega*dt - ------- * DF|
%       dt^2                 IY               \             VA        /
%

Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega.value = (Aircraft.Geometry.Horizontal.Movable.max_deflection.value*Aircraft.Certification.Aerodynamic_data.Horizontal.tau.value)*(1/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.total_time_of_deflection.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega.Attributes.unit = "deg/sec";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega_rad.Attributes.unit = "rad/sec";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.A0.value = (1/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.IY.value)*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.a_tail_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.qA.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.l_tail.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_tail.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.A0.Attributes.unit = "1/(rad*m*s^2)";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz.Attributes.unit = "radians";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz_deg.Attributes.unit = "degrees";

% SOLVING THE DIFFERENTIAL EQUATION 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.d2thetadt2.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.d2thetadt2.Attributes.unit = "rad/sec^2"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dthetadt.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dthetadt.Attributes.unit = "rad/sec"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz.Attributes.unit = "rad"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_theta.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_theta.Attributes.unit = "rad"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_v.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_v.Attributes.unit = "m/s"; 

 for i = 2:length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value)
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.d2thetadt2.value(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.A0.value*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz.value(i) - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_theta.value(i-1));
    disp(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.d2thetadt2.value(i))
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dthetadt.value(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dthetadt.value(i-1) + 0.5*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.d2thetadt2.value(i-1)+Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.d2thetadt2.value(i))*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_step.value);
    disp(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dthetadt.value(i));
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_v.value(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dthetadt.value(i)*Aircraft.Geometry.Horizontal.l.value;
    disp(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_v.value(i));
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_theta.value(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_v.value(i)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.damping_factor.value/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value);
    disp(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_theta.value(i));
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz.value(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz.value(i)  - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_theta.value(i);
    disp(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz.value(i));
 end

Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz_grad.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz_grad.Attributes.unit = "deg";

% LIMIT HORIZONTAL TAIL LOAD 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.DeltaLimitLTail.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz_grad.value(end)*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.a_tail_grad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_tail.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.qA.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.DeltaLimitLTail.Attributes.unit = "daN";

disp(" ++++ METHOD CS - VLA 423 (a) ++++");
disp(" ++++ UNCHECKED ++++");

% Horizontal tail loads increments
Increment = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.DeltaLimitLTail.value];
disp(" ++++++++++ Critical Horizontal Tail loads increments [daN] ++++++++++ ")
format = '%f          \n';
label  = 'VA                 \n';
fprintf(label);
fprintf(format, Increment.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

% A0 = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.qA.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_tail.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.a_tail_grad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.l_tail.value)/(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.IY.value);
% C  = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dvVA.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.damping_factor.value;
% B  =(Aircraft.Constants.g.value/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value - 1);
% t_span = [ 0.0, 50.0 ];
% f = @(t, y) [ B; y(2) - A0*(y(1)*t - C)]; % define function f(t,y)
% t0 = 0;  
% y0 = [ 17.0/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega.value; 17.0/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega.value ];                         % initial condition with vector y0
% % y0 = [ 0.0; 0.0 ]; 
% % func_hand = @(t, theta) [ theta(2) A0*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega.value*t - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dvVA.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.damping_factor.value) ]; 
% [t, theta_dot] = ode45(f, t_span, y0);
% theta0         = [ 0.0, 0.0 ]; 
% A_theta        = (theta_dot(2,1) - theta_dot(1,1))/(t(2) - t(1));
% f_theta        = @(t, theta_dot) - theta_dot + A_theta*t;
% [t_t, theta]   = ode45(f_theta, t_span, 0.0);
% % -------------------------------------------------------------------------
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.theta_dot_plot.value = figure();
% hold on 
% grid on; grid minor; 
% plot(t, theta_dot(:,1), '-k.', 'linewidth', 1.5);
% title("Angular speed vs Time", "interpreter", "latex");
% xlabel("Time - $t$ $[s]$", "interpreter", "latex");
% ylabel("Angular velocity - $\dot{\theta}$ $[deg/s]$", "interpreter", "latex");
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.theta_dot_plot.value, 'ThetaDot.pdf', 'ContentType', 'vector')
% % Saving figures inside correct folder
% fprintf('Saving in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile ThetaDot.pdf Output
% % -------------------------------------------------------------------------
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.theta_plot0.value = figure();
% hold on 
% grid on; grid minor; 
% plot(t_t, theta(:,1), '-r.', 'linewidth', 1.5);
% title("Pitch angle vs Time", "interpreter", "latex");
% xlabel("Time - $t$ $[s]$", "interpreter", "latex");
% ylabel("Pitch angle - $\theta$ $[deg]$", "interpreter", "latex");
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.theta_plot0.value, 'Theta.pdf', 'ContentType', 'vector')
% % Saving figures inside correct folder
% fprintf('Saving in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile Theta.pdf Output
% % -------------------------------------------------------------------------
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.theta_plot1.value = figure();
% hold on 
% grid on; grid minor; 
% plot(t_t, theta(:,1), '-r.', 'linewidth', 1.5);
% title("Pitch angle vs Time", "interpreter", "latex");
% xlabel("Time - $t$ $[s]$", "interpreter", "latex");
% ylabel("Pitch angle - $\theta$ $[deg]$", "interpreter", "latex");
% xlim([0.0 10.0]);
% ylim([0.0 10.0]);

% SAVING FIGURE INSIDE OUTPUT FOLDER
% exportgraphics(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.theta_plot1.value, 'ThetaZoom.pdf', 'ContentType', 'vector')
% % Saving figures inside correct folder
% fprintf('Saving in: ');
% fprintf('\n'); 
% fprintf('%s\n', SaveFolder);
% % Moving file inside correct folder
% movefile ThetaZoom.pdf Output

%% CS - VLA 423 - METHOD B - MANOEUVRING AIRSPEED VA 
% % NOSE UP PITCHING 
% 
% % L TAIL 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.l_tail.value = Aircraft.Geometry.Horizontal.l.value; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.l_tail.Attributes.unit = "meters";
% 
% % MOMENT OF INERTIA 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.IY.value = Aircraft.Weight.I_Level.IY.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.IY.Attributes.unit = "kg * m^2"; 
% 
% % V = VA
% % ANGULAR ACCELERATION PER CS VLA 423
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VA.pitch_up_acceleration.value = (20.1/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value - 1.5);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VA.pitch_up_acceleration.Attributes.unit = "rad/sec^2";
% 
% % PITCHING MOMENT AT VA
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VA.MA.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.IY.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VA.pitch_up_acceleration.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VA.MA.Attributes.unit = "N*m";
% 
% % TAIL LOAD
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VA.tail_load.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VA.MA.value/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.l_tail.value)*1e-1;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VA.tail_load.Attributes.unit = "daN";
% 
% % TOTAL LOAD 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VA.total_horizontal_load.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTail_A.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VA.tail_load.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VA.total_horizontal_load.Attributes.unit = "daN";
% 
% % V = VC
% % ANGULAR ACCELERATION PER CS VLA 423
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VC.pitch_up_acceleration.value = (20.1/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value - 1.5);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VC.pitch_up_acceleration.Attributes.unit = "rad/sec^2";
% 
% % PITCHING MOMENT AT VC
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VC.MC.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.IY.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VC.pitch_up_acceleration.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VC.MC.Attributes.unit = "N*m";
% 
% % TAIL LOAD
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VC.tail_load.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VC.MC.value/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.l_tail.value)*1e-1;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VC.tail_load.Attributes.unit = "daN";
% 
% % TOTAL LOAD 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VC.total_horizontal_load.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTail_C.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VC.tail_load.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VC.total_horizontal_load.Attributes.unit = "daN";
% 
% % V = VD
% % ANGULAR ACCELERATION PER CS VLA 423
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VD.pitch_up_acceleration.value = (20.1/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value - 1.5);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VD.pitch_up_acceleration.Attributes.unit = "rad/sec^2";
% 
% % PITCHING MOMENT AT VD
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VD.MD.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.IY.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VD.pitch_up_acceleration.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VD.MD.Attributes.unit = "N*m";
% 
% % TAIL LOAD
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VD.tail_load.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VD.MD.value/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.l_tail.value)*1e-1;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VD.tail_load.Attributes.unit = "daN";
% 
% % TOTAL LOAD 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VD.total_horizontal_load.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTail_D.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VD.tail_load.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VD.total_horizontal_load.Attributes.unit = "daN";
% 
% % NOSE DOWN PITCHING 
% 
% % V = VA
% % ANGULAR ACCELERATION PER CS VLA 423
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VA.pitch_down_acceleration.value = -(20.1/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value - 1.5);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VA.pitch_down_acceleration.Attributes.unit = "rad/sec^2";
% 
% % PITCHING MOMENT AT VA
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VA.MA.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.IY.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VA.pitch_down_acceleration.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VA.MA.Attributes.unit = "N*m";
% 
% % TAIL LOAD
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VA.tail_load.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VA.MA.value/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.l_tail.value)*1e-1;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VA.tail_load.Attributes.unit = "daN";
% 
% % TOTAL LOAD 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VA.total_horizontal_load.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTail_A.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VA.tail_load.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VA.total_horizontal_load.Attributes.unit = "daN";
% 
% % V = VC
% % ANGULAR ACCELERATION PER CS VLA 423
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VC.pitch_down_acceleration.value = -(20.1/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value - 1.5);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VC.pitch_down_acceleration.Attributes.unit = "rad/sec^2";
% 
% % PITCHING MOMENT AT VA
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VC.MC.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.IY.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VC.pitch_down_acceleration.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VC.MC.Attributes.unit = "N*m";
% 
% % TAIL LOAD
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VC.tail_load.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VC.MC.value/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.l_tail.value)*1e-1;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VC.tail_load.Attributes.unit = "daN";
% 
% % TOTAL LOAD 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VC.total_horizontal_load.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTail_C.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VC.tail_load.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VC.total_horizontal_load.Attributes.unit = "daN";
% 
% % V = VD
% % ANGULAR ACCELERATION PER CS VLA 423
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VD.pitch_down_acceleration.value = -(20.1/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value - 1.5);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VD.pitch_down_acceleration.Attributes.unit = "rad/sec^2";
% 
% % PITCHING MOMENT AT VA
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VD.MD.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.IY.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VD.pitch_down_acceleration.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VD.MD.Attributes.unit = "N*m";
% 
% % TAIL LOAD
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VD.tail_load.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VD.MD.value/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.l_tail.value)*1e-1;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VD.tail_load.Attributes.unit = "daN";
% 
% % TOTAL LOAD 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VD.total_horizontal_load.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTail_D.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VD.tail_load.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VD.total_horizontal_load.Attributes.unit = "daN";

%% PRINT RESULTS 

% disp(" ++++ METHOD CS - VLA 423 (b) ++++");
% disp(" ++++ PITCH UP ++++");
% 
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VA.tail_load.value, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VC.tail_load.value, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VD.tail_load.value];
% disp(" ++++++++++ Critical Horizontal Tail loads increments [daN] ++++++++++ ")
% format = '%f          %f          %f\n';
% label  = 'VA                 VC                VD\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % Total horizontal tail increment
% Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VA.total_horizontal_load.value, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VC.total_horizontal_load.value, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.VD.total_horizontal_load.value];
% disp(" +++++++++++++++++ Total Horizontal Tail loads [daN] +++++++++++++++++ ")
% format = '%f          %f          %f\n';
% label  = 'VA                 VC                VD\n';
% fprintf(label);
% fprintf(format, Total.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% disp(" ++++ METHOD CS - VLA 423 (b) ++++");
% disp(" ++++ PITCH DOWN ++++");
% 
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VA.tail_load.value, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VC.tail_load.value, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VD.tail_load.value];
% disp(" ++++++++++ Critical Horizontal Tail loads increments [daN] ++++++++++ ")
% format = '%f          %f          %f\n';
% label  = 'VA                 VC                VD\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % Total horizontal tail increment
% Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VA.total_horizontal_load.value, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VC.total_horizontal_load.value, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.VD.total_horizontal_load.value];
% disp(" +++++++++++++++++ Total Horizontal Tail loads [daN] +++++++++++++++++ ")
% format = '%f          %f          %f\n';
% label  = 'VA                 VC                VD\n';
% fprintf(label);
% fprintf(format, Total.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

%% CS - VLA 423 - METHOD B - PITCH UP
% LOAD FACTOR FROM C TO D
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.n_fromCtoD.value = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.value, ...
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value, ...
    length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value))';
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.n_fromCtoD.Attributes.unit = "g's";

% LOAD FACTOR FROM A TO C
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.n_fromAtoC.value = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.value, ...
length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value))';
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.n_fromAtoC.Attributes.unit = "g's";

% FULL LOAD FACTOR 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.full_load_factor_vector.value = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.n_fromAtoC.value; Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.n_fromCtoD.value];
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.full_load_factor_vector.Attributes.unit = "g's";

% UNIT LOAD FACTOR 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.unit_load_factor_vector.value = ones(2*length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value), 1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.unit_load_factor_vector.Attributes.unit = "g's";

% AIRSPEED VECTOR FOR CALCULATIONS
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.airspeed_vector.value = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value, ...
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value, ...
    2*length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value))';
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.airspeed_vector.Attributes.unit = "m/s";

% ANGULAR ACCELERATION CALCULATIONS 
v  = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.airspeed_vector.value;
nm = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.full_load_factor_vector.value;
lm = (20.1/v).*nm.*(nm - 1.5);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.ang_acc.value = lm;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.ang_acc.Attributes.unit = "rad/sec^2";
clear v nm lm

% MOMENT MY
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.MY.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.ang_acc.value), 1);

for i = 1:length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.ang_acc.value)
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.MY.value(i) = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.IY.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.ang_acc.value(i));
end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.MY.Attributes.unit = "N*m";

% DELTA TAIL AIRLOADS
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.delta_tail_airloads.value = (1e-1)*(-Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.MY.value)/(Aircraft.Geometry.Horizontal.l.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.delta_tail_airloads.Attributes.unit = "daN";

% BALANCING TAIL AIRLOADS
index_va = dsearchn(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.balancing_tail_airloads.value = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_unit_load_factor.value(index_va), ...
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_unit_load_factor.value(end), ...
    length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.delta_tail_airloads.value))';
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.balancing_tail_airloads.value = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_positivestall.value; Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromCtoD.value];
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.balancing_tail_airloads.Attributes.unit = "daN";

% TOTAL AIRLOADS ASSOCIATED WITH THE PITCH UP MANOEUVRE
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.total_tail_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.delta_tail_airloads.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.balancing_tail_airloads.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.total_tail_airloads.Attributes.unit = "daN";

% CRITICAL TAIL AIRLOADS VALUE 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.value = 0.0; 
Absolute_value    = abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.total_tail_airloads.value);
max_critical_load = max(Absolute_value);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.total_tail_airloads.value)
    x = (-max_critical_load) - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.total_tail_airloads.value;
    if abs(x) < 1e-2
        max_load_index = i;
        Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.total_tail_airloads.value(max_load_index);
        Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.Attributes.index = max_load_index;
    end
end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.Attributes.unit = "daN";

% Total horizontal tail increment
Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.value];
disp(" +++++++++++++++++ Total Horizontal Tail loads [daN] +++++++++++++++++ ")
format = '%f\n';
label  = 'VA\n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VA - POINT A OF THE FLIGHT ENVELOPE - CASE 1

% Normal load n at Point A of the flight envelope
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.nA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.nA.Attributes.unit = "g's";

% Normal load n at Point A1, indicated in figure 1 of the CS - VLA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.nA1.value = 1.0; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.nA1.Attributes.unit = "g's";

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.DeltaN.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.nA.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.nA1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.DeltaN.Attributes.unit = "g's";

% Evaluation of the product Mg * Delta N
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.DeltaN.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Manoeuvring_weight.Attributes.unit = "N";

% Evaluation of the ratio X_cg/lt 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.X_cg_lt_ratio.Attributes.unit = "Non dimensional";

% Defining the downwash factor
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Downwash_factor.value = 1 - Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Downwash_factor.Attributes.unit = "Non dimensional";

% Defining the ratio S_ht/S - surface_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Surface_ratio.Attributes.unit = "Non dimensional";

% Defining the ratio a_ht/a - Lift_curve_slope_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Downwash_factor.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Third_term.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE ONE +++ 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.parenthesis.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VA - POINT A OF THE FLIGHT ENVELOPE - CASE 2

% Normal load n at Point A of the flight envelope
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.nA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.nA.Attributes.unit = "g's";

% Normal load n at Point A1, indicated in figure 1 of the CS - VLA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.nA1.value = 1.0; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.nA1.Attributes.unit = "g's";

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.DeltaN.value = -Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.nA.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.nA1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.DeltaN.Attributes.unit = "g's";

% Evaluation of the product Mg * Delta N
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.DeltaN.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Manoeuvring_weight.Attributes.unit = "N";

% Evaluation of the ratio X_cg/lt 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.X_cg_lt_ratio.Attributes.unit = "Non dimensional";

% Defining the downwash factor
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Downwash_factor.value = 1 - Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Downwash_factor.Attributes.unit = "Non dimensional";

% Defining the ratio S_ht/S - surface_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Surface_ratio.Attributes.unit = "Non dimensional";

% Defining the ratio a_ht/a - Lift_curve_slope_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Downwash_factor.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Third_term.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE TWO +++ 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.parenthesis.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VA - POINT A OF THE FLIGHT ENVELOPE - CASE 3

% Normal load n at Point A of the flight envelope
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.nG.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.nG.Attributes.unit = "g's";

% Normal load n at Point A1, indicated in figure 1 of the CS - VLA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.nA1.value = 1.0; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.nA1.Attributes.unit = "g's";

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.DeltaN.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.nG.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.nA1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.DeltaN.Attributes.unit = "g's";

% Evaluation of the product Mg * Delta N
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.DeltaN.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Manoeuvring_weight.Attributes.unit = "N";

% Evaluation of the ratio X_cg/lt 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.X_cg_lt_ratio.Attributes.unit = "Non dimensional";

% Defining the downwash factor
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Downwash_factor.value = 1 - Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Downwash_factor.Attributes.unit = "Non dimensional";

% Defining the ratio S_ht/S - surface_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Surface_ratio.Attributes.unit = "Non dimensional";

% Defining the ratio a_ht/a - Lift_curve_slope_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Downwash_factor.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Third_term.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE THREE +++ 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.parenthesis.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VA - POINT A OF THE FLIGHT ENVELOPE - CASE 4

% Normal load n at Point A of the flight envelope
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.nG.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.nG.Attributes.unit = "g's";

% Normal load n at Point A1, indicated in figure 1 of the CS - VLA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.nA1.value = 1.0; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.nA1.Attributes.unit = "g's";

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.DeltaN.value = - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.nG.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.nA1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.DeltaN.Attributes.unit = "g's";

% Evaluation of the product Mg * Delta N
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.DeltaN.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Manoeuvring_weight.Attributes.unit = "N";

% Evaluation of the ratio X_cg/lt 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.X_cg_lt_ratio.Attributes.unit = "Non dimensional";

% Defining the downwash factor
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Downwash_factor.value = 1 - Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Downwash_factor.Attributes.unit = "Non dimensional";

% Defining the ratio S_ht/S - surface_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Surface_ratio.Attributes.unit = "Non dimensional";

% Defining the ratio a_ht/a - Lift_curve_slope_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Downwash_factor.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Third_term.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE FOUR +++ 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.parenthesis.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VD - POINT D OF THE FLIGHT ENVELOPE - CASE 1

% Normal load n at Point A of the flight envelope
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.nA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.nA.Attributes.unit = "g's";

% Normal load n at Point A1, indicated in figure 1 of the CS - VLA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.nA1.value = 1.0; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.nA1.Attributes.unit = "g's";

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.DeltaN.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.nA.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.nA1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.DeltaN.Attributes.unit = "g's";

% Evaluation of the product Mg * Delta N
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.DeltaN.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Manoeuvring_weight.Attributes.unit = "N";

% Evaluation of the ratio X_cg/lt 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.X_cg_lt_ratio.Attributes.unit = "Non dimensional";

% Defining the downwash factor
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Downwash_factor.value = 1 - Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Downwash_factor.Attributes.unit = "Non dimensional";

% Defining the ratio S_ht/S - surface_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Surface_ratio.Attributes.unit = "Non dimensional";

% Defining the ratio a_ht/a - Lift_curve_slope_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Downwash_factor.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Third_term.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE ONE +++ 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.parenthesis.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VD - POINT A OF THE FLIGHT ENVELOPE - CASE 2

% Normal load n at Point A of the flight envelope
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.nA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.nA.Attributes.unit = "g's";

% Normal load n at Point A1, indicated in figure 1 of the CS - VLA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.nA1.value = 1.0; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.nA1.Attributes.unit = "g's";

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.DeltaN.value = -Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.nA.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.nA1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.DeltaN.Attributes.unit = "g's";

% Evaluation of the product Mg * Delta N
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.DeltaN.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Manoeuvring_weight.Attributes.unit = "N";

% Evaluation of the ratio X_cg/lt 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.X_cg_lt_ratio.Attributes.unit = "Non dimensional";

% Defining the downwash factor
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Downwash_factor.value = 1 - Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Downwash_factor.Attributes.unit = "Non dimensional";

% Defining the ratio S_ht/S - surface_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Surface_ratio.Attributes.unit = "Non dimensional";

% Defining the ratio a_ht/a - Lift_curve_slope_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Downwash_factor.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Third_term.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE TWO +++ 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.parenthesis.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VD - POINT E OF THE FLIGHT ENVELOPE - CASE 3

% Normal load n at Point E of the flight envelope
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.nE.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.nE.Attributes.unit = "g's";

% Normal load n at Point D1, indicated in figure 1 of the CS - VLA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.nD1.value = 1.0; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.nD1.Attributes.unit = "g's";

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.DeltaN.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.nE.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.nD1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.DeltaN.Attributes.unit = "g's";

% Evaluation of the product Mg * Delta N
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.DeltaN.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Manoeuvring_weight.Attributes.unit = "N";

% Evaluation of the ratio X_cg/lt 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.X_cg_lt_ratio.Attributes.unit = "Non dimensional";

% Defining the downwash factor
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Downwash_factor.value = 1 - Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Downwash_factor.Attributes.unit = "Non dimensional";

% Defining the ratio S_ht/S - surface_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Surface_ratio.Attributes.unit = "Non dimensional";

% Defining the ratio a_ht/a - Lift_curve_slope_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Downwash_factor.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Third_term.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE THREE +++ 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.parenthesis.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VD - POINT E OF THE FLIGHT ENVELOPE - CASE 4

% Normal load n at Point A of the flight envelope
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.nE.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.nE.Attributes.unit = "g's";

% Normal load n at Point A1, indicated in figure 1 of the CS - VLA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.nD1.value = 1.0; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.nD1.Attributes.unit = "g's";

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.DeltaN.value = - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.nE.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.nD1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.DeltaN.Attributes.unit = "g's";

% Evaluation of the product Mg * Delta N
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.DeltaN.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Manoeuvring_weight.Attributes.unit = "N";

% Evaluation of the ratio X_cg/lt 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.X_cg_lt_ratio.Attributes.unit = "Non dimensional";

% Defining the downwash factor
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Downwash_factor.value = 1 - Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Downwash_factor.Attributes.unit = "Non dimensional";

% Defining the ratio S_ht/S - surface_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Surface_ratio.Attributes.unit = "Non dimensional";

% Defining the ratio a_ht/a - Lift_curve_slope_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Downwash_factor.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Third_term.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE FOUR +++ 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.parenthesis.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";

%% TOTAL AIRLOADS = L_tail + DELTA L_tail

% AIRSPEED VA - CASE 1 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTail_A.value;

% AIRSPEED VA - CASE 2 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTail_S.value;

% AIRSPEED VA - CASE 3 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTail_A.value;

% AIRSPEED VA - CASE 4 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTail_G.value;
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% AIRSPEED VD - CASE 1 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTail_F.value;

% AIRSPEED VD - CASE 2 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTail_E.value;

% AIRSPEED VD - CASE 3 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTail_E.value;

% AIRSPEED VD - CASE 4 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTail_D.value;

%% PRINT RESULTS 

% Horizontal tail loads increments
Increment = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Manoeuvring_Critical_Load_Increment.value, ...
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Manoeuvring_Critical_Load_Increment.value, ...
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Manoeuvring_Critical_Load_Increment.value, ...
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Manoeuvring_Critical_Load_Increment.value; ...
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Manoeuvring_Critical_Load_Increment.value, ...
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Manoeuvring_Critical_Load_Increment.value, ...
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Manoeuvring_Critical_Load_Increment.value, ...
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Manoeuvring_Critical_Load_Increment.value];
disp(" ++++++++++ Critical Horizontal Tail loads increments [daN] ++++++++++ ")
format = '%f          %f          %f          %f\n';
label  = 'Case 1              Case 2             Case 3             Case 4  \n';
fprintf(label);
fprintf(format, Increment.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

% Total horizontal tail increment
Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Total_airloads.value; ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Total_airloads.value];
disp(" +++++++++++++++++ Total Horizontal Tail loads [daN] +++++++++++++++++ ")
format = '%f          %f          %f          %f\n';
label  = 'Case 1              Case 2             Case 3             Case 4  \n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

%% METHOD 
