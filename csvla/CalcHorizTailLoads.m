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
Aircraft.Geometry.Horizontal.Movable.max_deflection.value = 15.0;
Aircraft.Geometry.Horizontal.Movable.max_deflection.Attributes.unit = "degrees";
Aircraft.Geometry.Fuselage.Non_dim_radius_of_gyration.value = 0.34;
Aircraft.Geometry.Fuselage.Non_dim_radius_of_gyration.Attributes.unit = "Non dimensional";
Aircraft.Geometry.Fuselage.Radius_of_gyration.value = Aircraft.Geometry.Fuselage.length.value*Aircraft.Geometry.Fuselage.Non_dim_radius_of_gyration.value*0.5;
Aircraft.Geometry.Fuselage.Radius_of_gyration.Attributes.unit = "meters";
Aircraft.Weight.I_Level.IY.value = (Aircraft.Weight.I_Level.W_maxTakeOff.value)*(Aircraft.Geometry.Fuselage.Radius_of_gyration.value)^2;
Aircraft.Weight.I_Level.IY.Attributes.unit = "kg * m^2"; 
Aircraft.Geometry.Horizontal.tau.value = 0.75; 
Aircraft.Geometry.Horizontal.tau.Attributes.unit = "Non dimensional";

% OPTION (a) OF THE CS VLA 423
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.l_tail.value = Aircraft.Geometry.Horizontal.l.value; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.l_tail.Attributes.unit = "meters";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.S_tail.value = Aircraft.Geometry.Horizontal.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.S_tail.Attributes.unit = "squared meters";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.qA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.qA.Attributes.unit = "Pa";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.a_tail_rad.value = 4.8;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.a_tail_rad.Attributes.unit = "1/radians";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.a_tail_grad.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.a_tail_rad.value*(pi/180);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.a_tail_grad.Attributes.unit = "1/degrees";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.time_step.value = 0.3;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.time_step.Attributes.unit = "seconds";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.omega.value = Aircraft.Geometry.Horizontal.Movable.max_deflection.value/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.time_step.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.omega.Attributes.unit = "deg/sec";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.damping_factor.value = 0.3; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.damping_factor.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.IY.value = Aircraft.Weight.I_Level.IY.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.IY.Attributes.unit = "kg * m^2"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.dvVA.value = 0.01*(180/pi); 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.dvVA.Attributes.unit = "degree";

% SOLVING THE EQUATION OF MOTION 
% 
%     d^2 theta   (q * S_tail * a_tail * d)   /           delta_v     \
%     --------- = ------------------------- *| omega*dt - ------- * DF|
%       dt^2                 IY               \             VA        /
%

A0 = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.qA.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.S_tail.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.a_tail_grad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.l_tail.value)/(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.IY.value);
C  = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.dvVA.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.damping_factor.value;
B  =(Aircraft.Constants.g.value/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value - 1);
t_span = [ 0.0, 50.0 ];
f = @(t, y) [ B; y(2) - A0*(y(1)*t - C)]; % define function f(t,y)
t0 = 0;  
y0 = [ 17.0/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.omega.value; 17.0/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.omega.value ];                         % initial condition with vector y0
% y0 = [ 0.0; 0.0 ]; 
% func_hand = @(t, theta) [ theta(2) A0*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.omega.value*t - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.dvVA.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.damping_factor.value) ]; 
[t, theta_dot] = ode45(f, t_span, y0);
theta0         = [ 0.0, 0.0 ]; 
A_theta        = (theta_dot(2,1) - theta_dot(1,1))/(t(2) - t(1));
f_theta        = @(t, theta_dot) - theta_dot + A_theta*t;
[t_t, theta]   = ode45(f_theta, t_span, 0.0);
% -------------------------------------------------------------------------
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.theta_dot_plot.value = figure();
hold on 
grid on; grid minor; 
plot(t, theta_dot(:,1), '-k.', 'linewidth', 1.5);
title("Angular speed vs Time", "interpreter", "latex");
xlabel("Time - $t$ $[s]$", "interpreter", "latex");
ylabel("Angular velocity - $\dot{\theta}$ $[deg/s]$", "interpreter", "latex");
% -------------------------------------------------------------------------
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.theta_plot0.value = figure();
hold on 
grid on; grid minor; 
plot(t_t, theta(:,1), '-r.', 'linewidth', 1.5);
title("Pitch angle vs Time", "interpreter", "latex");
xlabel("Time - $t$ $[s]$", "interpreter", "latex");
ylabel("Pitch angle - $\theta$ $[deg]$", "interpreter", "latex");
% -------------------------------------------------------------------------
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_a.theta_plot1.value = figure();
hold on 
grid on; grid minor; 
plot(t_t, theta(:,1), '-r.', 'linewidth', 1.5);
title("Pitch angle vs Time", "interpreter", "latex");
xlabel("Time - $t$ $[s]$", "interpreter", "latex");
ylabel("Pitch angle - $\theta$ $[deg]$", "interpreter", "latex");
xlim([0.0 10.0]);
ylim([0.0 10.0]);
% Aircraft.Geometry.Horizontal.S.value 

%% CS - VLA 423 - METHOD A - MANOEUVRING AIRSPEED VA 
% NOSE UP PITCHING 

% L TAIL 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.l_tail.value = Aircraft.Geometry.Horizontal.l.value; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.l_tail.Attributes.unit = "meters";

% MOMENT OF INERTIA 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.IY.value = Aircraft.Weight.I_Level.IY.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.IY.Attributes.unit = "kg * m^2"; 

% V = VA
% ANGULAR ACCELERATION PER CS VLA 423
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VA.pitch_up_acceleration.value = (20.1/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value - 1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VA.pitch_up_acceleration.Attributes.unit = "rad/sec^2";

% PITCHING MOMENT AT VA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VA.MA.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.IY.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VA.pitch_up_acceleration.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VA.MA.Attributes.unit = "N*m";

% TAIL LOAD
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VA.tail_load.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.MA.value/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.l_tail.value)*1e-1;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VA.tail_load.Attributes.unit = "daN";

% TOTAL LOAD 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VA.total_horizontal_load.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTail_A.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VAtail_load.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VA.total_horizontal_load.Attributes.unit = "daN";

% V = VC
% ANGULAR ACCELERATION PER CS VLA 423
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VC.pitch_up_acceleration.value = (20.1/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value - 1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VC.pitch_up_acceleration.Attributes.unit = "rad/sec^2";

% PITCHING MOMENT AT VA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VC.MC.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.IY.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VC.pitch_up_acceleration.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VC.MC.Attributes.unit = "N*m";

% TAIL LOAD
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VC.tail_load.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.MA.value/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.l_tail.value)*1e-1;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VC.tail_load.Attributes.unit = "daN";

% TOTAL LOAD 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VC.total_horizontal_load.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTail_C.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VC.tail_load.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_b.pitch_up.VC.total_horizontal_load.Attributes.unit = "daN";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VA - POINT A OF THE FLIGHT ENVELOPE - CASE 1

% Normal load n at Point A of the flight envelope
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.nA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.nA.Attributes.unit = "g's";

% Normal load n at Point A1, indicated in figure 1 of the CS - VLA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.nA1.value = 1.0; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.nA1.Attributes.unit = "g's";

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.DeltaN.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.nA.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.nA1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.DeltaN.Attributes.unit = "g's";

% Evaluation of the product Mg * Delta N
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.DeltaN.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Manoeuvring_weight.Attributes.unit = "N";

% Evaluation of the ratio X_cg/lt 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.X_cg_lt_ratio.Attributes.unit = "Non dimensional";

% Defining the downwash factor
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Downwash_factor.value = 1 - Aircraft.Geometry.Horizontal.DepsilonDalpha.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Downwash_factor.Attributes.unit = "Non dimensional";

% Defining the ratio S_ht/S - surface_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Surface_ratio.Attributes.unit = "Non dimensional";

% Defining the ratio a_ht/a - Lift_curve_slope_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Downwash_factor.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal_Tail_Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Third_term.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE ONE +++ 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.parenthesis.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VA - POINT A OF THE FLIGHT ENVELOPE - CASE 2

% Normal load n at Point A of the flight envelope
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.nA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.nA.Attributes.unit = "g's";

% Normal load n at Point A1, indicated in figure 1 of the CS - VLA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.nA1.value = 1.0; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.nA1.Attributes.unit = "g's";

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.DeltaN.value = -Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.nA.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.nA1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.DeltaN.Attributes.unit = "g's";

% Evaluation of the product Mg * Delta N
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.DeltaN.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Manoeuvring_weight.Attributes.unit = "N";

% Evaluation of the ratio X_cg/lt 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.X_cg_lt_ratio.Attributes.unit = "Non dimensional";

% Defining the downwash factor
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Downwash_factor.value = 1 - Aircraft.Geometry.Horizontal.DepsilonDalpha.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Downwash_factor.Attributes.unit = "Non dimensional";

% Defining the ratio S_ht/S - surface_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Surface_ratio.Attributes.unit = "Non dimensional";

% Defining the ratio a_ht/a - Lift_curve_slope_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Downwash_factor.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal_Tail_Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Third_term.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE TWO +++ 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.parenthesis.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VA - POINT A OF THE FLIGHT ENVELOPE - CASE 3

% Normal load n at Point A of the flight envelope
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.nG.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.nG.Attributes.unit = "g's";

% Normal load n at Point A1, indicated in figure 1 of the CS - VLA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.nA1.value = 1.0; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.nA1.Attributes.unit = "g's";

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.DeltaN.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.nG.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.nA1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.DeltaN.Attributes.unit = "g's";

% Evaluation of the product Mg * Delta N
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.DeltaN.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Manoeuvring_weight.Attributes.unit = "N";

% Evaluation of the ratio X_cg/lt 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.X_cg_lt_ratio.Attributes.unit = "Non dimensional";

% Defining the downwash factor
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Downwash_factor.value = 1 - Aircraft.Geometry.Horizontal.DepsilonDalpha.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Downwash_factor.Attributes.unit = "Non dimensional";

% Defining the ratio S_ht/S - surface_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Surface_ratio.Attributes.unit = "Non dimensional";

% Defining the ratio a_ht/a - Lift_curve_slope_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Downwash_factor.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal_Tail_Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Third_term.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE THREE +++ 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.parenthesis.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VA - POINT A OF THE FLIGHT ENVELOPE - CASE 4

% Normal load n at Point A of the flight envelope
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.nG.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.nG.Attributes.unit = "g's";

% Normal load n at Point A1, indicated in figure 1 of the CS - VLA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.nA1.value = 1.0; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.nA1.Attributes.unit = "g's";

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.DeltaN.value = - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.nG.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.nA1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.DeltaN.Attributes.unit = "g's";

% Evaluation of the product Mg * Delta N
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.DeltaN.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Manoeuvring_weight.Attributes.unit = "N";

% Evaluation of the ratio X_cg/lt 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.X_cg_lt_ratio.Attributes.unit = "Non dimensional";

% Defining the downwash factor
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Downwash_factor.value = 1 - Aircraft.Geometry.Horizontal.DepsilonDalpha.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Downwash_factor.Attributes.unit = "Non dimensional";

% Defining the ratio S_ht/S - surface_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Surface_ratio.Attributes.unit = "Non dimensional";

% Defining the ratio a_ht/a - Lift_curve_slope_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Downwash_factor.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal_Tail_Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Third_term.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE FOUR +++ 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.parenthesis.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VD - POINT D OF THE FLIGHT ENVELOPE - CASE 1

% Normal load n at Point A of the flight envelope
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.nA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.nA.Attributes.unit = "g's";

% Normal load n at Point A1, indicated in figure 1 of the CS - VLA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.nA1.value = 1.0; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.nA1.Attributes.unit = "g's";

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.DeltaN.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.nA.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.nA1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.DeltaN.Attributes.unit = "g's";

% Evaluation of the product Mg * Delta N
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.DeltaN.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Manoeuvring_weight.Attributes.unit = "N";

% Evaluation of the ratio X_cg/lt 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.X_cg_lt_ratio.Attributes.unit = "Non dimensional";

% Defining the downwash factor
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Downwash_factor.value = 1 - Aircraft.Geometry.Horizontal.DepsilonDalpha.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Downwash_factor.Attributes.unit = "Non dimensional";

% Defining the ratio S_ht/S - surface_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Surface_ratio.Attributes.unit = "Non dimensional";

% Defining the ratio a_ht/a - Lift_curve_slope_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Downwash_factor.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal_Tail_Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Third_term.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE ONE +++ 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.parenthesis.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VD - POINT A OF THE FLIGHT ENVELOPE - CASE 2

% Normal load n at Point A of the flight envelope
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.nA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.nA.Attributes.unit = "g's";

% Normal load n at Point A1, indicated in figure 1 of the CS - VLA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.nA1.value = 1.0; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.nA1.Attributes.unit = "g's";

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.DeltaN.value = -Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.nA.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.nA1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.DeltaN.Attributes.unit = "g's";

% Evaluation of the product Mg * Delta N
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.DeltaN.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Manoeuvring_weight.Attributes.unit = "N";

% Evaluation of the ratio X_cg/lt 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.X_cg_lt_ratio.Attributes.unit = "Non dimensional";

% Defining the downwash factor
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Downwash_factor.value = 1 - Aircraft.Geometry.Horizontal.DepsilonDalpha.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Downwash_factor.Attributes.unit = "Non dimensional";

% Defining the ratio S_ht/S - surface_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Surface_ratio.Attributes.unit = "Non dimensional";

% Defining the ratio a_ht/a - Lift_curve_slope_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Downwash_factor.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal_Tail_Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Third_term.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE TWO +++ 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.parenthesis.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VD - POINT E OF THE FLIGHT ENVELOPE - CASE 3

% Normal load n at Point E of the flight envelope
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.nE.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.nE.Attributes.unit = "g's";

% Normal load n at Point D1, indicated in figure 1 of the CS - VLA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.nD1.value = 1.0; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.nD1.Attributes.unit = "g's";

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.DeltaN.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.nE.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.nD1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.DeltaN.Attributes.unit = "g's";

% Evaluation of the product Mg * Delta N
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.DeltaN.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Manoeuvring_weight.Attributes.unit = "N";

% Evaluation of the ratio X_cg/lt 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.X_cg_lt_ratio.Attributes.unit = "Non dimensional";

% Defining the downwash factor
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Downwash_factor.value = 1 - Aircraft.Geometry.Horizontal.DepsilonDalpha.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Downwash_factor.Attributes.unit = "Non dimensional";

% Defining the ratio S_ht/S - surface_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Surface_ratio.Attributes.unit = "Non dimensional";

% Defining the ratio a_ht/a - Lift_curve_slope_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Downwash_factor.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal_Tail_Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Third_term.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE THREE +++ 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.parenthesis.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VD - POINT E OF THE FLIGHT ENVELOPE - CASE 4

% Normal load n at Point A of the flight envelope
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.nE.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.nE.Attributes.unit = "g's";

% Normal load n at Point A1, indicated in figure 1 of the CS - VLA
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.nD1.value = 1.0; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.nD1.Attributes.unit = "g's";

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.DeltaN.value = - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.nE.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.nD1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.DeltaN.Attributes.unit = "g's";

% Evaluation of the product Mg * Delta N
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.DeltaN.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Manoeuvring_weight.Attributes.unit = "N";

% Evaluation of the ratio X_cg/lt 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.X_cg_lt_ratio.Attributes.unit = "Non dimensional";

% Defining the downwash factor
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Downwash_factor.value = 1 - Aircraft.Geometry.Horizontal.DepsilonDalpha.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Downwash_factor.Attributes.unit = "Non dimensional";

% Defining the ratio S_ht/S - surface_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Surface_ratio.Attributes.unit = "Non dimensional";

% Defining the ratio a_ht/a - Lift_curve_slope_ratio
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Downwash_factor.value);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal_Tail_Normal_Force_Curve_Slope.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Third_term.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE FOUR +++ 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.parenthesis.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";

%% TOTAL AIRLOADS = L_tail + DELTA L_tail

% AIRSPEED VA - CASE 1 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTail_A.value;

% AIRSPEED VA - CASE 2 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTail_S.value;

% AIRSPEED VA - CASE 3 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTail_A.value;

% AIRSPEED VA - CASE 4 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTail_G.value;
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% AIRSPEED VD - CASE 1 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTail_F.value;

% AIRSPEED VD - CASE 2 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTail_E.value;

% AIRSPEED VD - CASE 3 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTail_E.value;

% AIRSPEED VD - CASE 4 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTail_D.value;

%% PRINT RESULTS 

% Horizontal tail loads increments
Increment = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Manoeuvring_Critical_Load_Increment.value, ...
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Manoeuvring_Critical_Load_Increment.value, ...
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Manoeuvring_Critical_Load_Increment.value, ...
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Manoeuvring_Critical_Load_Increment.value; ...
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Manoeuvring_Critical_Load_Increment.value, ...
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Manoeuvring_Critical_Load_Increment.value, ...
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Manoeuvring_Critical_Load_Increment.value, ...
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Manoeuvring_Critical_Load_Increment.value];
disp(" ++++++++++ Critical Horizontal Tail loads increments [daN] ++++++++++ ")
format = '%f          %f          %f          %f\n';
label  = 'Case 1              Case 2             Case 3             Case 4  \n';
fprintf(label);
fprintf(format, Increment.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

% Total horizontal tail increment
Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_one.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_two.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_three.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VA.case_four.Total_airloads.value; ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_one.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_two.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_three.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Option_d.VD.case_four.Total_airloads.value];
disp(" +++++++++++++++++ Total Horizontal Tail loads [daN] +++++++++++++++++ ")
format = '%f          %f          %f          %f\n';
label  = 'Case 1              Case 2             Case 3             Case 4  \n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

