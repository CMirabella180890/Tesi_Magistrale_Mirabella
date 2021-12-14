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

%% STRAIGHT FLIGHT
switch (Straight_flight_Case)
    % CASE 1: VA greater than the intercept
    case 'Case 1'
        % =================================================================
        if max(n_gust_cruise_plus) > nmax
        % =================================================================
        qA    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        VA    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
        LHTA  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.value;
        
        % LOAD FACTOR FROM C TO D
        n_fromCtoD = [n_fromCtoC2; n_fromC2toD];
        V_fromCtoD = [V_fromCtoC2; V_fromC2toD];
        
        % LOAD FACTOR FROM A TO C
        n_fromAtoC = [n_fromA1toC1; n_fromC1toC];
        V_fromAtoC = [V_fromA1toC1; V_fromC1toC];
        
        % FULL LOAD FACTOR VECTOR 
        full_load_factor_vector = [n_fromAtoC; n_fromCtoD];
        full_airspeed_vector    = [V_fromAtoC; V_fromCtoD];
        
        % UNIT LOAD FACTOR 
         V_unit_load_factor = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value;
         n_unit_load_factor = ones(length(V_unit_load_factor), 1);
         
        % =================================================================
        elseif max(n_gust_cruise_plus) < nmax
        % ================================================================= 
        qA   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        VA   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
        LHTA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.value;
        
        % LOAD FACTOR FROM C TO D
        % LOAD FACTOR FROM A TO C
        
        end
    % CASE 2: VA lower than the intercept
    case 'Case 2'
        qA   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        VA   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
        LHTA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.value;
        
        % LOAD FACTOR FROM C TO D
        n_fromCtoD = [n_fromCtoA2; n_fromA2toD];
        V_fromCtoD = [V_fromCtoA2; V_fromA2toD];
        
        % LOAD FACTOR FROM A TO C
        n_fromAtoC = n_fromA1toC; 
        V_fromAtoC = V_fromA1toC; 
        
        % FULL LOAD FACTOR VECTOR 
        full_load_factor_vector = [n_fromAtoC; n_fromCtoD];
        full_airspeed_vector    = [V_fromAtoC; V_fromCtoD];
        
        % UNIT LOAD FACTOR 
         V_unit_load_factor = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value;
         n_unit_load_factor = ones(length(V_unit_load_factor), 1);
end

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

% GENERAL VARIABLES FOR ALL THE CALCULATIONS
l_ht                  = Aircraft.Geometry.Horizontal.l.value;
S_ht                  = Aircraft.Geometry.Horizontal.S.value;
CLalfa_ht_rad         = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
CLalfa_ht_deg         = (CLalfa_ht_rad) * (pi/180);
CL_delta_elevator_rad = Aircraft.Certification.Aerodynamic_data.Horizontal.CL_delta_elevator.value;
CL_delta_elevator_deg = (CL_delta_elevator_rad) * (pi/180);
time_step             = 0.006;
time_interval         = 50;
time_vector           = linspace(0.0, time_interval * time_step, time_interval)';
damping_factor        = 0.3;
IY                    = Aircraft.Weight.I_Level.IY.value;
dvVA                  = 0.01*(180/pi);
MAC                   = Aircraft.Geometry.Wing.mac.value;
S                     = Aircraft.Geometry.Wing.S.value;

% RATIO L_TAIL / M.A.C. 
L_RATIO = l_ht / MAC;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.L_ratio.value = L_RATIO;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.L_ratio.Attributes.unit = "Non dimensional";

% RATIO S_TAIL / S_WING 
S_RATIO = S_ht / S;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_ratio.value = S_RATIO; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_ratio.Attributes.unit = "Non dimensional"; 

% HORIZONTAL TAIL VOLUME RATIO 
V_ht = L_RATIO * S_RATIO; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.Horizontal_Tail_Volume_Ratio.value = V_ht;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.Horizontal_Tail_Volume_Ratio.Attributes.unit = "Non dimensional";

% STORE INSIDE THE STRUCT VARIABLE 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_step.value = time_step;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_step.Attributes.unit = "seconds";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_interval_num.value = time_interval;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_interval_num.Attributes.unit = "Pure number";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value = time_vector;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.Attributes.unit = "seconds";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.damping_factor.value = damping_factor; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.damping_factor.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.IY.value = IY;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.IY.Attributes.unit = "kg * m^2"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dvVA.value = dvVA; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dvVA.Attributes.unit = "degree";

% AMC 23.423 ADVISED TOTAL DEFLECTION TIME INTERVAL 
% BE CAREFUL: AMC 23.423 Suggest the use of various total deflection time
% interval for different aircraft categories (Normal, Utility, Commuter,
% Acrobatic, ...). Pleas, take care of the definition of the time interval.
if (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Aerobatic") && (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag2 == "Stick")
    Aircraft.Geometry.Elevator.total_deflection_time.value = 0.1;
    Aircraft.Geometry.Elevator.total_deflection_time.Attributes = "seconds";
elseif (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Aerobatic") && (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag2 == "Wheel")
    Aircraft.Geometry.Elevator.total_deflection_time.value = 0.2;
    Aircraft.Geometry.Elevator.total_deflection_time.Attributes = "seconds";
elseif (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Normal") | (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Utility") | (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Normal")
    if Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag2 == "Stick"
        Aircraft.Geometry.Elevator.total_deflection_time.value = 0.2;
        Aircraft.Geometry.Elevator.total_deflection_time.Attributes = "seconds";
    elseif Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag2 == "Wheel"
        Aircraft.Geometry.Elevator.total_deflection_time.value = 0.3;
        Aircraft.Geometry.Elevator.total_deflection_time.Attributes = "seconds";
    end
elseif Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Normal"
    prmpt = "Enter total defl. time interval --> t_total_defl_time: ";
    Aircraft.Geometry.Elevator.total_deflection_time.value = input(prmpt);
    Aircraft.Geometry.Elevator.total_deflection_time.Attributes = "seconds";
end

total_deflection_time = Aircraft.Geometry.Elevator.total_deflection_time.value;

% SOLVING THE EQUATION OF MOTION 
% 
%     d^2 theta   (q * S_tail * a_tail * d)   /           delta_v     \
%     --------- = ------------------------- *| omega*dt - ------- * DF|
%       dt^2                 IY               \             VA        /
%
delta_elevator_max     = Aircraft.Geometry.Elevator.max_deflection.value;
tau                    = Aircraft.Certification.Aerodynamic_data.Horizontal.tau.value;
omega_deg              = (delta_elevator_max * tau)*(1/total_deflection_time);
omega_rad              = deg2rad(omega_deg);
A0                     = (1/IY) * CLalfa_ht_rad * qA * l_ht * S_ht;
alpha_prime_horiz_rad  = omega_rad * time_vector;
alpha_prime_horiz_deg  = rad2deg(alpha_prime_horiz_rad);

% STORE INSIDE THE STRUCT VARIABLE
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.omega.value = omega_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.omega.Attributes.unit = "deg/sec";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.omega_rad.value = omega_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.omega_rad.Attributes.unit = "rad/sec";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.A0.value = A0;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.A0.Attributes.unit = "1/(rad*m*s^2)";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.alpha_prime_horiz.value = alpha_prime_horiz_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.alpha_prime_horiz.Attributes.unit = "radians";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.alpha_prime_horiz_deg.value = alpha_prime_horiz_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.alpha_prime_horiz_deg.Attributes.unit = "degrees";

% SOLVING THE DIFFERENTIAL EQUATION - PITCH DOWN CASE
d2thetadt2          = zeros(length(time_vector), 1);
dthetadt            = zeros(length(time_vector), 1);
alpha_new_horiz_rad = zeros(length(time_vector), 1);
delta_theta_rad     = zeros(length(time_vector), 1);
delta_v             = zeros(length(time_vector), 1);

% SOLUTION SCHEME 
 for i = 2:length(time_vector)
    d2thetadt2(i)          = A0 * (alpha_prime_horiz_rad(i) - delta_theta_rad(i-1));
    dthetadt(i)            = dthetadt(i-1) + 0.5 * (d2thetadt2(i-1) + d2thetadt2(i))*(time_step);
    delta_v(i)             = dthetadt(i) * l_ht;
    delta_theta_rad(i)     = delta_v(i) * (damping_factor / VA);
    alpha_new_horiz_rad(i) = alpha_prime_horiz_rad(i)  - delta_theta_rad(i);
 end

% STORE INSIDE THE STRUCT VARIABLE
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.d2thetadt2.value = d2thetadt2;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.d2thetadt2.Attributes.unit = "rad/sec^2"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.dthetadt.value = dthetadt;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.dthetadt.Attributes.unit = "rad/sec"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.alpha_new_horiz.value = alpha_new_horiz_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.alpha_new_horiz.Attributes.unit = "rad"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.delta_theta.value = delta_theta_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.delta_theta.Attributes.unit = "rad"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.delta_v.value = delta_v;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.delta_v.Attributes.unit = "m/s";  
 
% ALPHA NEW 
alpha_new_horiz_deg = rad2deg(alpha_new_horiz_rad);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.alpha_new_horiz_grad.value = alpha_new_horiz_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.alpha_new_horiz_grad.Attributes.unit = "deg";

% LIMIT HORIZONTAL TAIL LOAD 
DeltaLimitLTail_pitch_down = alpha_new_horiz_deg(end) * CLalfa_ht_deg * S_ht * qA *(1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.DeltaLimitLTail.value = DeltaLimitLTail_pitch_down;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.DeltaLimitLTail.Attributes.unit = "daN";

% SOLVING THE DIFFERENTIAL EQUATION - PITCH UP CASE
% A WORD OF CAUTION: The maximum negative deflection angle of the elevator
% is delta_e_upward = - 20.0 deg = -0.8 * delta_max_down 
omega_deg             = (-0.8 * delta_elevator_max * tau )*(1/ total_deflection_time );
omega_rad             = deg2rad(omega_deg);
A0                    = (1/ IY )* CLalfa_ht_rad * qA * l_ht * S_ht;
alpha_prime_horiz_rad = omega_rad * time_vector;
alpha_prime_horiz_deg = rad2deg(alpha_prime_horiz_rad);

% STORE INSIDE THE STRUCT VARIABLE
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.omega.value = omega_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.omega.Attributes.unit = "deg/sec";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.omega_rad.value = omega_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.omega_rad.Attributes.unit = "rad/sec";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.A0.value = A0;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.A0.Attributes.unit = "1/(rad*m*s^2)";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.alpha_prime_horiz.value = alpha_prime_horiz_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.alpha_prime_horiz.Attributes.unit = "radians";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.alpha_prime_horiz_deg.value =alpha_prime_horiz_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.alpha_prime_horiz_deg.Attributes.unit = "degrees";

% SOLVING THE DIFFERENTIAL EQUATION - PITCH UP CASE
d2thetadt2          = zeros(length(time_vector), 1);
dthetadt            = zeros(length(time_vector), 1);
alpha_new_horiz_rad = zeros(length(time_vector), 1);
delta_theta_rad     = zeros(length(time_vector), 1);
delta_v             = zeros(length(time_vector), 1);

% SOLUTION SCHEME 
 for i = 2:length(time_vector)
    d2thetadt2(i)          = A0 * (alpha_prime_horiz_rad(i) - delta_theta_rad(i-1));
    dthetadt(i)            = dthetadt(i-1) + 0.5 * (d2thetadt2(i-1) + d2thetadt2(i))*(time_step);
    delta_v(i)             = dthetadt(i) * l_ht;
    delta_theta_rad(i)     = delta_v(i) * (damping_factor / VA);
    alpha_new_horiz_rad(i) = alpha_prime_horiz_rad(i)  - delta_theta_rad(i);
 end

% SOLVING THE DIFFERENTIAL EQUATION 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.d2thetadt2.value = d2thetadt2;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.d2thetadt2.Attributes.unit = "rad/sec^2"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.dthetadt.value = dthetadt;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.dthetadt.Attributes.unit = "rad/sec"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.alpha_new_horiz.value = alpha_new_horiz_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.alpha_new_horiz.Attributes.unit = "rad"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.delta_theta.value = delta_theta_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.delta_theta.Attributes.unit = "rad"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.delta_v.value = delta_v;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.delta_v.Attributes.unit = "m/s"; 

% ALPHA NEW DEG
alpha_new_horiz_deg = rad2deg(alpha_new_horiz_rad);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.alpha_new_horiz_grad.value = alpha_new_horiz_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.alpha_new_horiz_grad.Attributes.unit = "deg";

% LIMIT HORIZONTAL TAIL LOAD 
DeltaLimitLTail_pitch_up = alpha_new_horiz_deg(end) * CLalfa_ht_deg * S_ht * qA * (1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.DeltaLimitLTail.value = DeltaLimitLTail_pitch_up;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.DeltaLimitLTail.Attributes.unit = "daN";

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
disp( " ")
disp(" ++++ METHOD CS - VLA 423 (a) ++++");
disp(" ---------------------------------");
disp(" ++++ UNCHECKED ++++");

% Horizontal tail loads increments
Increment = [DeltaLimitLTail_pitch_up, ...
             DeltaLimitLTail_pitch_down];
disp(" ++++++++++ Critical Horizontal Tail loads increments [daN] ++++++++++ ")
format = '%6.6f          %6.6f\n';
label  = ' Pitch-up at VA     Pitch-down at VA\n';
fprintf(label);
fprintf(format, Increment.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

% CRITICAL LOAD FOR METHOD CS - VLA 423 (a) 
if abs(DeltaLimitLTail_pitch_up) > abs(DeltaLimitLTail_pitch_down)
    Total_critical_load = LHTA + DeltaLimitLTail_pitch_up;
    
elseif abs(DeltaLimitLTail_pitch_down) > abs(DeltaLimitLTail_pitch_up)
    Total_critical_load = LHTA + DeltaLimitLTail_pitch_down;
end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.Total_critical_load.value = Total_critical_load;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.Total_critical_load.Attributes.unit = "daN";

disp(" ")
% Horizontal tail loads increments
Increment = Total_critical_load;
disp(" ++++++++++ Total Critical Horizontal Tail loads - Method (a) [daN] ++++++++++ ")
format = ' %6.6f\n';
label  = ' Total critical load at VA\n';
fprintf(label);
fprintf(format, Increment.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

%% CS - VLA 423 - METHOD B - PITCH UP

% V UNIT LOAD FACTOR
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value = V_unit_load_factor;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.Attributes.unit = "m/s";

% FULL LOAD FACTOR VECTOR 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.full_load_factor_vector.value = full_load_factor_vector;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.full_load_factor_vector.Attributes.unit = "g's";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.full_airspeed_vector.value = full_airspeed_vector;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.full_airspeed_vector.Attributes.unit = "m/s";

% ANGULAR ACCELERATION CALCULATIONS 
angular_acceleration_pitch_up = zeros(length(full_load_factor_vector), 1);
for i = 1:length(full_airspeed_vector)
    lm = ( 20.1 / full_airspeed_vector(i) ) * full_load_factor_vector(i) * ( full_load_factor_vector(i) - 1.5 );
    angular_acceleration_pitch_up(i) = lm;
end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.ang_acc.value = angular_acceleration_pitch_up;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.ang_acc.Attributes.unit = "rad/sec^2";

% MOMENT MY
MY = zeros(length(angular_acceleration_pitch_up), 1);

for i = 1:length(angular_acceleration_pitch_up)
    MY(i) = ( IY ) * ( angular_acceleration_pitch_up(i) );
end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.MY.value = MY;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.MY.Attributes.unit = "N*m";

% DELTA TAIL AIRLOADS
delta_tail_airloads = (1e-1)*( - MY ) / ( l_ht);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.delta_tail_airloads.value = delta_tail_airloads; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.delta_tail_airloads.Attributes.unit = "daN";

% BALANCING TAIL AIRLOADS

LHT_unit_load_factor    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_unit_load_factor.value;
index_va                = dsearchn(V_unit_load_factor, VA);
balancing_tail_airloads = linspace(LHT_unit_load_factor(index_va), LHT_unit_load_factor(end), length(delta_tail_airloads))';
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.balancing_tail_airloads.value = balancing_tail_airloads;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.balancing_tail_airloads.Attributes.unit = "daN";

% TOTAL AIRLOADS ASSOCIATED WITH THE PITCH UP MANOEUVRE
total_tail_airloads = delta_tail_airloads + balancing_tail_airloads;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.total_tail_airloads.value = total_tail_airloads;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.total_tail_airloads.Attributes.unit = "daN";

% CRITICAL TAIL AIRLOADS VALUE 
critical_tail_airloads = -max(abs(total_tail_airloads)); 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.value = critical_tail_airloads;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.Attributes.unit = "daN";

disp(" ");
disp(" ++++ METHOD CS - VLA 423 (b) ++++");
disp(" ---------------------------------");
disp(" ++++ CHECKED ++++");

% Total horizontal tail increment
Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.value];
disp(" +++++++++++++++++ Total Horizontal Tail loads [daN] +++++++++++++++++ ")
format = ' %6.6f\n';
label  = ' At VA\n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

%% ========================================================================================================================        
% 
% % METHOD (a) OF THE CS VLA 423
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.l_tail.value = Aircraft.Geometry.Horizontal.l.value; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.l_tail.Attributes.unit = "meters";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_tail.value = Aircraft.Geometry.Horizontal.S.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_tail.Attributes.unit = "squared meters";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.qA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.qA.Attributes.unit = "Pa";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.a_tail_rad.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.a_tail_rad.Attributes.unit = "1/radians";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.a_tail_grad.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value*(pi/180);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.a_tail_grad.Attributes.unit = "1/degrees";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.CL_delta_elevator.value = Aircraft.Certification.Aerodynamic_data.Horizontal.CL_delta_elevator.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.CL_delta_elevator.Attributes.unit = "1/radians"; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.CL_delta_elevator_grad.value = Aircraft.Certification.Aerodynamic_data.Horizontal.CL_delta_elevator.value*(pi/180);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.CL_delta_elevator_grad.Attributes.unit = "1/degrees";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_step.value = 0.006;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_step.Attributes.unit = "seconds";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_interval_num.value = 50;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_interval_num.Attributes.unit = "Pure number";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value = linspace(0, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_interval_num.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_step.value, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_interval_num.value)';
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.Attributes.unit = "seconds";
% 
% % AMC 23.423 ADVISED TOTAL DEFLECTION TIME INTERVAL 
% % BE CAREFUL: AMC 23.423 Suggest the use of various total deflection time
% % interval for different aircraft categories (Normal, Utility, Commuter,
% % Acrobatic, ...). Pleas, take care of the definition of the time interval.
% if (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Aerobatic") && (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag2 == "Stick")
%     Aircraft.Geometry.Elevator.total_deflection_time.value = 0.1;
%     Aircraft.Geometry.Elevator.total_deflection_time.Attributes = "seconds";
% elseif (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Aerobatic") && (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag2 == "Wheel")
%     Aircraft.Geometry.Elevator.total_deflection_time.value = 0.2;
%     Aircraft.Geometry.Elevator.total_deflection_time.Attributes = "seconds";
% elseif (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Normal") | (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Utility") | (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Normal")
%     if Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag2 == "Stick"
%         Aircraft.Geometry.Elevator.total_deflection_time.value = 0.2;
%         Aircraft.Geometry.Elevator.total_deflection_time.Attributes = "seconds";
%     elseif Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag2 == "Wheel"
%         Aircraft.Geometry.Elevator.total_deflection_time.value = 0.3;
%         Aircraft.Geometry.Elevator.total_deflection_time.Attributes = "seconds";
%     end
% elseif Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Normal"
%     prmpt = "Enter total defl. time interval --> t_total_defl_time: ";
%     Aircraft.Geometry.Elevator.total_deflection_time.value = input(prmpt);
%     Aircraft.Geometry.Elevator.total_deflection_time.Attributes = "seconds";
% end
% 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.total_time_of_deflection.value = 0.3;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.total_time_of_deflection.Attributes.unit = "sec";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.damping_factor.value = 0.3; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.damping_factor.Attributes.unit = "Non dimensional";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.IY.value = Aircraft.Weight.I_Level.IY.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.IY.Attributes.unit = "kg * m^2"; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dvVA.value = 0.01*(180/pi); 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dvVA.Attributes.unit = "degree";
% 
% % RATIO L_TAIL / M.A.C. 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.L_ratio.value =  Aircraft.Geometry.Horizontal.l.value/Aircraft.Geometry.Wing.mac.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.L_ratio.Attributes.unit = "Non dimensional";
% 
% % RATIO S_TAIL / S_WING 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_ratio.Attributes.unit = "Non dimensional"; 
% 
% % HORIZONTAL TAIL VOLUME RATIO 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.Horizontal_Tail_Volume_Ratio.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.L_ratio.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_ratio.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.Horizontal_Tail_Volume_Ratio.Attributes.unit = "Non dimensional";
% 
% % SOLVING THE EQUATION OF MOTION 
% % 
% %     d^2 theta   (q * S_tail * a_tail * d)   /           delta_v     \
% %     --------- = ------------------------- *| omega*dt - ------- * DF|
% %       dt^2                 IY               \             VA        /
% %
% 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega.value = (Aircraft.Geometry.Elevator.max_deflection.value*Aircraft.Certification.Aerodynamic_data.Horizontal.tau.value)*(1/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.total_time_of_deflection.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega.Attributes.unit = "deg/sec";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega_rad.Attributes.unit = "rad/sec";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.A0.value = (1/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.IY.value)*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.a_tail_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.qA.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.l_tail.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_tail.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.A0.Attributes.unit = "1/(rad*m*s^2)";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz.Attributes.unit = "radians";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz_deg.Attributes.unit = "degrees";
% 
% % SOLVING THE DIFFERENTIAL EQUATION 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.d2thetadt2.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.d2thetadt2.Attributes.unit = "rad/sec^2"; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dthetadt.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dthetadt.Attributes.unit = "rad/sec"; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz.Attributes.unit = "rad"; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_theta.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_theta.Attributes.unit = "rad"; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_v.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_v.Attributes.unit = "m/s"; 
% 
%  for i = 2:length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value)
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.d2thetadt2.value(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.A0.value*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz.value(i) - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_theta.value(i-1));
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dthetadt.value(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dthetadt.value(i-1) + 0.5*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.d2thetadt2.value(i-1)+Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.d2thetadt2.value(i))*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_step.value);
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_v.value(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dthetadt.value(i)*Aircraft.Geometry.Horizontal.l.value;
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_theta.value(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_v.value(i)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.damping_factor.value/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value);
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz.value(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz.value(i)  - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_theta.value(i);
%  end
% 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz_grad.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz_grad.Attributes.unit = "deg";
% 
% % LIMIT HORIZONTAL TAIL LOAD 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.DeltaLimitLTail.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz_grad.value(end)*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.a_tail_grad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_tail.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.qA.value*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.DeltaLimitLTail.Attributes.unit = "daN";
% 
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % PITCH UP CASE 
% % A WORD OF CAUTION: The maximum negative deflection angle of the elevator
% % is delta_e_upward = - 20.0 deg = -0.8 * delta_max_down 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega.value = (-0.8*Aircraft.Geometry.Elevator.max_deflection.value*Aircraft.Certification.Aerodynamic_data.Horizontal.tau.value)*(1/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.total_time_of_deflection.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega.Attributes.unit = "deg/sec";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega_rad.Attributes.unit = "rad/sec";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.A0.value = (1/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.IY.value)*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.a_tail_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.qA.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.l_tail.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_tail.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.A0.Attributes.unit = "1/(rad*m*s^2)";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz.Attributes.unit = "radians";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.omega_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz_deg.Attributes.unit = "degrees";
% 
% % SOLVING THE DIFFERENTIAL EQUATION 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.d2thetadt2.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.d2thetadt2.Attributes.unit = "rad/sec^2"; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dthetadt.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dthetadt.Attributes.unit = "rad/sec"; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz.Attributes.unit = "rad"; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_theta.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_theta.Attributes.unit = "rad"; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_v.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_v.Attributes.unit = "m/s"; 
% 
%  for i = 2:length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value)
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.d2thetadt2.value(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.A0.value*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz.value(i) - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_theta.value(i-1));
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dthetadt.value(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dthetadt.value(i-1) + 0.5*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.d2thetadt2.value(i-1)+Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.d2thetadt2.value(i))*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_step.value);
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_v.value(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dthetadt.value(i)*Aircraft.Geometry.Horizontal.l.value;
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_theta.value(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_v.value(i)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.damping_factor.value/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value);
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz.value(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_prime_horiz.value(i)  - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.delta_theta.value(i);
%  end
% 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz_grad.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz_grad.Attributes.unit = "deg";
% 
% % LIMIT HORIZONTAL TAIL LOAD 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.DeltaLimitLTail.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.alpha_new_horiz_grad.value(end)*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.a_tail_grad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.S_tail.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.qA.value*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.DeltaLimitLTail.Attributes.unit = "daN";
% 
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
% disp( " ")
% disp(" ++++ METHOD CS - VLA 423 (a) ++++");
% disp(" ---------------------------------");
% disp(" ++++ UNCHECKED ++++");
% 
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.DeltaLimitLTail.value, ...
%              Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.DeltaLimitLTail.value];
% disp(" ++++++++++ Critical Horizontal Tail loads increments [daN] ++++++++++ ")
% format = '%6.6f          %6.6f\n';
% label  = ' Pitch-up at VA     Pitch-down at VA\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % CRITICAL LOAD FOR METHOD CS - VLA 423 (a) 
% if abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.DeltaLimitLTail.value) > abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.DeltaLimitLTail.value)
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.Total_critical_load.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.DeltaLimitLTail.value;
% elseif abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.DeltaLimitLTail.value) > abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.DeltaLimitLTail.value)
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.Total_critical_load.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.DeltaLimitLTail.value;
% end
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.Total_critical_load.Attributes.unit = "daN";
% 
% disp(" ")
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.Total_critical_load.value];
% disp(" ++++++++++ Total Critical Horizontal Tail loads - Method (a) [daN] ++++++++++ ")
% format = ' %6.6f\n';
% label  = ' Total critical load at VA\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% %% CS - VLA 423 - METHOD B - PITCH UP
% % LOAD FACTOR FROM C TO D
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.n_fromCtoD.value = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.value, ...
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nD.value, ...
%     length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value))';
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.n_fromCtoD.Attributes.unit = "g's";
% 
% % LOAD FACTOR FROM A TO C
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.n_fromAtoC.value = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value, ...
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nC.value, ...
% length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value))';
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.n_fromAtoC.Attributes.unit = "g's";
% 
% % FULL LOAD FACTOR 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.full_load_factor_vector.value = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.n_fromAtoC.value; Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.n_fromCtoD.value];
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.full_load_factor_vector.Attributes.unit = "g's";
% 
% % UNIT LOAD FACTOR 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.unit_load_factor_vector.value = ones(2*length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value), 1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.unit_load_factor_vector.Attributes.unit = "g's";
% 
% % AIRSPEED VECTOR FOR CALCULATIONS
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.airspeed_vector.value = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value, ...
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value, ...
%     2*length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value))';
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.airspeed_vector.Attributes.unit = "m/s";
% 
% % ANGULAR ACCELERATION CALCULATIONS 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.ang_acc.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.airspeed_vector.value), 1);
% v  = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.airspeed_vector.value;
% nm = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.full_load_factor_vector.value;
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.airspeed_vector.value)
%     lm = (20.1/v(i))*nm(i)*(nm(i) - 1.5);
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.ang_acc.value(i) = lm;
% end
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.ang_acc.Attributes.unit = "rad/sec^2";
% clear v nm lm
% 
% % MOMENT MY
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.MY.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.ang_acc.value), 1);
% 
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.ang_acc.value)
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.MY.value(i) = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.IY.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.ang_acc.value(i));
% end
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.MY.Attributes.unit = "N*m";
% 
% % DELTA TAIL AIRLOADS
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.delta_tail_airloads.value = (1e-1)*(-Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.MY.value)/(Aircraft.Geometry.Horizontal.l.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.delta_tail_airloads.Attributes.unit = "daN";
% 
% % BALANCING TAIL AIRLOADS
% index_va = dsearchn(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.balancing_tail_airloads.value = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_unit_load_factor.value(index_va), ...
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_unit_load_factor.value(end), ...
%     length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.delta_tail_airloads.value))';
% % Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.balancing_tail_airloads.value = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_positivestall.value; Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromCtoD.value];
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.balancing_tail_airloads.Attributes.unit = "daN";
% 
% % TOTAL AIRLOADS ASSOCIATED WITH THE PITCH UP MANOEUVRE
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.total_tail_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.delta_tail_airloads.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.balancing_tail_airloads.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.total_tail_airloads.Attributes.unit = "daN";
% 
% % CRITICAL TAIL AIRLOADS VALUE 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.value = -max(abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.total_tail_airloads.value)); 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.Attributes.unit = "daN";
% 
% disp(" ");
% disp(" ++++ METHOD CS - VLA 423 (b) ++++");
% disp(" ---------------------------------");
% disp(" ++++ CHECKED ++++");
% 
% % Total horizontal tail increment
% Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.value];
% disp(" +++++++++++++++++ Total Horizontal Tail loads [daN] +++++++++++++++++ ")
% format = ' %6.6f\n';
% label  = ' At VA\n';
% fprintf(label);
% fprintf(format, Total.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% %% CS - VLA 423 - METHOD B - PITCH DOWN
% 
% % ANGULAR ACCELERATION CALCULATIONS 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.ang_acc.value = -Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.ang_acc.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.ang_acc.Attributes.unit = "rad/sec^2";
% 
% % MOMENT MY
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.MY.value = -Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.MY.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.MY.Attributes.unit = "N*m";
% 
% % DELTA TAIL AIRLOADS
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.delta_tail_airloads.value = (1e-1)*(-Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.MY.value)/(Aircraft.Geometry.Horizontal.l.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.delta_tail_airloads.Attributes.unit = "daN";
% 
% % BALANCING TAIL AIRLOADS
% index_va = dsearchn(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.balancing_tail_airloads.value = linspace(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.value(index_va), ...
%     Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LW_from0toS_new.value(end), ...
%     length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.delta_tail_airloads.value))';
% % Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.balancing_tail_airloads.value = [Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_positivestall.value; Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromCtoD.value];
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.balancing_tail_airloads.Attributes.unit = "daN";
% 
% % TOTAL AIRLOADS ASSOCIATED WITH THE PITCH UP MANOEUVRE
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.total_tail_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.delta_tail_airloads.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.balancing_tail_airloads.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.total_tail_airloads.Attributes.unit = "daN";
% 
% % CRITICAL TAIL AIRLOADS VALUE 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.critical_tail_airloads.value = max(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.total_tail_airloads.value); 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.critical_tail_airloads.Attributes.unit = "daN";
% 
% % Total horizontal tail increment
% disp(" ");
% disp(" ++++ METHOD CS - VLA 423 (b) - PITCH DOWN ++++")
% disp(" ----------------------------------------------")
% Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.critical_tail_airloads.value];
% disp(" +++++++++++++++++ Horizontal Tail loads [daN] +++++++++++++++++ ")
% format = ' %6.6f\n';
% label  = ' At VA\n';
% fprintf(label);
% fprintf(format, Total.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% if abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.value) > abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.critical_tail_airloads.value)
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.Total_critical_loads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.value;
% elseif abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.critical_tail_airloads.value) > abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.value)
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.Total_critical_loads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.critical_tail_airloads.value;
% end
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.Total_critical_loads.Attributes.unit = "daN";
% 
% % Total horizontal tail increment
% disp(" ");
% disp(" ++++ METHOD CS - VLA 423 (b) - PITCH DOWN ++++")
% disp(" ----------------------------------------------")
% Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.Total_critical_loads.value];
% disp(" +++++++++++++++++ Total Horizontal Tail loads - Method (b) [daN] +++++++++++++++++ ")
% format = ' %6.6f\n';
% label  = ' At VA\n';
% fprintf(label);
% fprintf(format, Total.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
%     %% CS - VLA 423 - METHOD (A)+(B) - SWITCH-CASE TO ASSESS THE LOWEST (IN MODULE) LOAD ON THE HORIZ. TAIL 
% tl_0 = abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.DeltaLimitLTail.value);
% tl_1 = abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.DeltaLimitLTail.value);
% tl_2 = abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.value);
% tl_3 = abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.critical_tail_airloads.value);
% if (tl_0 < tl_1) && (tl_0 < tl_2) && (tl_0 < tl_3)
%     % disp(" ++++ CRITICAL CONDITION: METHOD (a) PITCH UP ++++")
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.Attributes.flag = "METHOD CS-VLA 423 (A) PITCH UP";
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.DeltaLimitLTail.value;    
% elseif (tl_1 < tl_2) && (tl_1 < tl_3) && (tl_1 < tl_0)
%     % disp(" ++++ CRITICAL CONDITION: METHOD (a) PITCH DOWN ++++")
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.Attributes.flag = "METHOD CS-VLA 423 (A) PITCH DOWN";
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.DeltaLimitLTail.value;
% elseif (tl_2 < tl_1) && (tl_2 < tl_3) && (tl_2 < tl_0)
%     % disp(" ++++ CRITICAL CONDITION: METHOD (b) PITCH UP ++++")
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.Attributes.flag = "METHOD CS-VLA 423 (B) PITCH UP";
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.value;
% elseif (tl_3 < tl_1) && (tl_3 < tl_2) && (tl_3 < tl_0)
%     % disp(" ++++ CRITICAL CONDITION: METHOD (b) PITCH DOWN ++++")
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.Attributes.flag = "METHOD CS-VLA 423 (B) PITCH DOWN";
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.critical_tail_airloads.value;
% end
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.Attributes.unit = "daN";
% 
% 
% % CRITICAL TAIL AIRLOADS
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.Total_critical_tail_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.Total_critical_tail_airloads.Attributes.unit = "daN";
% 
% % Total horizontal tail increment
% disp(" ");
% disp(" ++++ METHOD CS - VLA 423 (a)+(b) ++++")
% disp(" ----------------------------------------------")
% Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.Total_critical_tail_airloads.value];
% disp(" +++++++++++++++++ Total Horizontal Tail loads [daN] +++++++++++++++++ ")
% format = ' %6.6f\n';
% label  = ' At VA\n';
% fprintf(label);
% fprintf(format, Total.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% %% CS - VLA 423 - METHOD (C) 
% % ############### VA - MAX UPWARD DEFLECTION - NEGATIVE DELTA ELEVATOR ###############
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Airspeed_at_PointA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Airspeed_at_PointA.Attributes.unit = "m/sec";
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_unit_load_factor.value)   
%     V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value(i);
%     if abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Airspeed_at_PointA.value - V) < 1e-1
%         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.BalancingLoad.value = (0.5)*(V^2)* ... 
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_unit_load_factor.value(i))*(1E-1);   
%     end
% end
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.BalancingLoad.Attributes.unit = "daN";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.S_tail.value = Aircraft.Geometry.Horizontal.S.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.S_tail.Attributes.unit = "squared meters";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.qA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.qA.Attributes.unit = "Pa";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.a_tail_rad.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.a_tail_rad.Attributes.unit = "1/radians";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.a_tail_grad.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value*(pi/180);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.a_tail_grad.Attributes.unit = "1/degrees";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.total_time_of_deflection.value = 0.3;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.total_time_of_deflection.Attributes.unit = "sec";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.damping_factor.value = 0.3; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.damping_factor.Attributes.unit = "Non dimensional";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.IY.value = Aircraft.Weight.I_Level.IY.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.IY.Attributes.unit = "kg * m^2"; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.dvVA.value = 0.01*(180/pi); 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.dvVA.Attributes.unit = "degree";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_interval_num.value = 50;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_interval_num.Attributes.unit = "Pure number";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_step.value = 0.006;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_step.Attributes.unit = "seconds";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_vector.value = linspace(0, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_interval_num.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_step.value, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_interval_num.value)';
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_vector.Attributes.unit = "seconds";
% 
% % RATIO L_TAIL / M.A.C. 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.L_ratio.value =  Aircraft.Geometry.Horizontal.l.value/Aircraft.Geometry.Wing.mac.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.L_ratio.Attributes.unit = "Non dimensional";
% 
% % RATIO S_TAIL / S_WING 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.S_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.S_ratio.Attributes.unit = "Non dimensional"; 
% 
% % HORIZONTAL TAIL VOLUME RATIO 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Horizontal_Tail_Volume_Ratio.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.L_ratio.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.S_ratio.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Horizontal_Tail_Volume_Ratio.Attributes.unit = "Non dimensional";
% 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega.value = ((-0.8*Aircraft.Geometry.Elevator.max_deflection.value)*Aircraft.Certification.Aerodynamic_data.Horizontal.tau.value)*(1/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.total_time_of_deflection.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega.Attributes.unit = "deg/sec";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega_rad.Attributes.unit = "rad/sec";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.A0.value = (1/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.IY.value)*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.a_tail_rad.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.A0.Attributes.unit = "1/(rad*m*s^2)";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_vector.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz.Attributes.unit = "radians";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_vector.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz_deg.Attributes.unit = "degrees";
% 
% % SOLVING THE DIFFERENTIAL EQUATION 
% d2thetadt2 = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% dthetadt = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1); 
% alpha_new_horiz = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% delta_theta = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% delta_v = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% 
% % OUTPUT 
%  for i = 2:length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_vector.value)
%     d2thetadt2(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.A0.value*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz.value(i) - delta_theta(i-1));
%     dthetadt(i) = dthetadt(i-1) + 0.5*(d2thetadt2(i-1)+d2thetadt2(i))*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_step.value);
%     delta_v(i) = dthetadt(i)*Aircraft.Geometry.Horizontal.l.value;
%     delta_theta(i) = delta_v(i)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.damping_factor.value/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value);
%     alpha_new_horiz(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz.value(i)  - delta_theta(i);
%  end
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.Results.value = [d2thetadt2, dthetadt, delta_v, delta_theta, alpha_new_horiz];
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.Results.Attributes.contents = "d2thetadt2, dthetadt, delta_v, delta_theta, alpha_new_horiz";
% 
% % CONVERSIONE TO DEGRESS 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.alpha_new_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.Results.value(:,5));
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.alpha_new_deg.Attributes.unit = "deg";
% 
% % CALCULATION OF THE CRITICAL HORIZONTAL TAIL AIRLOADS AS 
% % L_BALANC AT VA + CRITICAL TAIL AIRLOADS DUE TO MANOEUVRE AT VA 
% % FIRSE WE NEED THE CRITICAL TAIL AIRLOADS 
% % LIMIT HORIZONTAL TAIL LOAD 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.DeltaLHorizoTail.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.alpha_new_deg.value(end)*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.a_tail_grad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.S_tail.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.qA.value*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.DeltaLHorizoTail.Attributes.unit = "daN";
% 
% % TOTAL AIRLOADS ACTING ON THE HORIZONTAL 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.TotalLoads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.BalancingLoad.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.DeltaLHorizoTail.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.TotalLoads.Attributes.unit = "daN";
% 
% % ############### VA - MAX DOWNWARD DEFLECTION - POSITIVE DELTA ELEVATOR ###############
% % POSITIVE MAXIMUM DEFLECTION ANGLE
% 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega.value = ((Aircraft.Geometry.Elevator.max_deflection.value)*Aircraft.Certification.Aerodynamic_data.Horizontal.tau.value)*(1/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.total_time_of_deflection.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega.Attributes.unit = "deg/sec";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega_rad.Attributes.unit = "rad/sec";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.A0.value = (1/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.IY.value)*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.a_tail_rad.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.A0.Attributes.unit = "1/(rad*m*s^2)";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_vector.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz.Attributes.unit = "radians";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_vector.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz_deg.Attributes.unit = "degrees";
% 
% % SOLVING THE DIFFERENTIAL EQUATION 
% d2thetadt2 = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% dthetadt = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1); 
% alpha_new_horiz = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% delta_theta = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% delta_v = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% 
% % OUTPUT 
%  for i = 2:length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_vector.value)
%     d2thetadt2(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.A0.value*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz.value(i) - delta_theta(i-1));
%     dthetadt(i) = dthetadt(i-1) + 0.5*(d2thetadt2(i-1)+d2thetadt2(i))*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_step.value);
%     delta_v(i) = dthetadt(i)*Aircraft.Geometry.Horizontal.l.value;
%     delta_theta(i) = delta_v(i)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.damping_factor.value/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value);
%     alpha_new_horiz(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz.value(i)  - delta_theta(i);
%  end
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.Results.value = [d2thetadt2, dthetadt, delta_v, delta_theta, alpha_new_horiz];
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.Results.Attributes.contents = "d2thetadt2, dthetadt, delta_v, delta_theta, alpha_new_horiz";
% 
% % CONVERSIONE TO DEGRESS 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.alpha_new_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.Results.value(:,5));
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.alpha_new_deg.Attributes.unit = "deg";
% 
% % CALCULATION OF THE CRITICAL HORIZONTAL TAIL AIRLOADS AS 
% % L_BALANC AT VA + CRITICAL TAIL AIRLOADS DUE TO MANOEUVRE AT VA 
% % FIRSE WE NEED THE CRITICAL TAIL AIRLOADS 
% % LIMIT HORIZONTAL TAIL LOAD 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.DeltaLHorizoTail.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.alpha_new_deg.value(end)*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.a_tail_grad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.S_tail.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.qA.value*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.DeltaLHorizoTail.Attributes.unit = "daN";
% 
% % TOTAL AIRLOADS ACTING ON THE HORIZONTAL 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.TotalLoads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.BalancingLoad.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.DeltaLHorizoTail.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.TotalLoads.Attributes.unit = "daN";
% 
% % ############### VD - MAX UPWARD DEFLECTION - NEGATIVE DEFLECTION - DELTA_E = -0.33 * DELTA_MAX ###############
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Airspeed_at_PointD.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Airspeed_at_PointD.Attributes.unit = "m/sec";
% for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_unit_load_factor.value)   
%     V = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value(i);
%     if abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Airspeed_at_PointD.value - V) < 1e-0
%         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.BalancingLoad.value = (0.5)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Airspeed_at_PointD.value^2)* ... 
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_unit_load_factor.value(end))*(1E-1);   
%     end
% end
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.BalancingLoad.Attributes.unit = "daN";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.qD.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.qD.Attributes.unit = "Pa";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega.value = ((-(1/3)*Aircraft.Geometry.Elevator.max_deflection.value)*Aircraft.Certification.Aerodynamic_data.Horizontal.tau.value)*(1/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.total_time_of_deflection.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega.Attributes.unit = "deg/sec";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega_rad.Attributes.unit = "rad/sec";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.A0.value = (1/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.IY.value)*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.a_tail_rad.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.A0.Attributes.unit = "1/(rad*m*s^2)";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_vector.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz.Attributes.unit = "radians";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_vector.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz_deg.Attributes.unit = "degrees";
% 
% % SOLVING THE DIFFERENTIAL EQUATION 
% d2thetadt2 = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% dthetadt = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1); 
% alpha_new_horiz = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% delta_theta = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% delta_v = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% 
% % OUTPUT 
%  for i = 2:length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_vector.value)
%     d2thetadt2(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.A0.value*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz.value(i) - delta_theta(i-1));
%     dthetadt(i) = dthetadt(i-1) + 0.5*(d2thetadt2(i-1)+d2thetadt2(i))*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_step.value);
%     delta_v(i) = dthetadt(i)*Aircraft.Geometry.Horizontal.l.value;
%     delta_theta(i) = delta_v(i)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.damping_factor.value/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value);
%     alpha_new_horiz(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz.value(i)  - delta_theta(i);
%  end
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.Results.value = [d2thetadt2, dthetadt, delta_v, delta_theta, alpha_new_horiz];
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.Results.Attributes.contents = "d2thetadt2, dthetadt, delta_v, delta_theta, alpha_new_horiz";
% 
% % CONVERSIONE TO DEGRESS 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.alpha_new_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.Results.value(:,5));
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.alpha_new_deg.Attributes.unit = "deg";
% 
% % CALCULATION OF THE CRITICAL HORIZONTAL TAIL AIRLOADS AS 
% % L_BALANC AT VA + CRITICAL TAIL AIRLOADS DUE TO MANOEUVRE AT VA 
% % FIRSE WE NEED THE CRITICAL TAIL AIRLOADS 
% % LIMIT HORIZONTAL TAIL LOAD 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.DeltaLHorizoTail.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.alpha_new_deg.value(end)*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.a_tail_grad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.S_tail.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.qD.value*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.DeltaLHorizoTail.Attributes.unit = "daN";
% 
% % TOTAL AIRLOADS ACTING ON THE HORIZONTAL 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.TotalLoads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.BalancingLoad.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.DeltaLHorizoTail.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.TotalLoads.Attributes.unit = "daN";
% 
% % ############### VD - MAX DOWNWARD DEFLECTION - POSITIVE DEFLECTION - DELTA_E = 0.33 * DELTA_MAX ###############
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.qD.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.qD.Attributes.unit = "Pa";
% 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega.value = (((1/3)*Aircraft.Geometry.Elevator.max_deflection.value)*Aircraft.Certification.Aerodynamic_data.Horizontal.tau.value)*(1/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.total_time_of_deflection.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega.Attributes.unit = "deg/sec";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega_rad.Attributes.unit = "rad/sec";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.A0.value = (1/Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.IY.value)*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.a_tail_rad.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.A0.Attributes.unit = "1/(rad*m*s^2)";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_vector.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz.Attributes.unit = "radians";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.omega_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_vector.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz_deg.Attributes.unit = "degrees";
% 
% % SOLVING THE DIFFERENTIAL EQUATION 
% d2thetadt2 = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% dthetadt = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1); 
% alpha_new_horiz = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% delta_theta = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% delta_v = zeros(length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value), 1);
% 
% % OUTPUT 
%  for i = 2:length(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_vector.value)
%     d2thetadt2(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.A0.value*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz.value(i) - delta_theta(i-1));
%     dthetadt(i) = dthetadt(i-1) + 0.5*(d2thetadt2(i-1)+d2thetadt2(i))*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.time_step.value);
%     delta_v(i) = dthetadt(i)*Aircraft.Geometry.Horizontal.l.value;
%     delta_theta(i) = delta_v(i)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.damping_factor.value/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value);
%     alpha_new_horiz(i) = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz.value(i)  - delta_theta(i);
%  end
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.Results.value = [d2thetadt2, dthetadt, delta_v, delta_theta, alpha_new_horiz];
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.Results.Attributes.contents = "d2thetadt2, dthetadt, delta_v, delta_theta, alpha_new_horiz";
% 
% % CONVERSIONE TO DEGRESS 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.alpha_new_deg.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.Results.value(:,5));
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.alpha_new_deg.Attributes.unit = "deg";
% 
% % CALCULATION OF THE CRITICAL HORIZONTAL TAIL AIRLOADS AS 
% % L_BALANC AT VA + CRITICAL TAIL AIRLOADS DUE TO MANOEUVRE AT VA 
% % FIRSE WE NEED THE CRITICAL TAIL AIRLOADS 
% % LIMIT HORIZONTAL TAIL LOAD 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.DeltaLHorizoTail.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.alpha_new_deg.value(end)*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.a_tail_grad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.S_tail.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.qD.value*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.DeltaLHorizoTail.Attributes.unit = "daN";
% 
% % TOTAL AIRLOADS ACTING ON THE HORIZONTAL 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.TotalLoads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.BalancingLoad.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.DeltaLHorizoTail.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.TotalLoads.Attributes.unit = "daN";
% 
% % Total horizontal tail increment
% disp(" ")
% disp(" ++++ METHOD CS - VLA 423 (c) - DeltaLtail ++++ ")
% disp(" ---------------------------------------------- ")
% Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.DeltaLHorizoTail.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.DeltaLHorizoTail.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.DeltaLHorizoTail.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.DeltaLHorizoTail.value];
% disp(" +++++++++++++++++ Delta Tail loads [daN] +++++++++++++++++ ")
% format = ' %6.6f         %6.6f            %6.6f         %6.6f\n';
% label  = '  Upward defl. VA   Downward defl. VA    Upward defl. VD  Downward defl. VD\n';
% fprintf(label);
% fprintf(format, Total.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % Total horizontal tail increment
% disp(" ++++ METHOD CS - VLA 423 (c) - BTL_1 + DeltaLtail ++++ ")
% disp(" ------------------------------------------------------ ")
% Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.TotalLoads.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.TotalLoads.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.TotalLoads.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.TotalLoads.value];
% disp(" +++++++++++++++++ Total Horizontal Tail loads [daN] +++++++++++++++++ ")
% format = ' %6.6f         %6.6f            %6.6f         %6.6f\n';
% label  = '  Upward defl. VA   Downward defl. VA     Upward defl. VD    Downward defl. VD\n';
% fprintf(label);
% fprintf(format, Total.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % DECISION ABOUT CRITICAL LOAD IN METHOD (c) 
% tl_0 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.TotalLoads.value;
% tl_1 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.TotalLoads.value;
% tl_2 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.TotalLoads.value;
% tl_3 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.TotalLoads.value;
% 
% if (abs(tl_0) > abs(tl_1)) && (abs(tl_0) > abs(tl_2)) && (abs(tl_0) > abs(tl_3))
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Total_critical_loads.value = tl_0;
% elseif (abs(tl_1) > abs(tl_0)) && (abs(tl_1) > abs(tl_2)) && abs(tl_1) > abs(tl_3)
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Total_critical_loads.value = tl_1;
% elseif (abs(tl_2) > abs(tl_0)) && (abs(tl_2) > abs(tl_1)) && abs(tl_2) > abs(tl_3)
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Total_critical_loads.value = tl_2;
% elseif (abs(tl_3) > abs(tl_0)) && (abs(tl_3) > abs(tl_2)) && abs(tl_3) > abs(tl_1)
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Total_critical_loads.value = tl_3;
% end
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Total_critical_loads.Attributes.unit = "daN";
% 
% %% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VA - POINT A OF THE FLIGHT ENVELOPE - CASE 1
% 
% % Normal load n at Point A of the flight envelope
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.nA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.nA.Attributes.unit = "g's";
% 
% % Normal load n at Point A1, indicated in figure 1 of the CS - VLA
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.nA1.value = 1.0; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.nA1.Attributes.unit = "g's";
% 
% % Evaluation of Delta N applied
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.DeltaN.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.nA.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.nA1.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.DeltaN.Attributes.unit = "g's";
% 
% % Evaluation of the product Mg * Delta N
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.DeltaN.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Manoeuvring_weight.Attributes.unit = "N";
% 
% % Evaluation of the ratio X_cg/lt 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.X_cg_lt_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the downwash factor
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Downwash_factor.value = 1 - Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Downwash_factor.Attributes.unit = "Non dimensional";
% 
% % Defining the ratio S_ht/S - surface_ratio
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Surface_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the ratio a_ht/a - Lift_curve_slope_ratio
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Downwash_factor.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Second_term.Attributes.unit = "Non dimensional";
% 
% % Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Third_term.Attributes.unit = "Non dimensional";
% 
% % Evaluating the parenthesis 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Third_term.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.parenthesis.Attributes.unit = "Non dimensional";
% 
% % +++ CRITICAL LOAD - CASE ONE +++ 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.parenthesis.value*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
% 
% %% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VA - POINT A OF THE FLIGHT ENVELOPE - CASE 2
% 
% % Normal load n at Point A of the flight envelope
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.nA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.nA.Attributes.unit = "g's";
% 
% % Normal load n at Point A1, indicated in figure 1 of the CS - VLA
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.nA1.value = 1.0; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.nA1.Attributes.unit = "g's";
% 
% % Evaluation of Delta N applied
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.DeltaN.value = -Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.nA.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.nA1.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.DeltaN.Attributes.unit = "g's";
% 
% % Evaluation of the product Mg * Delta N
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.DeltaN.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Manoeuvring_weight.Attributes.unit = "N";
% 
% % Evaluation of the ratio X_cg/lt 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.X_cg_lt_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the downwash factor
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Downwash_factor.value = 1 - Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Downwash_factor.Attributes.unit = "Non dimensional";
% 
% % Defining the ratio S_ht/S - surface_ratio
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Surface_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the ratio a_ht/a - Lift_curve_slope_ratio
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Downwash_factor.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Second_term.Attributes.unit = "Non dimensional";
% 
% % Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Third_term.Attributes.unit = "Non dimensional";
% 
% % Evaluating the parenthesis 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Third_term.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.parenthesis.Attributes.unit = "Non dimensional";
% 
% % +++ CRITICAL LOAD - CASE TWO +++ 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.parenthesis.value*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
% 
% %% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VA - POINT A OF THE FLIGHT ENVELOPE - CASE 3
% 
% % Normal load n at Point A of the flight envelope
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.nG.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.nG.Attributes.unit = "g's";
% 
% % Normal load n at Point A1, indicated in figure 1 of the CS - VLA
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.nA1.value = 1.0; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.nA1.Attributes.unit = "g's";
% 
% % Evaluation of Delta N applied
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.DeltaN.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.nG.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.nA1.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.DeltaN.Attributes.unit = "g's";
% 
% % Evaluation of the product Mg * Delta N
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.DeltaN.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Manoeuvring_weight.Attributes.unit = "N";
% 
% % Evaluation of the ratio X_cg/lt 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.X_cg_lt_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the downwash factor
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Downwash_factor.value = 1 - Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Downwash_factor.Attributes.unit = "Non dimensional";
% 
% % Defining the ratio S_ht/S - surface_ratio
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Surface_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the ratio a_ht/a - Lift_curve_slope_ratio
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Downwash_factor.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Second_term.Attributes.unit = "Non dimensional";
% 
% % Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Third_term.Attributes.unit = "Non dimensional";
% 
% % Evaluating the parenthesis 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Third_term.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.parenthesis.Attributes.unit = "Non dimensional";
% 
% % +++ CRITICAL LOAD - CASE THREE +++ 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.parenthesis.value*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
% 
% %% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VA - POINT A OF THE FLIGHT ENVELOPE - CASE 4
% 
% % Normal load n at Point A of the flight envelope
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.nG.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.nG.Attributes.unit = "g's";
% 
% % Normal load n at Point A1, indicated in figure 1 of the CS - VLA
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.nA1.value = 1.0; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.nA1.Attributes.unit = "g's";
% 
% % Evaluation of Delta N applied
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.DeltaN.value = - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.nG.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.nA1.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.DeltaN.Attributes.unit = "g's";
% 
% % Evaluation of the product Mg * Delta N
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.DeltaN.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Manoeuvring_weight.Attributes.unit = "N";
% 
% % Evaluation of the ratio X_cg/lt 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.X_cg_lt_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the downwash factor
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Downwash_factor.value = 1 - Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Downwash_factor.Attributes.unit = "Non dimensional";
% 
% % Defining the ratio S_ht/S - surface_ratio
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Surface_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the ratio a_ht/a - Lift_curve_slope_ratio
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Downwash_factor.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Second_term.Attributes.unit = "Non dimensional";
% 
% % Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Third_term.Attributes.unit = "Non dimensional";
% 
% % Evaluating the parenthesis 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Third_term.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.parenthesis.Attributes.unit = "Non dimensional";
% 
% % +++ CRITICAL LOAD - CASE FOUR +++ 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.parenthesis.value*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% %% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VD - POINT D OF THE FLIGHT ENVELOPE - CASE 1
% 
% % Normal load n at Point A of the flight envelope
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.nA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.nA.Attributes.unit = "g's";
% 
% % Normal load n at Point A1, indicated in figure 1 of the CS - VLA
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.nA1.value = 1.0; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.nA1.Attributes.unit = "g's";
% 
% % Evaluation of Delta N applied
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.DeltaN.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.nA.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.nA1.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.DeltaN.Attributes.unit = "g's";
% 
% % Evaluation of the product Mg * Delta N
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.DeltaN.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Manoeuvring_weight.Attributes.unit = "N";
% 
% % Evaluation of the ratio X_cg/lt 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.X_cg_lt_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the downwash factor
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Downwash_factor.value = 1 - Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Downwash_factor.Attributes.unit = "Non dimensional";
% 
% % Defining the ratio S_ht/S - surface_ratio
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Surface_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the ratio a_ht/a - Lift_curve_slope_ratio
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Downwash_factor.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Second_term.Attributes.unit = "Non dimensional";
% 
% % Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Third_term.Attributes.unit = "Non dimensional";
% 
% % Evaluating the parenthesis 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Third_term.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.parenthesis.Attributes.unit = "Non dimensional";
% 
% % +++ CRITICAL LOAD - CASE ONE +++ 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.parenthesis.value*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
% 
% %% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VD - POINT A OF THE FLIGHT ENVELOPE - CASE 2
% 
% % Normal load n at Point A of the flight envelope
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.nA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.nA.Attributes.unit = "g's";
% 
% % Normal load n at Point A1, indicated in figure 1 of the CS - VLA
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.nA1.value = 1.0; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.nA1.Attributes.unit = "g's";
% 
% % Evaluation of Delta N applied
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.DeltaN.value = -Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.nA.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.nA1.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.DeltaN.Attributes.unit = "g's";
% 
% % Evaluation of the product Mg * Delta N
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.DeltaN.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Manoeuvring_weight.Attributes.unit = "N";
% 
% % Evaluation of the ratio X_cg/lt 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.X_cg_lt_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the downwash factor
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Downwash_factor.value = 1 - Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Downwash_factor.Attributes.unit = "Non dimensional";
% 
% % Defining the ratio S_ht/S - surface_ratio
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Surface_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the ratio a_ht/a - Lift_curve_slope_ratio
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Downwash_factor.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Second_term.Attributes.unit = "Non dimensional";
% 
% % Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Third_term.Attributes.unit = "Non dimensional";
% 
% % Evaluating the parenthesis 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Third_term.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.parenthesis.Attributes.unit = "Non dimensional";
% 
% % +++ CRITICAL LOAD - CASE TWO +++ 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.parenthesis.value*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
% 
% %% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VD - POINT E OF THE FLIGHT ENVELOPE - CASE 3
% 
% % Normal load n at Point E of the flight envelope
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.nE.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.nE.Attributes.unit = "g's";
% 
% % Normal load n at Point D1, indicated in figure 1 of the CS - VLA
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.nD1.value = 1.0; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.nD1.Attributes.unit = "g's";
% 
% % Evaluation of Delta N applied
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.DeltaN.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.nE.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.nD1.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.DeltaN.Attributes.unit = "g's";
% 
% % Evaluation of the product Mg * Delta N
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.DeltaN.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Manoeuvring_weight.Attributes.unit = "N";
% 
% % Evaluation of the ratio X_cg/lt 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.X_cg_lt_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the downwash factor
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Downwash_factor.value = 1 - Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Downwash_factor.Attributes.unit = "Non dimensional";
% 
% % Defining the ratio S_ht/S - surface_ratio
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Surface_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the ratio a_ht/a - Lift_curve_slope_ratio
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Downwash_factor.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Second_term.Attributes.unit = "Non dimensional";
% 
% % Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Third_term.Attributes.unit = "Non dimensional";
% 
% % Evaluating the parenthesis 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Third_term.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.parenthesis.Attributes.unit = "Non dimensional";
% 
% % +++ CRITICAL LOAD - CASE THREE +++ 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.parenthesis.value*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
% 
% %% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VD - POINT E OF THE FLIGHT ENVELOPE - CASE 4
% 
% % Normal load n at Point A of the flight envelope
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.nE.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Man_load_factor_nE.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.nE.Attributes.unit = "g's";
% 
% % Normal load n at Point A1, indicated in figure 1 of the CS - VLA
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.nD1.value = 1.0; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.nD1.Attributes.unit = "g's";
% 
% % Evaluation of Delta N applied
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.DeltaN.value = - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.nE.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.nD1.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.DeltaN.Attributes.unit = "g's";
% 
% % Evaluation of the product Mg * Delta N
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Manoeuvring_weight.value = Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.DeltaN.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Manoeuvring_weight.Attributes.unit = "N";
% 
% % Evaluation of the ratio X_cg/lt 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.X_cg_lt_ratio.value = Aircraft.Weight.I_Level.X_cg.value/Aircraft.Geometry.Horizontal.l.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.X_cg_lt_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the downwash factor
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Downwash_factor.value = 1 - Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Downwash_factor.Attributes.unit = "Non dimensional";
% 
% % Defining the ratio S_ht/S - surface_ratio
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Surface_ratio.value = Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Wing.S.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Surface_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the ratio a_ht/a - Lift_curve_slope_ratio
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Lift_curve_slope_ratio.value = Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value/Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";
% 
% % Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Second_term.value = (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Surface_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Lift_curve_slope_ratio.value)*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Downwash_factor.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Second_term.Attributes.unit = "Non dimensional";
% 
% % Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Third_term.value = (1/Aircraft.Weight.I_Level.W_maxTakeOff.value)*0.5*Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value*Aircraft.Geometry.Horizontal.l.value*Aircraft.Geometry.Horizontal.S.value*Aircraft.Certification.Aerodynamic_data.Horizontal.Horizontal_Tail_Normal_Force_Curve_Slope.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Third_term.Attributes.unit = "Non dimensional";
% 
% % Evaluating the parenthesis 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.parenthesis.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.X_cg_lt_ratio.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Second_term.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Third_term.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.parenthesis.Attributes.unit = "Non dimensional";
% 
% % +++ CRITICAL LOAD - CASE FOUR +++ 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Manoeuvring_Critical_Load_Increment.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Manoeuvring_weight.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.parenthesis.value*(1e-1);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
% 
% %% TOTAL AIRLOADS = L_tail + DELTA L_tail
% 
% % AIRSPEED VA - CASE 1 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Total_airloads.Attributes.unit = "daN";
% 
% % AIRSPEED VA - CASE 2 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTS.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Total_airloads.Attributes.unit = "daN";
% 
% % AIRSPEED VA - CASE 3 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Total_airloads.Attributes.unit = "daN";
% 
% % AIRSPEED VA - CASE 4 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Total_airloads.Attributes.unit = "daN";
% 
% % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % AIRSPEED VD - CASE 1 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Total_airloads.Attributes.unit = "daN";
% 
% % AIRSPEED VD - CASE 2 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Total_airloads.Attributes.unit = "daN";
% 
% % AIRSPEED VD - CASE 3 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Total_airloads.Attributes.unit = "daN";
% 
% % AIRSPEED VD - CASE 4 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Total_airloads.Attributes.unit = "daN";
% 
% %% PRINT RESULTS 
% disp(" ")
% disp(" ++++++++++ CS - VLA 423 - METHOD D ++++++++++ ")
% disp(" --------------------------------------------- ")
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Manoeuvring_Critical_Load_Increment.value, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Manoeuvring_Critical_Load_Increment.value, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Manoeuvring_Critical_Load_Increment.value, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Manoeuvring_Critical_Load_Increment.value; ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Manoeuvring_Critical_Load_Increment.value, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Manoeuvring_Critical_Load_Increment.value, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Manoeuvring_Critical_Load_Increment.value, ...
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Manoeuvring_Critical_Load_Increment.value];
% disp(" ++++++++++ Critical Horizontal Tail loads increments [daN] ++++++++++ ")
% format = ' %6.6f          %6.6f          %6.6f          %6.6f\n';
% label  = '  Case 1             Case 2             Case 3              Case 4  \n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % Total horizontal tail increment
% Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Total_airloads.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Total_airloads.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Total_airloads.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Total_airloads.value; ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Total_airloads.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Total_airloads.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Total_airloads.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Total_airloads.value];
% disp(" +++++++++++++++++ Total Horizontal Tail loads [daN] +++++++++++++++++ ")
% format = ' %6.6f          %6.6f         %6.6f          %6.6f\n';
% label  = ' Case 1              Case 2             Case 3             Case 4  \n';
% fprintf(label);
% fprintf(format, Total.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % METHOD CS - VLA (d) CRITICAL LOADS at VA
% tl_0 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_one.Total_airloads.value;
% tl_1 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_two.Total_airloads.value;
% tl_2 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_three.Total_airloads.value;
% tl_3 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.Total_airloads.value; 
% if (abs(tl_0) > abs(tl_1)) && (abs(tl_0) > abs(tl_2)) && (abs(tl_0) > abs(tl_3))
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.Total_critical_loads.value = tl_0;
% elseif (abs(tl_1) > abs(tl_0)) && (abs(tl_1) > abs(tl_2)) && abs(tl_1) > abs(tl_3)
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.Total_critical_loads.value = tl_1;
% elseif (abs(tl_2) > abs(tl_0)) && (abs(tl_2) > abs(tl_1)) && abs(tl_2) > abs(tl_3)
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.Total_critical_loads.value = tl_2;
% elseif (abs(tl_3) > abs(tl_0)) && (abs(tl_3) > abs(tl_2)) && abs(tl_3) > abs(tl_1)
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.Total_critical_loads.value = tl_3;
% end
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.Total_critical_loads.Attributes.unit = "daN";
% 
% % METHOD CS - VLA (d) CRITICAL LOADS at VD
% tl_0 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_one.Total_airloads.value;
% tl_1 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_two.Total_airloads.value;
% tl_2 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_three.Total_airloads.value;
% tl_3 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VD.case_four.Total_airloads.value;
% if (abs(tl_0) > abs(tl_1)) && (abs(tl_0) > abs(tl_2)) && (abs(tl_0) > abs(tl_3))
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.Total_critical_loads.value = tl_0;
% elseif (abs(tl_1) > abs(tl_0)) && (abs(tl_1) > abs(tl_2)) && abs(tl_1) > abs(tl_3)
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.Total_critical_loads.value = tl_1;
% elseif (abs(tl_2) > abs(tl_0)) && (abs(tl_2) > abs(tl_1)) && abs(tl_2) > abs(tl_3)
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.Total_critical_loads.value = tl_2;
% elseif (abs(tl_3) > abs(tl_0)) && (abs(tl_3) > abs(tl_2)) && abs(tl_3) > abs(tl_1)
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.Total_critical_loads.value = tl_3;
% end
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.Total_critical_loads.Attributes.unit = "daN";
% 
% % MAXIMUM VALUE
% max1 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.Total_critical_loads.value;
% max2 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.Total_critical_loads.value;
% if abs(max1) > abs(max2) 
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.Tot_crit_loads.value = max1;
% elseif abs(max2) > abs(max1) 
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.Tot_crit_loads.value = max2;
% end
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.Tot_crit_loads.Attributes.unit = "daN";
% 
% % Horizontal tail loads increments
% Increment = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.Tot_crit_loads.value];
% disp(" ++++++++++ Critical conditions - CS - VLA 423 Method (d) [daN] ++++++++++ ")
% format = ' %6.6f\n';
% label  = 'Critical conditions\n';
% fprintf(label);
% fprintf(format, Increment.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% %% GUST LOAD - CS - VLA 425 
% % CS - VLA 425 Gust loads 
% %  (a) Each horizontal tail surface must be designed for loads resulting
% %      from
% %      (1) Gust velocities specified in CS - VLA 333 (c) with flaps
% %          retracted; and 
% %      (2) Positive and negative gusts of 7.62 m/s nominal intensity at VF
% %          corresponding to the flight conditions specified in CS - VLA 345
% %          (a)(2). 
% %  (b) The average loadings in figure B3 and the distribution of figure B8 
% %      may be used to determine the incremental gust loads for the
% %      requiremenets of subparagraph (a) applied as both up and down
% %      increments for subparagraph (c). 
% %  (c) When determining the total load on the horizontal tail for the
% %      conditions specified in subparagraph (a) of this paragraph, the 
% %      initial balancing tail loads for steady unaccelerated flight at the
% %      pertinent design speeds VF, VC and VD must first be determined. The
% %      incremental tail load resulting from the gusts must be added to the
% %      initial balancing tailload to obtain the total tail load.
% %  (d) In the abscence of a more rational analysis, the incremental tail
% %      load due to the gust, must be computed as follows: 
% %  
% %                   K_g * U_de * V * a_ht * S_ht   /    d epsilon \
% %      Delta L_ht = ---------------------------- * |1 - --------- | 
% %                              16 * 3              \     d alpha  /
% % 
% %      where 
% %      Delta L_ht       = incremental horizontal tail load (daN) 
% %      K_g              = gust alleviation factor defined in CS - VLA 341
% %      U_de             = derived gust velocity (m/s) 
% %      V                = aeroplane equivalent speed (m/s) 
% %      a_ht             = slope of horizonta tail lift curve per radian 
% %      S_ht             = area of horizontal tail (m^2) 
% %      
% %      /    d epsilon \
% %      |1 - --------- | = downwash factor
% %      \     d alpha  /
% %
% % REMIND THAT: 
% % ++++++++++++++++++++
% %       (0.88) * mu_g |
% % K_g = --------------| 
% %         5.3 + mu_g  |
% % ++++++++++++++++++++|
% %          2 * (M/S)  |
% % mu_g = -------------|
% %        rho * MAC * a|
% % ++++++++++++++++++++|
% 
% % MEAN AERODYNAMIC CHORD CALCULATOR OF THE HORIZONTAL TAIL 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.taper_ratio_ht.value = Aircraft.Geometry.Horizontal.ctip.value/Aircraft.Geometry.Horizontal.croot.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.taper_ratio_ht.Attributes.unit = "Non dimensional";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.y_ht.value = linspace(0, 0.5*Aircraft.Geometry.Horizontal.b.value, 1000)';
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.y_ht.Attributes.unit = "meters";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.eta_ht.value = 2*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.y_ht.value/Aircraft.Geometry.Horizontal.b.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.eta_ht.Attributes.unit = "Non dimensional";
% 
% % CHORD DISTRIBUTION
% chord_distr_ht = @(eta) (2*(Aircraft.Geometry.Horizontal.S.value/Aircraft.Geometry.Horizontal.b.value)*(1/(1+Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.taper_ratio_ht.value)))*(1 - ((1-Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.taper_ratio_ht.value)/Aircraft.Geometry.Horizontal.b.value)*abs(eta));
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Chord_distribution_ht.value = chord_distr_ht(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.eta_ht.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Chord_distribution_ht.Attributes.unit = "meters";
% 
% % MEAN AERODYNAMIC CHORD 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.MAC_ht.value = (2/Aircraft.Geometry.Horizontal.S.value)*trapz( Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.y_ht.value, (Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Chord_distribution_ht.value).^2);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.MAC_ht.Attributes.unit = "meters";
% 
% % DATA REQUIRED FOR GUST CALCULATION
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.a_tail_rad.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.a_tail_rad.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.a_tail_rad.Attributes.unit = "1/rad"; 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.mu_g.value = 2*(Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value/Aircraft.Constants.g.value)/(Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.a_tail_rad.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.MAC_ht.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.mu_g.Attributes.unit = "Non dimensional";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.K_g.value = (0.88*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.mu_g.value)/(5.3 + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.mu_g.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.K_g.Attributes.unit = "Non dimensional";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.downwashfactor.value = 1 - Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.downwashfactor.Attributes.unit = "Non dimensional";
% % POINT F
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VF.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VF.Attributes.unit = "daN";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.VF.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.VF.Attributes.unit = "m/s";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Ude_F.value = 7.62;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Ude_F.Attributes.unit = "m/s";
% % POINT C
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VC.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VC.Attributes.unit = "daN";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.VC.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.VC.Attributes.unit = "m/s";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Ude_C.value = 15.20;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Ude_C.Attributes.unit = "m/s";
% % POINT D
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VD.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VD.Attributes.unit = "daN";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.VD.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.VD.Attributes.unit = "m/s";
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Ude_D.value = 7.62;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Ude_D.Attributes.unit = "m/s";
% 
% % FUNCTION TO EVALUEATE DELTA L_ht
% DeltaL_ht = @(V, U_de) (1/16.3)*V*U_de*(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.K_g.value * Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.a_tail_rad.value * Aircraft.Geometry.Horizontal.S.value)*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.downwashfactor.value;
% 
% % DELTA L_HT AT VF
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VF.value = DeltaL_ht(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.VF.value, Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Ude_F.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VF.Attributes.unit = "daN";
% % DELTA L_HT AT VC
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VC.value = DeltaL_ht(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.VC.value, Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Ude_C.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VC.Attributes.unit = "daN";
% % DELTA L_HT AT VD
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VD.value = DeltaL_ht(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.VD.value, Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Ude_D.value);
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VD.Attributes.unit = "daN";
% 
% % PRINT PARTIAL RESULTS
% disp(" ")
% disp(" AIRWORTHINESS RULES: CS - VLA 425 GUST AIRLOADS")
% disp(" --------------------------------- ")
% % Equilibrium balancing loads
% Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VF.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VC.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VD.value];
% disp(" +++++++++++++++++ Total Horizontal Tail loads [daN] +++++++++++++++++ ")
% format = ' %6.6f          %6.6f         %6.6f\n';
% label  = '  L_ht VF             L_ht VC            L_ht VD\n';
% fprintf(label);
% fprintf(format, Total.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % DeltaL_ht horizontal tail increment
% Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VF.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VC.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VD.value];
% disp(" +++++++++++++++++ Delta Horizontal Tail loads [daN] +++++++++++++++++ ")
% % format = '%f          %f         %f\n';
% label  = ' DeltaL_ht VF       DeltaL_ht VC      DeltaL_ht VD\n';
% fprintf(label);
% fprintf(format, Total.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % CALCULATION OF TOTAL GUST LOADS ON THE HORIZONTA TAILPLANE - POSITIVE
% % POINT F
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VF_plus.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VF.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VF.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VF_plus.Attributes.unit = "daN";
% % POINT C
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VC_plus.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VC.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VC.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VC_plus.Attributes.unit = "daN";
% % POINT D
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VD_plus.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VD.value + Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VD.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VD_plus.Attributes.unit = "daN";
% 
% % PRINT TOTAL GUST LOADS
% Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VF_plus.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VC_plus.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VD_plus.value];
% disp(" +++++++++++++++++ Total Horizontal Tail loads - Positive Gust [daN] +++++++++++++++++ ")
% % format = '%f          %f         %f\n';
% label  = '  Gust + L_ht VF    Gust + L_ht VC     Gust + L_ht VD\n';
% fprintf(label);
% fprintf(format, Total.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% % CALCULATION OF TOTAL GUST LOADS ON THE HORIZONTA TAILPLANE - NEGATIVE
% % POINT F
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VF_minus.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VF.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VF.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VF_minus.Attributes.unit = "daN";
% % POINT C
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VC_minus.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VC.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VC.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VC_minus.Attributes.unit = "daN";
% % POINT D
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VD_minus.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VD.value - Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VD.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VD_minus.Attributes.unit = "daN";
% 
% % PRINT TOTAL GUST LOADS
% Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VF_minus.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VC_minus.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VD_minus.value];
% disp(" +++++++++++++++++ Total Horizontal Tail loads - Negative Gust [daN] +++++++++++++++++ ")
% % format = '%f          %f         %f\n';
% label  = '  Gust - L_ht VF      Gust - L_ht VC     Gust - L_ht VD\n';
% fprintf(label);
% fprintf(format, Total.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% %% CRITICAL LOADS FOR SYMMETRICAL CONDITIONS
% % METHOD CS - VLA 423 (a)
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_a.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.Total_critical_load.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_a.Attributes.unit = "daN";
% % METHOD CS - VLA 423 (b)
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_b.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.Total_critical_loads.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_b.Attributes.unit = "daN";
% % METHOD CS - VLA 423 (a)+(b)
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_a_plus_b.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.Total_critical_tail_airloads.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_a_plus_b.Attributes.unit = "daN";
% % METHOD CS - VLA 423 (c)
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_c.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Total_critical_loads.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_c.Attributes.unit = "daN";
% % METHOD CS - VLA 423 (d) 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_d.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.Tot_crit_loads.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_d.Attributes.unit = "daN";
% 
% tl_0 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_a.value;
% tl_1 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_b.value;
% tl_2 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_a_plus_b.value;
% tl_3 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_c.value;
% tl_4 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_d.value;
% if (abs(tl_0) > abs(tl_1)) && (abs(tl_0) > abs(tl_2)) && (abs(tl_0) > abs(tl_3)) && (abs(tl_0) > abs(tl_4))
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.value = tl_0;
% elseif (abs(tl_1) > abs(tl_0)) && (abs(tl_1) > abs(tl_2)) && (abs(tl_1) > abs(tl_3)) && (abs(tl_1) > abs(tl_4))
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.value = tl_1;
% elseif (abs(tl_2) > abs(tl_1)) && (abs(tl_2) > abs(tl_0)) && (abs(tl_2) > abs(tl_3)) && (abs(tl_2) > abs(tl_4))
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.value = tl_2;
% elseif (abs(tl_3) > abs(tl_1)) && (abs(tl_3) > abs(tl_2)) && (abs(tl_3) > abs(tl_0)) && (abs(tl_3) > abs(tl_4))
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.value = tl_3;
% elseif (abs(tl_4) > abs(tl_1)) && (abs(tl_4) > abs(tl_2)) && (abs(tl_4) > abs(tl_3)) && (abs(tl_4) > abs(tl_0))
%     Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.value = tl_4;
% end
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.Attributes.unit = "daN";
% 
% % PRINT TOTAL GUST LOADS
% disp(" ") 
% disp(" Final results from CS - VLA 423 Airworthiness prescriptions ")
% Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.value];
% disp(" +++++++++++++++++ Total Horizontal Tail loads - Maximum Symmetrical load [daN] +++++++++++++++++ ")
% format = ' %f\n';
% label  = '  Maximum symmetrical load\n';
% fprintf(label);
% fprintf(format, Total.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
% 
% %% UNSYMMETRICAL LOADS 
% 
% % CS - VLA 427 Unsymmetrical loads 
% %   (a) Horizontal tail surfaces and their supporting structure must be
% %       designed for unsymmetrical loads arising from yawing and slipstrem
% %       effects in combination with the loads prescribe for the flight
% %       conditions set forth in CS - VLA 421 to 425. 
% %   (b) In the absence of more rational data for aeroplanes that are
% %       conventional in regard to location of the engine, wings, tail
% %       surfaces and fuselage shape 
% %       (1) 100% of the maximum loading from the symmetrical flight
% %           conditions may be assumed on the surface on one side of the
% %           plane of symmetry; and 
% %       (2) The following percentage of that loading must be applied to the
% %           opposite side:
% %           % = 100 -10 * (n - 1)
% %           where n is the specified positive manoeuvring load factor, but
% %           this value may not be more than 80%.
% 
% % NOTA: Il 100% del carico agente sul piano di coda va moltiplicato per 1/2
% %       sul lato del piano di coda interessato dal 100% del carico 
% %       simmetrico. 
% 
% % EVALUATION OF THE PERCENTAGE TO APPLY 
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Percentage_load.value = (100 - 10 * (Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value - 1))*1e-2;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Percentage_load.Attributes.unit = "g's";
% 
% % AIRLOADS ACTING ON THE HORIZ. TAIL IN UNSYMM. CONDITIONS -- FULL LOAD
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Full_load_side.value = 0.5*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Full_load_side.Attributes.unit = "daN";
% % AIRLOADS ACTING ON THE HORIZ. TAIL IN UNSYMM. CONDITIONS -- PARTIAL LOAD
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Partial_load_side.value = 0.5*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Percentage_load.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.value;
% Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Partial_load_side.Attributes.unit = "daN";
% 
% disp(" ")
% disp(" UNSYMMETRICAL LOADS PER CS - VLA 427 ")
% % PRINT TOTAL UNSYMM. LOADS
% Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Full_load_side.value, ...
%          Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Partial_load_side.value];
% disp(" +++++++++++++++++ Total Horizontal Tail loads - Unsymmetrical conditions [daN] +++++++++++++++++ ")
% format = ' %6.6f          %6.6f\n';
% label  = '  Full load side      Partial load side\n';
% fprintf(label);
% fprintf(format, Total.');
% disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
