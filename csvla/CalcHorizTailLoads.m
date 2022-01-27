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
            qD   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value;
            VA    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
            LHTA  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.value;
            LHTD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value;
            
            % CLHT UNIT LOAD FACTOR
            CLHT_unit_load_factor = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_unit_load_factor.value;            
            
            % LHT FROM 0 TO S
            LHT_from0toS = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_from0toS.value;
            
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
            qD   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value;
            VA   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
            LHTA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.value;
            LHTD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value;
            
            % CLHT UNIT LOAD FACTOR
            CLHT_unit_load_factor = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_unit_load_factor.value;
            
            % LHT FROM 0 TO S
            LHT_from0toS = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_from0toS.value;

            % LOAD FACTOR FROM C TO D
            n_fromCtoD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromCtoD.value;
            V_fromCtoD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromCtoD.value;

            % LOAD FACTOR FROM A TO C
            n_fromAtoC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromA1toC.value ;
            V_fromAtoC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromA1toC.value ;

            % FULL LOAD FACTOR VECTOR 
            full_load_factor_vector = [n_fromAtoC; n_fromCtoD];
            full_airspeed_vector    = [V_fromAtoC; V_fromCtoD];

            % UNIT LOAD FACTOR 
            V_unit_load_factor = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.V_unit_load_factor.value;
            n_unit_load_factor = ones(length(V_unit_load_factor), 1);

        end
    % CASE 2: VA lower than the intercept
    case 'Case 2'
        qA   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
        qD   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value;
        VA   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
        LHTA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTA.value;
        LHTD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value;
        
        % CLHT UNIT LOAD FACTOR
        CLHT_unit_load_factor = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CLHT_unit_load_factor.value;
        
        % LHT FROM 0 TO S
        LHT_from0toS = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_from0toS.value;
        
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
Aircraft.Certification.Aerodynamic_data.Horizontal.tau.Attributes.unit = "Non dimensional";

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
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_step.Attributes.unit = "s";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_interval_num.value = time_interval;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_interval_num.Attributes.unit = "Pure number";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.value = time_vector;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.time_vector.Attributes.unit = "s";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.damping_factor.value = damping_factor; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.damping_factor.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.IY.value = IY;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.IY.Attributes.unit = "kg * m^2"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dvVA.value = dvVA; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.dvVA.Attributes.unit = "deg";

% AMC 23.423 ADVISED TOTAL DEFLECTION TIME INTERVAL 
% BE CAREFUL: AMC 23.423 Suggest the use of various total deflection time
% interval for different aircraft categories (Normal, Utility, Commuter,
% Acrobatic, ...). Pleas, take care of the definition of the time interval.
if (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Aerobatic") && (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag2 == "Stick")
    Aircraft.Geometry.Elevator.total_deflection_time.value = 0.1;
    Aircraft.Geometry.Elevator.total_deflection_time.Attributes.unit = "s";
elseif (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Aerobatic") && (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag2 == "Wheel")
    Aircraft.Geometry.Elevator.total_deflection_time.value = 0.2;
    Aircraft.Geometry.Elevator.total_deflection_time.Attributes.unit = "s";
elseif (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Normal") | (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Utility") | (Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Normal")
    if Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag2 == "Stick"
        Aircraft.Geometry.Elevator.total_deflection_time.value = 0.2;
        Aircraft.Geometry.Elevator.total_deflection_time.Attributes.unit = "s";
    elseif Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag2 == "Wheel"
        Aircraft.Geometry.Elevator.total_deflection_time.value = 0.3;
        Aircraft.Geometry.Elevator.total_deflection_time.Attributes.unit = "s";
    end
elseif Aircraft.Geometry.Elevator.total_deflection_time.Attributes.flag1 == "Normal"
    prmpt = "Enter total defl. time interval --> t_total_defl_time: ";
    Aircraft.Geometry.Elevator.total_deflection_time.value = input(prmpt);
    Aircraft.Geometry.Elevator.total_deflection_time.Attributes = "s";
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
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.omega.Attributes.unit = "deg/s";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.omega_rad.value = omega_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.omega_rad.Attributes.unit = "rad/s";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.A0.value = A0;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.A0.Attributes.unit = "1/(rad*m*s^2)";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.alpha_prime_horiz.value = alpha_prime_horiz_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.alpha_prime_horiz.Attributes.unit = "rad";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.alpha_prime_horiz_deg.value = alpha_prime_horiz_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.alpha_prime_horiz_deg.Attributes.unit = "deg";

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
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.dthetadt.Attributes.unit = "rad/s"; 
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
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.omega.Attributes.unit = "deg/s";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.omega_rad.value = omega_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.omega_rad.Attributes.unit = "rad/s";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.A0.value = A0;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.A0.Attributes.unit = "1/(rad*m*s^2)";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.alpha_prime_horiz.value = alpha_prime_horiz_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.alpha_prime_horiz.Attributes.unit = "rad";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.alpha_prime_horiz_deg.value =alpha_prime_horiz_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.alpha_prime_horiz_deg.Attributes.unit = "deg";

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
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.dthetadt.Attributes.unit = "rad/s"; 
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
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.DeltaLimitLTail.Attributes.cs = " 423(a) ";

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
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.Total_critical_load.Attributes.cs = " 423(a) ";

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
MY_pitch_up = zeros(length(angular_acceleration_pitch_up), 1);

for i = 1:length(angular_acceleration_pitch_up)
    MY_pitch_up(i) = ( IY ) * ( angular_acceleration_pitch_up(i) );
end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.MY.value = MY_pitch_up;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.MY.Attributes.unit = "N*m";

% DELTA TAIL AIRLOADS
delta_tail_airloads_pitch_up = (1e-1)*( - MY_pitch_up ) / ( l_ht);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.delta_tail_airloads.value = delta_tail_airloads_pitch_up; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.delta_tail_airloads.Attributes.unit = "daN";

% BALANCING TAIL AIRLOADS

LHT_unit_load_factor    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.LHT_unit_load_factor.value;
index_va                = dsearchn(V_unit_load_factor, VA);
balancing_tail_airloads_pitch_up = linspace(LHT_unit_load_factor(index_va), LHT_unit_load_factor(end), length(delta_tail_airloads_pitch_up))';
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.balancing_tail_airloads.value = balancing_tail_airloads_pitch_up;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.balancing_tail_airloads.Attributes.unit = "daN";

% TOTAL AIRLOADS ASSOCIATED WITH THE PITCH UP MANOEUVRE
total_tail_airloads_pitch_up = delta_tail_airloads_pitch_up + balancing_tail_airloads_pitch_up;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.total_tail_airloads.value = total_tail_airloads_pitch_up;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.total_tail_airloads.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.total_tail_airloads.Attributes.cs = " 423(b) ";

% CRITICAL TAIL AIRLOADS VALUE 
critical_tail_airloads_pitch_up = -max(abs(total_tail_airloads_pitch_up)); 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.value = critical_tail_airloads_pitch_up;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.Attributes.cs = " 423(b) ";

disp(" ");
disp(" ++++ METHOD CS - VLA 423 (b) ++++");
disp(" ---------------------------------");
disp(" ++++ CHECKED ++++");

% Total horizontal tail increment
Total = critical_tail_airloads_pitch_up;
disp(" +++++++++++++++++ Total Horizontal Tail loads [daN] +++++++++++++++++ ")
format = ' %6.6f\n';
label  = ' At VA\n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

%% CS - VLA 423 - METHOD B - PITCH DOWN

% ANGULAR ACCELERATION CALCULATIONS 
angular_acceleration_pitch_down = -Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.ang_acc.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.ang_acc.value = angular_acceleration_pitch_down;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.ang_acc.Attributes.unit = "rad/sec^2";

% MOMENT MY
MY_pitch_down = -Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.MY.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.MY.value = MY_pitch_down;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.MY.Attributes.unit = "N*m";

% DELTA TAIL AIRLOADS
delta_tail_airloads_pitch_down = ( 1e-1 )*( -MY_pitch_down )/( l_ht );
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.delta_tail_airloads.value = delta_tail_airloads_pitch_down;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.delta_tail_airloads.Attributes.unit = "daN";

% BALANCING TAIL AIRLOADS
index_va = dsearchn(V_unit_load_factor, VA);
balancing_tail_airloads_pitch_down = linspace(LHT_from0toS(index_va), LHT_from0toS(end), length(delta_tail_airloads_pitch_down))';
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.balancing_tail_airloads.value = balancing_tail_airloads_pitch_down;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.balancing_tail_airloads.Attributes.unit = "daN";

% TOTAL AIRLOADS ASSOCIATED WITH THE PITCH UP MANOEUVRE
total_tail_airloads_pitch_down = delta_tail_airloads_pitch_down + balancing_tail_airloads_pitch_down;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.total_tail_airloads.value = total_tail_airloads_pitch_down;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.total_tail_airloads.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.total_tail_airloads.Attributes.cs = " 423(b) ";

% CRITICAL TAIL AIRLOADS VALUE 
critical_tail_airloads_pitch_down = max(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.total_tail_airloads.value); 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.critical_tail_airloads.value = critical_tail_airloads_pitch_down;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.critical_tail_airloads.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.critical_tail_airloads.Attributes.cs = " 423(b) ";

% Total horizontal tail increment
disp(" ");
disp(" ++++ METHOD CS - VLA 423 (b) - PITCH DOWN ++++")
disp(" ----------------------------------------------")
Total = critical_tail_airloads_pitch_down;
disp(" +++++++++++++++++ Horizontal Tail loads [daN] +++++++++++++++++ ")
format = ' %6.6f\n';
label  = ' At VA\n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

if abs(critical_tail_airloads_pitch_up) > abs(critical_tail_airloads_pitch_down)
    Total_critical_loads = critical_tail_airloads_pitch_up;
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.Total_critical_loads.value = Total_critical_loads;
elseif abs(critical_tail_airloads_pitch_down) > abs(critical_tail_airloads_pitch_up)
    Total_critical_loads = critical_tail_airloads_pitch_down;
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.Total_critical_loads.value = Total_critical_loads;
end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.Total_critical_loads.Attributes.unit = "daN";

% Total horizontal tail increment
disp(" ");
disp(" ++++ METHOD CS - VLA 423 (b) - PITCH DOWN ++++")
disp(" ----------------------------------------------")
Total = Total_critical_loads;
disp(" +++++++++++++++++ Total Horizontal Tail loads - Method (b) [daN] +++++++++++++++++ ")
format = ' %6.6f\n';
label  = ' At VA\n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

%% CS - VLA 423 - METHOD (A)+(B) - SWITCH-CASE TO ASSESS THE LOWEST (IN MODULE) LOAD ON THE HORIZ. TAIL 
tl_0 = abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.DeltaLimitLTail.value);
tl_1 = abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.DeltaLimitLTail.value);
tl_2 = abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.value);
tl_3 = abs(Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.critical_tail_airloads.value);
if (tl_0 < tl_1) && (tl_0 < tl_2) && (tl_0 < tl_3)
    % disp(" ++++ CRITICAL CONDITION: METHOD (a) PITCH UP ++++")
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.Attributes.flag = "METHOD CS-VLA 423 (A) PITCH UP";
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_up.DeltaLimitLTail.value;    
elseif (tl_1 < tl_2) && (tl_1 < tl_3) && (tl_1 < tl_0)
    % disp(" ++++ CRITICAL CONDITION: METHOD (a) PITCH DOWN ++++")
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.Attributes.flag = "METHOD CS-VLA 423 (A) PITCH DOWN";
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.pitch_down.DeltaLimitLTail.value;
elseif (tl_2 < tl_1) && (tl_2 < tl_3) && (tl_2 < tl_0)
    % disp(" ++++ CRITICAL CONDITION: METHOD (b) PITCH UP ++++")
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.Attributes.flag = "METHOD CS-VLA 423 (B) PITCH UP";
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_up.critical_tail_airloads.value;
elseif (tl_3 < tl_1) && (tl_3 < tl_2) && (tl_3 < tl_0)
    % disp(" ++++ CRITICAL CONDITION: METHOD (b) PITCH DOWN ++++")
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.Attributes.flag = "METHOD CS-VLA 423 (B) PITCH DOWN";
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.pitch_down.critical_tail_airloads.value;
end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.Attributes.cs = " 423(a)/423(b) ";

% CRITICAL TAIL AIRLOADS
Total_critical_tail_airloads = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.critical_tail_airloads.value + LHTA;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.Total_critical_tail_airloads.value = Total_critical_tail_airloads;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.Total_critical_tail_airloads.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.Total_critical_tail_airloads.Attributes.cs = " 423(a)/423(b) ";

% Total horizontal tail increment
disp(" ");
disp(" ++++ METHOD CS - VLA 423 (a)+(b) ++++")
disp(" ----------------------------------------------")
Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.Total_critical_tail_airloads.value];
disp(" +++++++++++++++++ Total Horizontal Tail loads [daN] +++++++++++++++++ ")
format = ' %6.6f\n';
label  = ' At VA\n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

%% CS - VLA 423 - METHOD (C) 
% ############### VA - MAX UPWARD DEFLECTION - NEGATIVE DELTA ELEVATOR ###############
VA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;

for i = 1:length(CLHT_unit_load_factor)   
    V = V_unit_load_factor(i);
    if abs(VA - V) < 1e-1
         LHT_at_VA = (0.5) * (V^2) * ( S ) * ( rho0 ) * (CLHT_unit_load_factor(i))*(1e-1);   
    end
end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.LHT_at_VA.value = LHT_at_VA;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.LHT_at_VA.Attributes.unit = "daN";

% RATIO L_TAIL / M.A.C. 
L_ratio = l_ht / MAC;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.L_ratio.value =  L_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.L_ratio.Attributes.unit = "Non dimensional";

% RATIO S_TAIL / S_WING 
S_ratio = S_ht / S;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.S_ratio.value = S_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.S_ratio.Attributes.unit = "Non dimensional"; 

% HORIZONTAL TAIL VOLUME RATIO 
Horizontal_Tail_Volume_Ratio = L_ratio * S_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.Horizontal_Tail_Volume_Ratio.value = Horizontal_Tail_Volume_Ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.Horizontal_Tail_Volume_Ratio.Attributes.unit = "Non dimensional";

% PRELIMINARY CALCULATIONS
omega_deg             = ( ( -0.8 * delta_elevator_max ) * tau ) * ( 1 / total_deflection_time );
omega_rad             = deg2rad(omega_deg);
A0                    = ( 1 / IY ) * CLalfa_ht_rad * qA * l_ht * S_ht;
alpha_prime_horiz_rad = omega_rad * time_vector;
alpha_prime_horiz_deg = rad2deg(alpha_prime_horiz_rad);

% STORE INSIDE THE STRUCT VARIABLE
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.omega.value = omega_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.omega.Attributes.unit = "deg/s";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.omega_rad.value = omega_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.omega_rad.Attributes.unit = "rad/s";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.A0.value = A0;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.A0.Attributes.unit = "1/(rad*m*s^2)";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.alpha_prime_horiz.value = alpha_prime_horiz_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.alpha_prime_horiz.Attributes.unit = "rad";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.alpha_prime_horiz_deg.value = alpha_prime_horiz_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.alpha_prime_horiz_deg.Attributes.unit = "deg";

% SOLVING THE DIFFERENTIAL EQUATION 
d2thetadt2      = zeros(length(time_vector), 1);
dthetadt        = zeros(length(time_vector), 1); 
alpha_new_horiz_rad = zeros(length(time_vector), 1);
delta_theta     = zeros(length(time_vector), 1);
delta_v         = zeros(length(time_vector), 1);

% OUTPUT 
 for i = 2:length(time_vector)
    d2thetadt2(i) = A0 * (alpha_prime_horiz_rad(i) - delta_theta(i-1));
    dthetadt(i) = dthetadt(i-1) + 0.5 * ( d2thetadt2(i-1) + d2thetadt2(i)) * (time_step);
    delta_v(i) = dthetadt(i) * l_ht;
    delta_theta(i) = delta_v(i) * ( damping_factor / VA );
    alpha_new_horiz_rad(i) = alpha_prime_horiz_rad(i)  - delta_theta(i);
 end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.Results.value = [d2thetadt2, dthetadt, delta_v, delta_theta, alpha_new_horiz_rad];
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.Results.Attributes.contents = "d2thetadt2, dthetadt, delta_v, delta_theta, alpha_new_horiz";

% CONVERSIONE TO DEGRESS 
alpha_new_horiz_deg = rad2deg(alpha_new_horiz_rad);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.alpha_new_deg.value = alpha_new_horiz_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.alpha_new_deg.Attributes.unit = "deg";

% CALCULATION OF THE CRITICAL HORIZONTAL TAIL AIRLOADS AS 
% L_BALANC AT VA + CRITICAL TAIL AIRLOADS DUE TO MANOEUVRE AT VA 
% FIRSE WE NEED THE CRITICAL TAIL AIRLOADS 
% LIMIT HORIZONTAL TAIL LOAD 
DeltaLHorizoTail_at_VA_upward = alpha_new_horiz_deg(end )* CLalfa_ht_deg * S_ht * qA * (1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.DeltaLHorizoTail.value = DeltaLHorizoTail_at_VA_upward;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.DeltaLHorizoTail.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.DeltaLHorizoTail.Attributes.cs = " 423(c) ";

% TOTAL AIRLOADS ACTING ON THE HORIZONTAL 
Total_loads_at_VA_upward = LHT_at_VA + DeltaLHorizoTail_at_VA_upward;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.TotalLoads.value = Total_loads_at_VA_upward;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.TotalLoads.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.TotalLoads.Attributes.cs = " 423(c) ";

% ############### VA - MAX DOWNWARD DEFLECTION - POSITIVE DELTA ELEVATOR ###############
% POSITIVE MAXIMUM DEFLECTION ANGLE
omega_deg             = ( ( delta_elevator_max ) * tau ) * ( 1 / total_deflection_time);
omega_rad             = deg2rad(omega_deg);
A0                    = (1/ IY )* CLalfa_ht_rad * qA * l_ht * S_ht;
alpha_prime_horiz_rad =  omega_rad * time_vector;
alpha_prime_horiz_deg = rad2deg(alpha_prime_horiz_rad);

Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.omega.value = omega_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.omega.Attributes.unit = "deg/s";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.omega_rad.value = omega_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.omega_rad.Attributes.unit = "rad/s";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.A0.value = A0;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.A0.Attributes.unit = "1/(rad*m*s^2)";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.alpha_prime_horiz.value = alpha_prime_horiz_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.alpha_prime_horiz.Attributes.unit = "rad";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.alpha_prime_horiz_deg.value = alpha_prime_horiz_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.alpha_prime_horiz_deg.Attributes.unit = "deg";

% SOLVING THE DIFFERENTIAL EQUATION 
d2thetadt2          = zeros(length(time_vector), 1);
dthetadt            = zeros(length(time_vector), 1); 
alpha_new_horiz_rad = zeros(length(time_vector), 1);
delta_theta         = zeros(length(time_vector), 1);
delta_v             = zeros(length(time_vector), 1);

% OUTPUT 
 for i = 2:length(time_vector)
    d2thetadt2(i)          = A0 * (alpha_prime_horiz_rad(i) - delta_theta(i-1));
    dthetadt(i)            = dthetadt(i-1) + 0.5 * ( d2thetadt2(i-1) + d2thetadt2(i) ) * (time_step);
    delta_v(i)             = dthetadt(i) * l_ht;
    delta_theta(i)         = delta_v(i) * ( damping_factor / VA );
    alpha_new_horiz_rad(i) = alpha_prime_horiz_rad(i)  - delta_theta(i);
 end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.Results.value = [d2thetadt2, dthetadt, delta_v, delta_theta, alpha_new_horiz_rad];
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.Results.Attributes.contents = "d2thetadt2, dthetadt, delta_v, delta_theta, alpha_new_horiz";

% CONVERSIONE TO DEGRESS 
alpha_new_horiz_deg = rad2deg(alpha_new_horiz_rad);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.alpha_new_deg.value = alpha_new_horiz_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.alpha_new_deg.Attributes.unit = "deg";

% CALCULATION OF THE CRITICAL HORIZONTAL TAIL AIRLOADS AS 
% L_BALANC AT VA + CRITICAL TAIL AIRLOADS DUE TO MANOEUVRE AT VA 
% FIRSE WE NEED THE CRITICAL TAIL AIRLOADS 
% LIMIT HORIZONTAL TAIL LOAD 
DeltaLHorizoTail_at_VA_downward = alpha_new_horiz_deg(end) * CLalfa_ht_deg * S_ht * qA * (1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.DeltaLHorizoTail.value = DeltaLHorizoTail_at_VA_downward;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.DeltaLHorizoTail.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.DeltaLHorizoTail.Attributes.cs = " 423(c) ";

% TOTAL AIRLOADS ACTING ON THE HORIZONTAL 
TotalLoads_at_VA_downward = DeltaLHorizoTail_at_VA_downward + LHTA;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.TotalLoads.value = TotalLoads_at_VA_downward;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.TotalLoads.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.TotalLoads.Attributes.cs = " 423(c) ";

% ############### VD - MAX UPWARD DEFLECTION - NEGATIVE DEFLECTION - DELTA_E = -0.33 * DELTA_MAX ###############
VD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;

for i = 1:length(CLHT_unit_load_factor)   
    V = V_unit_load_factor(i);
    if abs(VD - V) < 1e-1
        LHT_at_VD = (0.5) * (VD^2) * ( S ) * ( rho0 )*( CLHT_unit_load_factor(i) ) * (1e-1);   
    end
end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.BalancingLoad.value = LHT_at_VD;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.BalancingLoad.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.BalancingLoad.Attributes.cs = " 423(c) ";

% PRELIMINARY CALCULATIONS
omega_deg             = ( ( - ( 1/3 ) * delta_elevator_max )* tau )*( 1 / total_deflection_time );
omega_rad             = deg2rad(omega_deg);
A0                    = ( 1 / IY )* CLalfa_ht_rad * qD * l_ht * S_ht; 
alpha_prime_horiz_rad = omega_rad * time_vector;
alpha_prime_horiz_deg = rad2deg(alpha_prime_horiz_rad);

% STORE INSIDE THE STRUCT VARIABLE
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.omega.value = omega_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.omega.Attributes.unit = "deg/s";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.omega_rad.value = omega_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.omega_rad.Attributes.unit = "rad/s";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.A0.value = A0;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.A0.Attributes.unit = "1/(rad*m*s^2)";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.alpha_prime_horiz.value = alpha_prime_horiz_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.alpha_prime_horiz.Attributes.unit = "rad";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.alpha_prime_horiz_deg.value = alpha_prime_horiz_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.alpha_prime_horiz_deg.Attributes.unit = "deg";

% SOLVING THE DIFFERENTIAL EQUATION 
d2thetadt2          = zeros(length(time_vector), 1);
dthetadt            = zeros(length(time_vector), 1); 
alpha_new_horiz_rad = zeros(length(time_vector), 1);
delta_theta         = zeros(length(time_vector), 1);
delta_v             = zeros(length(time_vector), 1);

% OUTPUT 
 for i = 2:length(time_vector)
    d2thetadt2(i)      = A0 * ( alpha_prime_horiz_rad(i) - delta_theta(i-1) );
    dthetadt(i)        = dthetadt(i-1) + 0.5 * ( d2thetadt2(i-1) + d2thetadt2(i) ) * ( time_step );
    delta_v(i)         = dthetadt(i) * l_ht;
    delta_theta(i)     = delta_v(i) * ( damping_factor / VD );
    alpha_new_horiz_rad(i) = alpha_prime_horiz_rad(i)  - delta_theta(i);
 end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.Results.value = [d2thetadt2, dthetadt, delta_v, delta_theta, alpha_new_horiz_rad];
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.Results.Attributes.contents = "d2thetadt2, dthetadt, delta_v, delta_theta, alpha_new_horiz";

% CONVERSIONE TO DEGRESS
alpha_new_horiz_deg = rad2deg(alpha_new_horiz_rad);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.alpha_new_deg.value = alpha_new_horiz_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.alpha_new_deg.Attributes.unit = "deg";

% CALCULATION OF THE CRITICAL HORIZONTAL TAIL AIRLOADS AS 
% L_BALANC AT VA + CRITICAL TAIL AIRLOADS DUE TO MANOEUVRE AT VA 
% FIRSE WE NEED THE CRITICAL TAIL AIRLOADS 
% LIMIT HORIZONTAL TAIL LOAD 
DeltaLHorizoTail_at_VD_upward = alpha_new_horiz_deg(end) * CLalfa_ht_deg * S_ht * qD * (1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.DeltaLHorizoTail.value = DeltaLHorizoTail_at_VD_upward;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.DeltaLHorizoTail.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.DeltaLHorizoTail.Attributes.cs = " 423(c) ";

% TOTAL AIRLOADS ACTING ON THE HORIZONTAL 
TotalLoads_at_VD_upward = LHT_at_VD + DeltaLHorizoTail_at_VD_upward;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.TotalLoads.value = TotalLoads_at_VD_upward;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.TotalLoads.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.TotalLoads.Attributes.cs = " 423(c) ";

% ############### VD - MAX DOWNWARD DEFLECTION - POSITIVE DEFLECTION - DELTA_E = 0.33 * DELTA_MAX ###############
omega_deg = ( ( ( 1 / 3 ) * delta_elevator_max ) * tau ) * ( 1 / total_deflection_time );
omega_rad = deg2rad(omega_deg);
A0        = ( 1 / IY ) * CLalfa_ht_rad * qD * l_ht * S_ht;
alpha_prime_horiz_rad = omega_rad * time_vector;
alpha_prime_horiz_deg = rad2deg(alpha_prime_horiz_rad);

Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.omega.value = omega_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.omega.Attributes.unit = "deg/s";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.omega_rad.value = omega_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.omega_rad.Attributes.unit = "rad/s";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.A0.value = A0;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.A0.Attributes.unit = "1/(rad*m*s^2)";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.alpha_prime_horiz.value = alpha_prime_horiz_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.alpha_prime_horiz.Attributes.unit = "rad";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.alpha_prime_horiz_deg.value = alpha_prime_horiz_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.alpha_prime_horiz_deg.Attributes.unit = "deg";

% SOLVING THE DIFFERENTIAL EQUATION 
d2thetadt2          = zeros(length(time_vector), 1);
dthetadt            = zeros(length(time_vector), 1); 
alpha_new_horiz_rad = zeros(length(time_vector), 1);
delta_theta         = zeros(length(time_vector), 1);
delta_v             = zeros(length(time_vector), 1);

% OUTPUT 
 for i = 2:length(time_vector)
    d2thetadt2(i)          = A0 * ( alpha_prime_horiz_rad(i) - delta_theta(i-1) );
    dthetadt(i)            = dthetadt(i-1) + 0.5 * ( d2thetadt2(i-1) + d2thetadt2(i) ) * ( time_step );
    delta_v(i)             = dthetadt(i) * l_ht;
    delta_theta(i)         = delta_v(i) * ( damping_factor / VD);
    alpha_new_horiz_rad(i) = alpha_prime_horiz_rad(i)  - delta_theta(i);
 end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.Results.value = [d2thetadt2, dthetadt, delta_v, delta_theta, alpha_new_horiz_rad];
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.Results.Attributes.contents = "d2thetadt2, dthetadt, delta_v, delta_theta, alpha_new_horiz";

% CONVERSIONE TO DEGRESS 
alpha_new_horiz_deg = rad2deg(alpha_new_horiz_rad);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.alpha_new_deg.value = alpha_new_horiz_deg;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.alpha_new_deg.Attributes.unit = "deg";

% CALCULATION OF THE CRITICAL HORIZONTAL TAIL AIRLOADS AS 
% L_BALANC AT VA + CRITICAL TAIL AIRLOADS DUE TO MANOEUVRE AT VA 
% FIRSE WE NEED THE CRITICAL TAIL AIRLOADS 
% LIMIT HORIZONTAL TAIL LOAD 
DeltaLHorizoTail_at_VD_downward = alpha_new_horiz_deg(end)* CLalfa_ht_deg * S_ht * qD * (1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.DeltaLHorizoTail.value = DeltaLHorizoTail_at_VD_downward;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.DeltaLHorizoTail.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.DeltaLHorizoTail.Attributes.cs = " 423(c) ";

% TOTAL AIRLOADS ACTING ON THE HORIZONTAL 
TotalLoads_at_VD_downward = LHT_at_VD + DeltaLHorizoTail_at_VD_downward;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.TotalLoads.value = TotalLoads_at_VD_downward;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.TotalLoads.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.TotalLoads.Attributes.cs = " 423(c) ";

% Total horizontal tail increment
disp(" ")
disp(" ++++ METHOD CS - VLA 423 (c) - DeltaLtail ++++ ")
disp(" ---------------------------------------------- ")
Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.DeltaLHorizoTail.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.DeltaLHorizoTail.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.DeltaLHorizoTail.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.DeltaLHorizoTail.value];
disp(" +++++++++++++++++ Delta Tail loads [daN] +++++++++++++++++ ")
format = ' %6.6f         %6.6f            %6.6f         %6.6f\n';
label  = '  Upward defl. VA   Downward defl. VA    Upward defl. VD  Downward defl. VD\n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

% Total horizontal tail increment
disp(" ++++ METHOD CS - VLA 423 (c) - BTL_1 + DeltaLtail ++++ ")
disp(" ------------------------------------------------------ ")
Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.TotalLoads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.TotalLoads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.TotalLoads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.TotalLoads.value];
disp(" +++++++++++++++++ Total Horizontal Tail loads [daN] +++++++++++++++++ ")
format = ' %6.6f         %6.6f            %6.6f         %6.6f\n';
label  = '  Upward defl. VA   Downward defl. VA     Upward defl. VD    Downward defl. VD\n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

% DECISION ABOUT CRITICAL LOAD IN METHOD (c) 
tl_0 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.upward_defl.TotalLoads.value;
tl_1 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.va.downward_defl.TotalLoads.value;
tl_2 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.upward_defl.TotalLoads.value;
tl_3 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.vd.downward_defl.TotalLoads.value;

if (abs(tl_0) > abs(tl_1)) && (abs(tl_0) > abs(tl_2)) && (abs(tl_0) > abs(tl_3))
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Total_critical_loads.value = tl_0;
elseif (abs(tl_1) > abs(tl_0)) && (abs(tl_1) > abs(tl_2)) && abs(tl_1) > abs(tl_3)
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Total_critical_loads.value = tl_1;
elseif (abs(tl_2) > abs(tl_0)) && (abs(tl_2) > abs(tl_1)) && abs(tl_2) > abs(tl_3)
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Total_critical_loads.value = tl_2;
elseif (abs(tl_3) > abs(tl_0)) && (abs(tl_3) > abs(tl_2)) && abs(tl_3) > abs(tl_1)
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Total_critical_loads.value = tl_3;
end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Total_critical_loads.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Total_critical_loads.Attributes.cs = " 423(c) ";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VA - POINT A OF THE FLIGHT ENVELOPE - CASE 1

% SOME CONSTANTS INVOLVED IN THE FOLLOWING CALCULATIONS
Mass            = Aircraft.Weight.I_Level.W_maxTakeOff.value;
g               = Aircraft.Constants.g.value;
Xcg             = Aircraft.Geometry.General.X_cg.value;
DepsilonDalfa   = Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
Downwash_factor = 1.0 - DepsilonDalfa;
S_ratio         = S_ht / S;
CLalfa_rad      = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
a_ratio         = CLalfa_ht_rad / CLalfa_rad;
% -------------------------------------------------------------------------
Delta_load_factor   = nA - 1.0;
DeltaN_times_Weight = Mass * g * Delta_load_factor;
Xcg_l_ht_ratio      = Xcg / l_ht;

% PRODUCT 2: (surface ratio) x (lift curve slope ratio) x (downwash factor)
product2 = S_ratio * a_ratio * Downwash_factor;

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
product3 = ( 1 / Mass ) * 0.5 * rho0 * l_ht * S_ht * CLalfa_ht_rad;

% Evaluating the parenthesis 
parenthesis = Xcg_l_ht_ratio - product2 - product3;

% STORE INSIDE THE STRUCT VARIABLE
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.DeltaN.value = Delta_load_factor;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.DeltaN.Attributes.unit = "g's";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Manoeuvring_weight.value = DeltaN_times_Weight;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Manoeuvring_weight.Attributes.unit = "N";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.X_cg_lt_ratio.value = Xcg_l_ht_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.X_cg_lt_ratio.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Downwash_factor.value = Downwash_factor;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Downwash_factor.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Lift_curve_slope_ratio.value = a_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Second_term.value = product2;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Third_term.value = product3;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.parenthesis.value = parenthesis;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE ONE +++ 
Crit_Load_at_VA_case1 = DeltaN_times_Weight * parenthesis * (1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Manoeuvring_Critical_Load_Increment.value = Crit_Load_at_VA_case1;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Manoeuvring_Critical_Load_Increment.Attributes.cs = " 423(d) ";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VA - POINT A OF THE FLIGHT ENVELOPE - CASE 2

% SOME CONSTANTS INVOLVED IN THE FOLLOWING CALCULATIONS
Mass            = Aircraft.Weight.I_Level.W_maxTakeOff.value;
g               = Aircraft.Constants.g.value;
Xcg             = Aircraft.Geometry.General.X_cg.value;
DepsilonDalfa   = Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
Downwash_factor = 1.0 - DepsilonDalfa;
S_ratio         = S_ht / S;
CLalfa_rad      = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
a_ratio         = CLalfa_ht_rad / CLalfa_rad;
% -------------------------------------------------------------------------
Delta_load_factor   = - nA + 1.0;
DeltaN_times_Weight = Mass * g * Delta_load_factor;
Xcg_l_ht_ratio      = Xcg / l_ht;

% PRODUCT 2: (surface ratio) x (lift curve slope ratio) x (downwash factor)
product2 = S_ratio * a_ratio * Downwash_factor;

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
product3 = ( 1 / Mass ) * 0.5 * rho0 * l_ht * S_ht * CLalfa_ht_rad;

% Evaluating the parenthesis 
parenthesis = Xcg_l_ht_ratio - product2 - product3;

% STORE INSIDE THE STRUCT VARIABLE
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.DeltaN.value = Delta_load_factor;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.DeltaN.Attributes.unit = "g's";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Manoeuvring_weight.value = DeltaN_times_Weight;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Manoeuvring_weight.Attributes.unit = "N"; 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.X_cg_lt_ratio.value = Xcg_l_ht_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.X_cg_lt_ratio.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Downwash_factor.value = Downwash_factor;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Downwash_factor.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Lift_curve_slope_ratio.value = a_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Second_term.value = product2;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Third_term.value = product3;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.parenthesis.value = parenthesis;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE TWO +++ 
Crit_Load_at_VA_case2 = DeltaN_times_Weight * parenthesis * (1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Manoeuvring_Critical_Load_Increment.value = Crit_Load_at_VA_case2;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Manoeuvring_Critical_Load_Increment.Attributes.cs = " 423(d) ";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VA - POINT G OF THE FLIGHT ENVELOPE - CASE 3

% SOME CONSTANTS INVOLVED IN THE FOLLOWING CALCULATIONS
Mass            = Aircraft.Weight.I_Level.W_maxTakeOff.value;
g               = Aircraft.Constants.g.value;
Xcg             = Aircraft.Geometry.General.X_cg.value;
DepsilonDalfa   = Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
Downwash_factor = 1.0 - DepsilonDalfa;
S_ratio         = S_ht / S;
CLalfa_rad      = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
a_ratio         = CLalfa_ht_rad / CLalfa_rad;
% -------------------------------------------------------------------------
Delta_load_factor   = nmin - 1.0;
DeltaN_times_Weight = Mass * g * Delta_load_factor;
Xcg_l_ht_ratio      = Xcg / l_ht;

% PRODUCT 2: (surface ratio) x (lift curve slope ratio) x (downwash factor)
product2 = S_ratio * a_ratio * Downwash_factor;

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
product3 = ( 1 / Mass ) * 0.5 * rho0 * l_ht * S_ht * CLalfa_ht_rad;

% Evaluating the parenthesis 
parenthesis = Xcg_l_ht_ratio - product2 - product3;

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.DeltaN.value = Delta_load_factor;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.DeltaN.Attributes.unit = "g's";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Manoeuvring_weight.value = DeltaN_times_Weight;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Manoeuvring_weight.Attributes.unit = "N";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.X_cg_lt_ratio.value = Xcg_l_ht_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.X_cg_lt_ratio.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Downwash_factor.value = Downwash_factor;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Downwash_factor.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Lift_curve_slope_ratio.value = a_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Second_term.value = product2;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Third_term.value = product3;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.parenthesis.value = parenthesis;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE THREE +++ 
Crit_Load_at_VA_case3 = DeltaN_times_Weight * parenthesis * (1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Manoeuvring_Critical_Load_Increment.value = Crit_Load_at_VA_case3;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Manoeuvring_Critical_Load_Increment.Attributes.cs = " 423(d) ";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VA - POINT A OF THE FLIGHT ENVELOPE - CASE 4

% SOME CONSTANTS INVOLVED IN THE FOLLOWING CALCULATIONS
Mass            = Aircraft.Weight.I_Level.W_maxTakeOff.value;
g               = Aircraft.Constants.g.value;
Xcg             = Aircraft.Geometry.General.X_cg.value;
DepsilonDalfa   = Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
Downwash_factor = 1.0 - DepsilonDalfa;
S_ratio         = S_ht / S;
CLalfa_rad      = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
a_ratio         = CLalfa_ht_rad / CLalfa_rad;
% -------------------------------------------------------------------------
Delta_load_factor   = - nmin + 1.0;
DeltaN_times_Weight = Mass * g * Delta_load_factor;
Xcg_l_ht_ratio      = Xcg / l_ht;

% PRODUCT 2: (surface ratio) x (lift curve slope ratio) x (downwash factor)
product2 = S_ratio * a_ratio * Downwash_factor;

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
product3 = ( 1 / Mass ) * 0.5 * rho0 * l_ht * S_ht * CLalfa_ht_rad;

% Evaluating the parenthesis 
parenthesis = Xcg_l_ht_ratio - product2 - product3;

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.DeltaN.value = Delta_load_factor;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.VA.case_four.DeltaN.Attributes.unit = "g's";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Manoeuvring_weight.value = DeltaN_times_Weight;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Manoeuvring_weight.Attributes.unit = "N";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.X_cg_lt_ratio.value = Xcg_l_ht_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.X_cg_lt_ratio.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Downwash_factor.value = Downwash_factor;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Downwash_factor.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Lift_curve_slope_ratio.value = a_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Second_term.value = product2;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Third_term.value = product3;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.parenthesis.value = parenthesis;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE FOUR +++ 
Crit_Load_at_VA_case4 = DeltaN_times_Weight * parenthesis * (1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Manoeuvring_Critical_Load_Increment.value = Crit_Load_at_VA_case4;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Manoeuvring_Critical_Load_Increment.Attributes.cs = " 423(d) ";
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VD - POINT D OF THE FLIGHT ENVELOPE - CASE 1

% SOME CONSTANTS INVOLVED IN THE FOLLOWING CALCULATIONS
Mass            = Aircraft.Weight.I_Level.W_maxTakeOff.value;
g               = Aircraft.Constants.g.value;
Xcg             = Aircraft.Geometry.General.X_cg.value;
DepsilonDalfa   = Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
Downwash_factor = 1.0 - DepsilonDalfa;
S_ratio         = S_ht / S;
CLalfa_rad      = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
a_ratio         = CLalfa_ht_rad / CLalfa_rad;
% -------------------------------------------------------------------------
Delta_load_factor   = nmax - 1.0;
DeltaN_times_Weight = Mass * g * Delta_load_factor;
Xcg_l_ht_ratio      = Xcg / l_ht;

% PRODUCT 2: (surface ratio) x (lift curve slope ratio) x (downwash factor)
product2 = S_ratio * a_ratio * Downwash_factor;

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
product3 = ( 1 / Mass ) * 0.5 * rho0 * l_ht * S_ht * CLalfa_ht_rad;

% Evaluating the parenthesis 
parenthesis = Xcg_l_ht_ratio - product2 - product3;

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.DeltaN.value = Delta_load_factor;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.DeltaN.Attributes.unit = "g's";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Manoeuvring_weight.value = DeltaN_times_Weight;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Manoeuvring_weight.Attributes.unit = "N";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.X_cg_lt_ratio.value = Xcg_l_ht_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.X_cg_lt_ratio.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Downwash_factor.value = Downwash_factor;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Downwash_factor.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Lift_curve_slope_ratio.value = a_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Second_term.value = product2;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Third_term.value = product3;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.parenthesis.value = parenthesis;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE ONE +++ 
Crit_Load_at_VD_case1 = DeltaN_times_Weight * parenthesis * (1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Manoeuvring_Critical_Load_Increment.value = Crit_Load_at_VD_case1;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Manoeuvring_Critical_Load_Increment.Attributes.cs = " 423(d) ";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VD - POINT A OF THE FLIGHT ENVELOPE - CASE 2

% SOME CONSTANTS INVOLVED IN THE FOLLOWING CALCULATIONS
Mass            = Aircraft.Weight.I_Level.W_maxTakeOff.value;
g               = Aircraft.Constants.g.value;
Xcg             = Aircraft.Geometry.General.X_cg.value;
DepsilonDalfa   = Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
Downwash_factor = 1.0 - DepsilonDalfa;
S_ratio         = S_ht / S;
CLalfa_rad      = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
a_ratio         = CLalfa_ht_rad / CLalfa_rad;
% -------------------------------------------------------------------------
Delta_load_factor   = - nmax + 1.0;
DeltaN_times_Weight = Mass * g * Delta_load_factor;
Xcg_l_ht_ratio      = Xcg / l_ht;

% PRODUCT 2: (surface ratio) x (lift curve slope ratio) x (downwash factor)
product2 = S_ratio * a_ratio * Downwash_factor;

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
product3 = ( 1 / Mass ) * 0.5 * rho0 * l_ht * S_ht * CLalfa_ht_rad;

% Evaluating the parenthesis 
parenthesis = Xcg_l_ht_ratio - product2 - product3;

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.DeltaN.value = Delta_load_factor;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.DeltaN.Attributes.unit = "g's";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Manoeuvring_weight.value = DeltaN_times_Weight;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Manoeuvring_weight.Attributes.unit = "N";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.X_cg_lt_ratio.value = Xcg_l_ht_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.X_cg_lt_ratio.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Downwash_factor.value = Downwash_factor;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Downwash_factor.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Lift_curve_slope_ratio.value = a_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Second_term.value = product2;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Third_term.value = product3;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.parenthesis.value = parenthesis;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE TWO +++ 
Crit_Load_at_VD_case2 = DeltaN_times_Weight * parenthesis * (1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Manoeuvring_Critical_Load_Increment.value = Crit_Load_at_VD_case2;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Manoeuvring_Critical_Load_Increment.Attributes.cs = " 423(d) ";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VD - POINT E OF THE FLIGHT ENVELOPE - CASE 3

% SOME CONSTANTS INVOLVED IN THE FOLLOWING CALCULATIONS
Mass            = Aircraft.Weight.I_Level.W_maxTakeOff.value;
g               = Aircraft.Constants.g.value;
Xcg             = Aircraft.Geometry.General.X_cg.value;
DepsilonDalfa   = Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
Downwash_factor = 1.0 - DepsilonDalfa;
S_ratio         = S_ht / S;
CLalfa_rad      = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
a_ratio         = CLalfa_ht_rad / CLalfa_rad;
% -------------------------------------------------------------------------
Delta_load_factor   = nmin - 1.0;
DeltaN_times_Weight = Mass * g * Delta_load_factor;
Xcg_l_ht_ratio      = Xcg / l_ht;

% PRODUCT 2: (surface ratio) x (lift curve slope ratio) x (downwash factor)
product2 = S_ratio * a_ratio * Downwash_factor;

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
product3 = ( 1 / Mass ) * 0.5 * rho0 * l_ht * S_ht * CLalfa_ht_rad;

% Evaluating the parenthesis 
parenthesis = Xcg_l_ht_ratio - product2 - product3;

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.DeltaN.value = Delta_load_factor;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.DeltaN.Attributes.unit = "g's";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Manoeuvring_weight.value = DeltaN_times_Weight;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Manoeuvring_weight.Attributes.unit = "N";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.X_cg_lt_ratio.value = Xcg_l_ht_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.X_cg_lt_ratio.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Downwash_factor.value = Downwash_factor;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Downwash_factor.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Lift_curve_slope_ratio.value = a_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Second_term.value = product2;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Third_term.value = product3;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.parenthesis.value = parenthesis;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE THREE +++ 
Crit_Load_at_VD_case3 = DeltaN_times_Weight * parenthesis * (1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Manoeuvring_Critical_Load_Increment.value = Crit_Load_at_VD_case3;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Manoeuvring_Critical_Load_Increment.Attributes.cs = " 423(d) ";

%% CS - VLA 423 - METHOD D - MANOEUVRING AIRSPEED VD - POINT E OF THE FLIGHT ENVELOPE - CASE 4


% SOME CONSTANTS INVOLVED IN THE FOLLOWING CALCULATIONS
Mass            = Aircraft.Weight.I_Level.W_maxTakeOff.value;
g               = Aircraft.Constants.g.value;
Xcg             = Aircraft.Geometry.General.X_cg.value;
DepsilonDalfa   = Aircraft.Certification.Aerodynamic_data.Horizontal.DepsilonDalpha.value;
Downwash_factor = 1.0 - DepsilonDalfa;
S_ratio         = S_ht / S;
CLalfa_rad      = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
a_ratio         = CLalfa_ht_rad / CLalfa_rad;
% -------------------------------------------------------------------------
Delta_load_factor   = - nmin + 1.0;
DeltaN_times_Weight = Mass * g * Delta_load_factor;
Xcg_l_ht_ratio      = Xcg / l_ht;

% PRODUCT 2: (surface ratio) x (lift curve slope ratio) x (downwash factor)
product2 = S_ratio * a_ratio * Downwash_factor;

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
product3 = ( 1 / Mass ) * 0.5 * rho0 * l_ht * S_ht * CLalfa_ht_rad;

% Evaluating the parenthesis 
parenthesis = Xcg_l_ht_ratio - product2 - product3;

% Evaluation of Delta N applied
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.DeltaN.value = Delta_load_factor;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.DeltaN.Attributes.unit = "g's";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Manoeuvring_weight.value = DeltaN_times_Weight;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Manoeuvring_weight.Attributes.unit = "N";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.X_cg_lt_ratio.value = Xcg_l_ht_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.X_cg_lt_ratio.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Downwash_factor.value = Downwash_factor;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Downwash_factor.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Lift_curve_slope_ratio.value = a_ratio;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Lift_curve_slope_ratio.Attributes.unit = "Non dimensional";

% Defining the product (surface ratio) x (lift curve slope ratio) x (downwash factor) - Second term
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Second_term.value = product2;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Second_term.Attributes.unit = "Non dimensional";

% Defining the product 0.5 * (1/M) * rho0 * lt * S_ht * a_ht
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Third_term.value = product3;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Third_term.Attributes.unit = "Non dimensional";

% Evaluating the parenthesis 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.parenthesis.value = parenthesis;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.parenthesis.Attributes.unit = "Non dimensional";

% +++ CRITICAL LOAD - CASE FOUR +++ 
Crit_Load_at_VD_case4 = DeltaN_times_Weight * parenthesis * (1e-1);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Manoeuvring_Critical_Load_Increment.value = Crit_Load_at_VD_case4;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Manoeuvring_Critical_Load_Increment.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Manoeuvring_Critical_Load_Increment.Attributes.cs = " 423(d) ";

%% TOTAL AIRLOADS = L_tail + DELTA L_tail

% AIRSPEED VA - CASE 1 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Total_airloads.Attributes.unit = "daN";

% AIRSPEED VA - CASE 2 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTS.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Total_airloads.Attributes.unit = "daN";

% AIRSPEED VA - CASE 3 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.LHTA1.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Total_airloads.Attributes.unit = "daN";

% AIRSPEED VA - CASE 4 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTG.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Total_airloads.Attributes.unit = "daN";

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% AIRSPEED VD - CASE 1 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Total_airloads.Attributes.unit = "daN";

% AIRSPEED VD - CASE 2 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Total_airloads.Attributes.unit = "daN";

% AIRSPEED VD - CASE 3 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTE.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Total_airloads.Attributes.unit = "daN";

% AIRSPEED VD - CASE 4 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Total_airloads.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Manoeuvring_Critical_Load_Increment.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Total_airloads.Attributes.unit = "daN";

%% PRINT RESULTS 
disp(" ")
disp(" ++++++++++ CS - VLA 423 - METHOD D ++++++++++ ")
disp(" --------------------------------------------- ")
% Horizontal tail loads increments
Increment = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Manoeuvring_Critical_Load_Increment.value, ...
            Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Manoeuvring_Critical_Load_Increment.value, ...
            Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Manoeuvring_Critical_Load_Increment.value, ...
            Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Manoeuvring_Critical_Load_Increment.value; ...
            Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Manoeuvring_Critical_Load_Increment.value, ...
            Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Manoeuvring_Critical_Load_Increment.value, ...
            Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Manoeuvring_Critical_Load_Increment.value, ...
            Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Manoeuvring_Critical_Load_Increment.value];
disp(" ++++++++++ Critical Horizontal Tail loads increments [daN] ++++++++++ ")
format = ' %6.6f          %6.6f          %6.6f          %6.6f\n';
label  = '  Case 1             Case 2             Case 3              Case 4  \n';
fprintf(label);
fprintf(format, Increment.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

% Total horizontal tail increment
Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Total_airloads.value; ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Total_airloads.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Total_airloads.value];
disp(" +++++++++++++++++ Total Horizontal Tail loads [daN] +++++++++++++++++ ")
format = ' %6.6f          %6.6f         %6.6f          %6.6f\n';
label  = ' Case 1              Case 2             Case 3             Case 4  \n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

% METHOD CS - VLA (d) CRITICAL LOADS at VA
tl_0 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_one.Total_airloads.value;
tl_1 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_two.Total_airloads.value;
tl_2 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_three.Total_airloads.value;
tl_3 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.va.case_four.Total_airloads.value; 
if (abs(tl_0) > abs(tl_1)) && (abs(tl_0) > abs(tl_2)) && (abs(tl_0) > abs(tl_3))
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.va.Total_critical_loads.value = tl_0;
elseif (abs(tl_1) > abs(tl_0)) && (abs(tl_1) > abs(tl_2)) && abs(tl_1) > abs(tl_3)
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.va.Total_critical_loads.value = tl_1;
elseif (abs(tl_2) > abs(tl_0)) && (abs(tl_2) > abs(tl_1)) && abs(tl_2) > abs(tl_3)
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.va.Total_critical_loads.value = tl_2;
elseif (abs(tl_3) > abs(tl_0)) && (abs(tl_3) > abs(tl_2)) && abs(tl_3) > abs(tl_1)
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.va.Total_critical_loads.value = tl_3;
end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.va.Total_critical_loads.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.va.Total_critical_loads.Attributes.cs = " 423(d) ";

% METHOD CS - VLA (d) CRITICAL LOADS at VD
tl_0 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_one.Total_airloads.value;
tl_1 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_two.Total_airloads.value;
tl_2 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_three.Total_airloads.value;
tl_3 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.vd.case_four.Total_airloads.value;
if (abs(tl_0) > abs(tl_1)) && (abs(tl_0) > abs(tl_2)) && (abs(tl_0) > abs(tl_3))
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.vd.Total_critical_loads.value = tl_0;
elseif (abs(tl_1) > abs(tl_0)) && (abs(tl_1) > abs(tl_2)) && abs(tl_1) > abs(tl_3)
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.vd.Total_critical_loads.value = tl_1;
elseif (abs(tl_2) > abs(tl_0)) && (abs(tl_2) > abs(tl_1)) && abs(tl_2) > abs(tl_3)
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.vd.Total_critical_loads.value = tl_2;
elseif (abs(tl_3) > abs(tl_0)) && (abs(tl_3) > abs(tl_2)) && abs(tl_3) > abs(tl_1)
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.vd.Total_critical_loads.value = tl_3;
end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.vd.Total_critical_loads.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.vd.Total_critical_loads.Attributes.cs = " 423(d) ";

% MAXIMUM VALUE
max1 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.va.Total_critical_loads.value;
max2 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.vd.Total_critical_loads.value;
if abs(max1) > abs(max2) 
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.Tot_crit_loads.value = max1;
elseif abs(max2) > abs(max1) 
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.Tot_crit_loads.value = max2;
end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.Tot_crit_loads.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.Tot_crit_loads.Attributes.cs = " 423(d) ";

% Horizontal tail loads increments
Increment = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.Tot_crit_loads.value];
disp(" ++++++++++ Critical conditions - CS - VLA 423 Method (d) [daN] ++++++++++ ")
format = ' %6.6f\n';
label  = 'Critical conditions\n';
fprintf(label);
fprintf(format, Increment.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.general_rules.cs = " 423 "; 
% COMMENTO IMPORTANTE PER LE PROSSIME SESSIONI: 
% In questo momento, per come sono stati definiti i fattori di carico, non
% vi sono differenze tra il caso V = VA ed il caso V = VD; questo può
% essere facilmente cambiato intervenendo sulla variabile DeltaN,
% scegliendo altri punti dell'inviluppo di volo.

%% GUST LOAD - CS - VLA 425 
% CS - VLA 425 Gust loads 
%  (a) Each horizontal tail surface must be designed for loads resulting
%      from
%      (1) Gust velocities specified in CS - VLA 333 (c) with flaps
%          retracted; and 
%      (2) Positive and negative gusts of 7.62 m/s nominal intensity at VF
%          corresponding to the flight conditions specified in CS - VLA 345
%          (a)(2). 
%  (b) The average loadings in figure B3 and the distribution of figure B8 
%      may be used to determine the incremental gust loads for the
%      requiremenets of subparagraph (a) applied as both up and down
%      increments for subparagraph (c). 
%  (c) When determining the total load on the horizontal tail for the
%      conditions specified in subparagraph (a) of this paragraph, the 
%      initial balancing tail loads for steady unaccelerated flight at the
%      pertinent design speeds VF, VC and VD must first be determined. The
%      incremental tail load resulting from the gusts must be added to the
%      initial balancing tailload to obtain the total tail load.
%  (d) In the abscence of a more rational analysis, the incremental tail
%      load due to the gust, must be computed as follows: 
%  
%                   K_g * U_de * V * a_ht * S_ht   /    d epsilon \
%      Delta L_ht = ---------------------------- * |1 - --------- | 
%                              16 * 3              \     d alpha  /
% 
%      where 
%      Delta L_ht       = incremental horizontal tail load (daN) 
%      K_g              = gust alleviation factor defined in CS - VLA 341
%      U_de             = derived gust velocity (m/s) 
%      V                = aeroplane equivalent speed (m/s) 
%      a_ht             = slope of horizonta tail lift curve per radian 
%      S_ht             = area of horizontal tail (m^2) 
%      
%      /    d epsilon \
%      |1 - --------- | = downwash factor
%      \     d alpha  /
%
% REMIND THAT: 
% ++++++++++++++++++++
%       (0.88) * mu_g |
% K_g = --------------| 
%         5.3 + mu_g  |
% ++++++++++++++++++++|
%          2 * (M/S)  |
% mu_g = -------------|
%        rho * MAC * a|
% ++++++++++++++++++++|
numb           = 1e3;
b              = Aircraft.Geometry.Horizontal.b.value;
b_half         = 0.5*Aircraft.Geometry.Horizontal.b.value;
ctip           = Aircraft.Geometry.Horizontal.ctip.value; 
croot          = Aircraft.Geometry.Horizontal.croot.value;
taper_ratio_ht = ctip / croot; 
y_ht           = linspace(0, b_half, numb)';
eta_ht         = 2 * ( y_ht / b);
WS             = Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value;

% MEAN AERODYNAMIC CHORD CALCULATOR OF THE HORIZONTAL TAIL 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.taper_ratio_ht.value = taper_ratio_ht;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.taper_ratio_ht.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.y_ht.value = y_ht;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.y_ht.Attributes.unit = "m";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.eta_ht.value = eta_ht;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.eta_ht.Attributes.unit = "Non dimensional";

% CHORD DISTRIBUTION
chord_distr_ht        = @(eta) ( 2 * ( S_ht / b )*( 1 /(1 + taper_ratio_ht) ) ) * ( 1 - ( ( 1 - taper_ratio_ht ) / b ) * abs(eta) );
Chord_distribution_ht = chord_distr_ht(eta_ht);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Chord_distribution_ht.value = Chord_distribution_ht;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Chord_distribution_ht.Attributes.unit = "m";

% MEAN AERODYNAMIC CHORD 
MAC_ht = (2 / S_ht) * trapz( y_ht, ( Chord_distribution_ht ).^2);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.MAC_ht.value = MAC_ht;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.MAC_ht.Attributes.unit = "m";

% DATA REQUIRED FOR GUST CALCULATION
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.a_tail_rad.value = CLalfa_ht_rad;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.a_tail_rad.Attributes.unit = "1/rad"; 

% MU GUST
mu_g = 2 * ( WS / g ) / ( rho0 * CLalfa_ht_rad * MAC_ht);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.mu_g.value = mu_g;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.mu_g.Attributes.unit = "Non dimensional";

% GUST ALLEVIATION FACTOR
K_g = ( 0.88 * mu_g ) / ( 5.3 + mu_g );                                                 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.K_g.value = K_g;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.K_g.Attributes.unit = "Non dimensional";

% DOWNWASH FACTOR
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.downwashfactor.value = Downwash_factor;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.downwashfactor.Attributes.unit = "Non dimensional";

% POINT F
L_ht_VF = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTF.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VF.value = L_ht_VF;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VF.Attributes.unit = "daN";
VF = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.VF.value = VF;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.VF.Attributes.unit = "m/s";
Ude_F = 7.62;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Ude_F.value = Ude_F;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Ude_F.Attributes.unit = "m/s";

% POINT C
L_ht_VC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTC.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VC.value = L_ht_VC;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VC.Attributes.unit = "daN";
VC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.VC.value = VC;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.VC.Attributes.unit = "m/s";
Ude_C = 15.20;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Ude_C.value = Ude_C;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Ude_C.Attributes.unit = "m/s";

% POINT D
L_ht_VD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTD.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VD.value = L_ht_VD;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.L_ht_VD.Attributes.unit = "daN";
VD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.VD.value = VD;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.VD.Attributes.unit = "m/s";
Ude_D = 7.62;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Ude_D.value = Ude_D;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Ude_D.Attributes.unit = "m/s";

% FUNCTION TO EVALUEATE DELTA L_ht
DeltaL_ht = @(V, U_de) ( 1 / 16.3 ) * V * U_de *( K_g * CLalfa_ht_rad * S_ht )* Downwash_factor;

% DELTA L_HT AT VF
DeltaL_ht_VF = DeltaL_ht( VF, Ude_F );
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VF.value = DeltaL_ht_VF;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VF.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VF.Attributes.cs = " 425 ";
% DELTA L_HT AT VC
DeltaL_ht_VC = DeltaL_ht( VC, Ude_C);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VC.value = DeltaL_ht_VC;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VC.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VC.Attributes.cs = " 425 ";
% DELTA L_HT AT VD
DeltaL_ht_VD = DeltaL_ht( VD, Ude_D);
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VD.value = DeltaL_ht_VD;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VD.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.DeltaL_ht_VD.Attributes.cs = " 425 ";

% PRINT PARTIAL RESULTS
disp(" ")
disp(" AIRWORTHINESS RULES: CS - VLA 425 GUST AIRLOADS")
disp(" --------------------------------- ")
% Equilibrium balancing loads
Total = [L_ht_VF, ...
         L_ht_VC, ...
         L_ht_VD];
disp(" +++++++++++++++++ Total Horizontal Tail loads [daN] +++++++++++++++++ ")
format = ' %6.6f          %6.6f         %6.6f\n';
label  = '  L_ht VF             L_ht VC            L_ht VD\n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

% DeltaL_ht horizontal tail increment
Total = [DeltaL_ht_VF, ...
         DeltaL_ht_VC, ...
         DeltaL_ht_VD];
disp(" +++++++++++++++++ Delta Horizontal Tail loads [daN] +++++++++++++++++ ")
% format = '%f          %f         %f\n';
label  = ' DeltaL_ht VF       DeltaL_ht VC      DeltaL_ht VD\n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

% CALCULATION OF TOTAL GUST LOADS ON THE HORIZONTA TAILPLANE - POSITIVE
% POINT F
Total_gust_at_VF_plus = L_ht_VF + DeltaL_ht_VF;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VF_plus.value = Total_gust_at_VF_plus;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VF_plus.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VF_plus.Attributes.cs = " 425 ";
% POINT C
Total_gust_at_VC_plus = L_ht_VC + DeltaL_ht_VC;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VC_plus.value = Total_gust_at_VC_plus;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VC_plus.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VC_plus.Attributes.cs = " 425 ";
% POINT D
Total_gust_at_VD_plus = L_ht_VD + DeltaL_ht_VD;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VD_plus.value = Total_gust_at_VD_plus;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VD_plus.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VD_plus.Attributes.cs = " 425 ";

% PRINT TOTAL GUST LOADS
Total = [Total_gust_at_VF_plus, ...
         Total_gust_at_VC_plus, ...
         Total_gust_at_VD_plus];
disp(" +++++++++++++++++ Total Horizontal Tail loads - Positive Gust [daN] +++++++++++++++++ ")
% format = '%f          %f         %f\n';
label  = '  Gust + L_ht VF    Gust + L_ht VC     Gust + L_ht VD\n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

% CALCULATION OF TOTAL GUST LOADS ON THE HORIZONTA TAILPLANE - NEGATIVE
% POINT F
Total_gust_at_VF_minus = L_ht_VF - DeltaL_ht_VF;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VF_minus.value = Total_gust_at_VF_minus;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VF_minus.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VF_minus.Attributes.cs = " 425 ";
% POINT C
Total_gust_at_VC_minus = L_ht_VC - DeltaL_ht_VC;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VC_minus.value = Total_gust_at_VC_minus;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VC_minus.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VC_minus.Attributes.cs = " 425 ";
% POINT D
Total_gust_at_VD_minus = L_ht_VD - DeltaL_ht_VD;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VD_minus.value = Total_gust_at_VD_minus;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VD_minus.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Gustloads.Total_gust_at_VD_minus.Attributes.cs = " 425 ";

% PRINT TOTAL GUST LOADS
Total = [Total_gust_at_VF_minus, ...
         Total_gust_at_VC_minus, ...
         Total_gust_at_VD_minus];
disp(" +++++++++++++++++ Total Horizontal Tail loads - Negative Gust [daN] +++++++++++++++++ ")
% format = '%f          %f         %f\n';
label  = '  Gust - L_ht VF      Gust - L_ht VC     Gust - L_ht VD\n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

%% CRITICAL LOADS FOR SYMMETRICAL CONDITIONS
% METHOD CS - VLA 423 (a)
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_a.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a.Total_critical_load.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_a.Attributes.unit = "daN";
% METHOD CS - VLA 423 (b)
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_b.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_b.Total_critical_loads.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_b.Attributes.unit = "daN";
% METHOD CS - VLA 423 (a)+(b)
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_a_plus_b.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_a_plus_b.Total_critical_tail_airloads.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_a_plus_b.Attributes.unit = "daN";
% METHOD CS - VLA 423 (c)
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_c.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_c.Total_critical_loads.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_c.Attributes.unit = "daN";
% METHOD CS - VLA 423 (d) 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_d.value = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.Method_d.critical_case.Tot_crit_loads.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_d.Attributes.unit = "daN";

tl_0 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_a.value;
tl_1 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_b.value;
tl_2 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_a_plus_b.value;
tl_3 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_c.value;
tl_4 = Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Method_d.value;
if (abs(tl_0) > abs(tl_1)) && (abs(tl_0) > abs(tl_2)) && (abs(tl_0) > abs(tl_3)) && (abs(tl_0) > abs(tl_4))
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.value = tl_0;
elseif (abs(tl_1) > abs(tl_0)) && (abs(tl_1) > abs(tl_2)) && (abs(tl_1) > abs(tl_3)) && (abs(tl_1) > abs(tl_4))
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.value = tl_1;
elseif (abs(tl_2) > abs(tl_1)) && (abs(tl_2) > abs(tl_0)) && (abs(tl_2) > abs(tl_3)) && (abs(tl_2) > abs(tl_4))
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.value = tl_2;
elseif (abs(tl_3) > abs(tl_1)) && (abs(tl_3) > abs(tl_2)) && (abs(tl_3) > abs(tl_0)) && (abs(tl_3) > abs(tl_4))
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.value = tl_3;
elseif (abs(tl_4) > abs(tl_1)) && (abs(tl_4) > abs(tl_2)) && (abs(tl_4) > abs(tl_3)) && (abs(tl_4) > abs(tl_0))
    Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.value = tl_4;
end
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.Attributes.unit = "daN";

% PRINT TOTAL GUST LOADS
disp(" ") 
disp(" Final results from CS - VLA 423 Airworthiness prescriptions ")
Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.value];
disp(" +++++++++++++++++ Total Horizontal Tail loads - Maximum Symmetrical load [daN] +++++++++++++++++ ")
format = ' %f\n';
label  = '  Maximum symmetrical load\n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")

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
% 
%% UNSYMMETRICAL LOADS 

% CS - VLA 427 Unsymmetrical loads 
%   (a) Horizontal tail surfaces and their supporting structure must be
%       designed for unsymmetrical loads arising from yawing and slipstrem
%       effects in combination with the loads prescribe for the flight
%       conditions set forth in CS - VLA 421 to 425. 
%   (b) In the absence of more rational data for aeroplanes that are
%       conventional in regard to location of the engine, wings, tail
%       surfaces and fuselage shape 
%       (1) 100% of the maximum loading from the symmetrical flight
%           conditions may be assumed on the surface on one side of the
%           plane of symmetry; and 
%       (2) The following percentage of that loading must be applied to the
%           opposite side:
%           % = 100 -10 * (n - 1)
%           where n is the specified positive manoeuvring load factor, but
%           this value may not be more than 80%.

% NOTA: Il 100% del carico agente sul piano di coda va moltiplicato per 1/2
%       sul lato del piano di coda interessato dal 100% del carico 
%       simmetrico. 

% EVALUATION OF THE PERCENTAGE TO APPLY 
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Percentage_load.value = (100 - 10 * (Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value - 1))*1e-2;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Percentage_load.Attributes.unit = "g's";

% AIRLOADS ACTING ON THE HORIZ. TAIL IN UNSYMM. CONDITIONS -- FULL LOAD
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Full_load_side.value = 0.5*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Full_load_side.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Full_load_side.Attributes.cs = " 427 ";
% AIRLOADS ACTING ON THE HORIZ. TAIL IN UNSYMM. CONDITIONS -- PARTIAL LOAD
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Partial_load_side.value = 0.5*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Percentage_load.value*Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.CriticalLoads.Maximum.value;
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Partial_load_side.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Partial_load_side.Attributes.cs = " 427 ";

disp(" ")
disp(" UNSYMMETRICAL LOADS PER CS - VLA 427 ")
% PRINT TOTAL UNSYMM. LOADS
Total = [Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Full_load_side.value, ...
         Aircraft.Certification.Regulation.SubpartC.HorizontalTailLoads.UnsymmetricalLoads.Partial_load_side.value];
disp(" +++++++++++++++++ Total Horizontal Tail loads - Unsymmetrical conditions [daN] +++++++++++++++++ ")
format = ' %6.6f          %6.6f\n';
label  = '  Full load side      Partial load side\n';
fprintf(label);
fprintf(format, Total.');
disp(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ")


