%% ENGINE TORQUE AND ENGINE MOUNT SUPPORT STRUCTURES 
%  -------------------------------------------------
%  
%  In this calculator we attack the engine structural loads to size its
%  support structures. In the following, the proper paragraphs from CS -
%  VLA airworthiness rules will be reported. 
%
%% CS - VLA 361 ENGINE TORQUE 
%  (a) The engine mount and its supporting structure must be designed for
%      the effects of
%
%      (1) a limit engine torque corresponding to take-off power and
%          propeller speed acting simultaneously with 75 % of the limit
%          loads from flight condition A of CS - VLA 333 (d); 
%
%      (2) the limiting engine torque as specified in CS - VLA 361 (b)
%          acting simultaneously with the limit loads from flight condition
%          A of CS - VLA 333 (d);
%
%  (b) The limit engine torque to be considered under subparagraph (a)(2)
%      of this paragraph must be obtained by multiplying the mean torque
%      for maximum continuous power by a factore determined as follows: 
%
%      (1) For four-stroke engines 
% 
%          (i) 1.33 for engines with five or more cylinders; 
%
%         (ii) 2, 3, 4 or 8 for engines with four, three, two or one
%              cylinders, respectively.
%
%      (2) For two-stroke engines
%
%          (i) 2 for engines with three or more cylinders;
%
%         (ii) 3 or 6 for engines with two or one cylinder respectively.
%
%% CS - VLA 363 SIDE LOAD ON ENGINE MOUNT 
%
%  (a) The engine mount and its supporting structure must be designed for a
%      limit load factor in a lateral direction for the side load on the
%      engine mount of not less than 1.33. 
%
%  (b) The side load prescribed in subparagraph (a) of this paragraph may 
%      assumed to be independent of other flight conditions. 
%
% ------------------------------
% OTHER APPLICABLE CONSIDERATION
% ------------------------------
%  Inside the reference material, other airworthiness rules has been
%  usefule to correctly define all the loads associated with the engine
%  mount structure and its support elements. 
%
%% CS 23.371 GYROSCOPIC AND AERODYNAMIC LOADS 
%
%  (a) Each engine mount and its supporting structure must be designed for
%      the gyroscopic, inertial and aerodynamic loads that result, with the
%      engine (or engines) and propeller (or propellers), if applicable, 
%      at max continuous RPM under either 
%
%      (1) the conditions prescribed in CS 23.351 and 23.423; or 
%
%      (2) All possible combinations of the following: 
%
%          (i) a yaw velocity of 2.5 rad/sec; 
%
%         (ii) a pitch velocity of 1.0 rad/sec; 
%
%        (iii) a normal load factor of 2.5; and 
%
%         (iv) max continuous thrust or power. 
%
%  (b) For aeroplanes approved for aerobatic manoeuvres each engine mount
%      and its supporting structure must meet the requirements of
%      subparagraph (a) and be designed to withstand the load factors
%      expected during combined maximum yaw and pitch velocities.
%
%  (c) For aeroplanes certificated in the commuter category, each engine
%  mount and its supporting structure must meet the requirements of
%  sub-paragraph (a) and the gust conditions specified in CS 23.341. 
%
%% ACCEPTABLE MEANS OF COMPLIANCE 
%  Also, we add two related acceptable means of compliance. 
%% AMC 23.371 METHOD OF EVALUATION OF GYROSCOPIC LOADS 
%  For a two bladed propeller, the maximum gyroscopic couple in N*m is
%  given by 
%                      _______________________________
%                      | 2 * I_p * omega_1 * omega_2 |
%                      -------------------------------
%  For three or more evenly spaced blades, the gyroscopic couples is 
%                          ___________________________
%                          | I_p * omega_1 * omega_2 |
%                          ---------------------------
%  In those formulas, we have 
% 
%  I_p     = Is the polar moment of inertia of the propeller in kg * m^2
%
%  omega_1 = Is the propeller rotation in rad/sec. 
%
%  omega_2 = Is the rate of pitch or yaw in rad/sec.
%
%% AMC 23.371 (a) GYROSCOPIC AND AERODYNAMIC LOADS 
%  The aerodynamic loads specified in CS 23.371 include asymmetric flow
%  through the propeller disc. Experience has shown that the effects of
%  this asymmetric flow on the engine mount and its supporting structure
%  are relatively small and may be discounted, if propellers are installed
%  having diameters of 2.74 m ( 9.00 ft) or less. 

% SWITCH CASE TO ASSIGN THE ADEQUATE CORRECTION FACTOR
switch (Aircraft.Engine.Correction_factor.Attributes.flag1)
    case 'FOUR STROKE' 
        if Aircraft.Engine.Correction_factor.Attributes.flag2 >= 5
            Aircraft.Engine.Correction_factor.value = 1.33;
        elseif Aircraft.Engine.Correction_factor.Attributes.flag2 == 4
            Aircraft.Engine.Correction_factor.value = 2;
        elseif Aircraft.Engine.Correction_factor.Attributes.flag2 == 3
            Aircraft.Engine.Correction_factor.value = 3;
        elseif Aircraft.Engine.Correction_factor.Attributes.flag2 == 2
            Aircraft.Engine.Correction_factor.value = 4;
        elseif Aircraft.Engine.Correction_factor.Attributes.flag2 == 1
            Aircraft.Engine.Correction_factor.value = 8;
        end
    case 'TWO STROKE'
        if Aircraft.Engine.Correction_factor.Attributes.flag2 >= 3
            Aircraft.Engine.Correction_factor.value = 2;
        elseif Aircraft.Engine.Correction_factor.Attributes.flag2 == 2
            Aircraft.Engine.Correction_factor.value = 3;
        elseif Aircraft.Engine.Correction_factor.Attributes.flag2 == 1
            Aircraft.Engine.Correction_factor.value = 6;
        end
end

% LOCAL VARIABLES TO PERFORM CALCULATIONS
correction_factor     = Aircraft.Engine.Correction_factor.value;
takeoff_power         = Aircraft.Engine.Takeoff.Power.value;
takeoff_rpm           = Aircraft.Engine.Takeoff.RPM.value;
max_continous_power   = Aircraft.Engine.Max_Continous.Power.value;
max_continous_rpm     = Aircraft.Engine.Max_Continous.RPM.value;
reduction_ratio       = Aircraft.Engine.Reduction_ratio.value;
g                     = Aircraft.Constants.g.value; 
Limit_side_load       = Aircraft.Engine.Limit_side_load.value;
Gust_limit_load       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.nA1.value;

% ENGINE BLOCK MASS 
Engine_mount_mass       = Aircraft.Engine.Engine_mount_mass.value;
Engine_accessories_mass = Aircraft.Engine.Engine_accessories_mass.value;
Propeller_spinner_mass  = Aircraft.Engine.Propeller_spinner_mass.value;
Engine_block_mass       = Engine_mount_mass + Engine_accessories_mass + Propeller_spinner_mass;
Aircraft.Engine.Engine_block_mass.value           = Engine_block_mass;
Aircraft.Engine.Engine_block_mass.Attributes.unit = "kg";

% CONVERTS MASS TO WEIGHT EXPRESSED IN NEWTON [daN]
Engine_block_weight = Engine_block_mass * g * (1e-1);
Aircraft.Engine.Engine_block_weight.value = Engine_block_weight;
Aircraft.Engine.Engine_block_weight.Attribute.unit = "daN";

% PROPELLER ROTATIONAL SPEED 
prop_rotational_speed = (takeoff_rpm) / (reduction_ratio);
Aircraft.Engine.Propeller_rotational_speed.value = prop_rotational_speed;
Aircraft.Engine.Propeller_rotational_speed.Attributes.unit = "RPM";

% CONVERTS TO ROTATION PER SECONDS 
prop_rotational_speed_sec = prop_rotational_speed * ( 1/60);
Aircraft.Engine.Propeller_rotational_speed_sec.value = prop_rotational_speed_sec; 
Aircraft.Engine.Propeller_rotational_speed_sec.Attributes.unit = "RPS";

% -------------------------------------------------------------------------
% ++++++++++++++++++++++++++++++++ TAKEOFF ++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------
% MEAN ENGINE TORQUE AT MAX TAKEOFF POWER. REMEMBER: COMBINE THIS LOAD WITH
% 0.75 THE LOAD ASSOCIATED WITH POINT A FLIGHT CONDITION
Takeoff_mean_engine_torque = (takeoff_power) * ( (1e3) / (2 * pi * prop_rotational_speed_sec) );
Aircraft.Engine.Takeoff.Mean_torque.value = Takeoff_mean_engine_torque;
Aircraft.Engine.Takeoff.Mean_torque.Attributes.unit = "N * m";

% LIMIT TORQUE AT TAKEOFF POWER
Takeoff_limit_engine_torque = correction_factor * Takeoff_mean_engine_torque;
Aircraft.Engine.Takeoff.Limit_torque.value = Takeoff_limit_engine_torque;
Aircraft.Engine.Takeoff.Limit_torque.Attributes.unit = "N * m";
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++ MAX CONTINOUS +++++++++++++++++++++++++++++
% -------------------------------------------------------------------------
% MEAN ENGINE TORQUE AT MAX CONTINOUS POWER. REMEMBER: COMBINE THIS LOAD WITH
% 1.00 THE LOAD ASSOCIATED WITH POINT A FLIGHT CONDITION
Max_continous_mean_engine_torque = (max_continous_power) * ( (1e3) / (2 * pi * prop_rotational_speed_sec) );
Aircraft.Engine.Max_Continous.Mean_torque.value = Max_continous_mean_engine_torque;
Aircraft.Engine.Max_Continous.Mean_torque.Attributes.unit = "N * m";

% LIMIT TORQUE AT MAX CONTINOUS POWER
Max_continous_limit_engine_torque = correction_factor * Max_continous_mean_engine_torque;
Aircraft.Engine.Max_Continous.Limit_torque.value = Max_continous_limit_engine_torque;
Aircraft.Engine.Max_Continous.Limit_torque.Attributes.unit = "N * m";
% -------------------------------------------------------------------------

% TOTAL SIDE LOAD 
total_side_load = Engine_block_weight * Limit_side_load;
Aircraft.Engine.Total_side_load.value = total_side_load;
Aircraft.Engine.Total_side_load.Attributes.unit = "daN";

% INERTIA LOAD ON ENGINE MOUNT 
% The inertia load is equal to the max limit load factor times the engine
% group weight (in this case point A1 of the diagram)
inertia_load_on_engine_mount = Gust_limit_load * total_side_load; 
Aircraft.Engine.Inertia_load_on_engine_mount.value = inertia_load_on_engine_mount;
Aircraft.Engine.Inertia_load_on_engine_mount.Attributes.unit = "daN";
