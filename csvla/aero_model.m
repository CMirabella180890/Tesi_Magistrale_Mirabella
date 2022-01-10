classdef aero_model
    %% aero_model - Tools for aerodynamic data interpolation 
    %   This class files are mainly devoted to define functions and methods
    %   that are able to conveniently interpolate all the aerodynamic data
    %   necessary to evluate horizontal tail balancing loads. Following
    %   airworthiness rules, it is possible to define balancing loads as
    %   follow (CS - VLA 421 Balancing loads, pag. 49/190): 
    %   
    %   (a) A horizontal tail balancing load is a load necessary to
    %       maintain equilibrium in any specified flight condition with no
    %       pitching acceleration. 
    %   (b) Horizontal tail surfaces must be designed for the balancing
    %       loads occurring at any point ON THE LIMIT MANOEUVRING ENVELOPE
    %       and IN THE FLAP CONDITIONS specified in CS - VLA 345. 
    %    
    %   This definition can be interpreted as a requirements about the
    %   structural sizing and design of the horizontal empennage in any
    %   possible flight condition that shows to be critical for structural
    %   integrity of the aircraft in static, non - accelerated flight
    %   conditions. Effects on structural sizing of accelerated flight are
    %   not taken into account in this particular case. 
%% Property 1    
    properties
        aero_interp
        % Property defined to correctly interpolate aerodynamic data.  
        % Users can choose the most convenient way to interpolate the data.
    end
    
    methods
        % -----------------------------------------------------------------
        function obj = CLMax(obj, a, b, c)
            % METHOD 1
            % Maximum Lift Coefficient = CLMax(a, b, c)
            %   
            %   This method is aimed to evaluate a usable maximum lift
            %   coefficient, able to take into account the fact that wings
            %   must compensate for the extra lift produced by the
            %   horizontal tail to achieve trimmed flight in any possible
            %   conditions. The function construct a second order
            %   polynomial to evaluate CLMax at high angle of attack.
            %   Applicability of this particular function is strictly
            %   limited to the case of the CS - VLA certification process
            %   for a MTOW = 100 kg drone examined as a test case.
            %   The equation used is the following:
            %   
            %   CLMax = a*(alpha)^2 + b*(alpha) + c
            %
            %   INPUTS 
            %   a     --> Coefficient of the second order term inside the used
            %             equation
            %   b     --> Coefficient of the linear term inside the used
            %             equation
            %   c     --> Constant coefficient inside the used equation 
            %   
            %   OUTPUTS 
            %   CLMax --> Estimated maximum lift coefficient
            %   
            %   NOTE: It must be noticed that is not strictly necessary to
            %         define a wide range of AoA for the purpose of this
            %         function. This is done here for illustration
            %         purposes.
            
            AOA_aux = linspace(-20.0, 20.0, 250);
            vect    = a.*(AOA_aux.^2) + b.*AOA_aux + c;
            clmax   = max(vect);
            obj     = clmax;
        end
        % -----------------------------------------------------------------
        function obj = cd_calc(obj, CD0, CL, AR, e, k1, k2)
            % METHOD 2
            % Drag Coefficient = CD(CD0, CL, AR, e, k1, k2)
            % 
            %   This function is able to evaluate the drag coefficient of
            %   the aircraft in the prescribed flight condition. It must be
            %   noticed that 
            %   - to apply this function, certain aerodynamic and
            %     geometrical parameters of the aircraft must be known; and
            %   - the applied formula for the non linear part can be not so
            %     accurate and checks about the results must be performed.
            %   The function can be recursively applied to store drag
            %   coefficient values inside a vector. For instance, if a
            %   vector which contains CL values is defined, is possible to
            %   write
            %   
            %   CD = zeros(length(CL), 1);
            %   for i = 1:length(CL) 
            %          CD(i) = cd_calc(CD0, CL(i), AR, e, k1, k2); 
            %   end
            %   To approximate the drag coefficient, the following scheme
            %   is used:
            %   IF CL < 1.2                         (CL)^2
            %   CONVENTIONAL DRAG POLAR = CD0 + ---------------
            %                                     PI * AR * e
            %   IF CL >= 1.2
            %   The conventional drag polar is corrected by two terms
            %   LINEAR CORRECTION    --> k1*(CL - 1.2)
            %   QUADRATIC CORRECTION --> k2*(CL - 1.2)^2
            %
            %   INPUTS
            %   CD0 --> Zero - lift drag coefficient 
            %   CL  --> Input lift coefficient, related to a particular
            %           flight condition which is completely defined in
            %           terms of airspeed V and load factor n (aircraft
            %           weight is a known quantity). It must be remembered
            %           that the CL used here is relative to a wing - body
            %           configuration (wing - fuselage CL, without
            %           horizontal empennage)
            %   AR  --> Aircraft aspect ratio 
            %   e   --> Oswald's efficiency factor
            %   k1  --> Coefficient of the linear term of the formula used
            %           to approximate more exactly the drag coefficient
            %           inside the non - linear portion of Lift vs AoA
            %           curve.
            %   k2  --> Coefficient of the quadratic term of the formula
            %           used to approximate more exactly the drag 
            %           coefficient inside the non - linear portion of 
            %           Lift vs AoA curve.
            %   OUTPUTS 
            %   CD  --> Evaluated drag coefficient
            
            cdcalc = 0.0;
            if CL < 1.2 
                cdcalc = CD0 + ((CL^2)/(pi*AR*e));
            elseif CL >= 1.2 
                cdcalc = CD0 + ((CL^2)/(pi*AR*e)) + k2*(CL - 1.2)^2 + ...
                    k1*(CL - 1.2);
            end
            obj = cdcalc;
        end
        % -----------------------------------------------------------------
        function obj = alpha_calc(obj, CL, CL0, CL_star, CLalpha, p)
            % METHOD 3
            % Angle of attack = ALPHA(obj, CL, CL0, CL_star, CLalpha, a, b, c)
            % 
            %   This function is able to evaluate the angle of attack of
            %   the aircraft in the prescribed flight condition. It must be
            %   noticed that 
            %   - to apply this function, certain aerodynamic and
            %     geometrical parameters of the aircraft must be known; and
            %   - the applied formula for the non linear part can be not so
            %     accurate and checks about the results must be performed.
            %   The function can be recursively applied to store AoA
            %   values inside a vector. For instance, if a vector which 
            %   contains CL values is defined, is possible to write
            %   
            %   alpha = zeros(length(CL), 1);
            %   for i = 1:length(CL) 
            %          alpha(i) = alpha_calc(obj, CL, CL0, CL_star, CLalpha, a, b, c); 
            %   end
            %   
            %   The inputs must be checked to the correct application of
            %   the function, especially with respect to the normal -
            %   force/lift - force curve slope CLalpha. This must be
            %   defined in [1/deg] if the angle of attack have to be 
            %   expressed in degrees. 
            %
            %   INPUTS
            %   CL      --> Input lift coefficient, related to a particular
            %               flight condition which is completely defined in
            %               terms of airspeed V and load factor n (aircraft
            %               weight is a known quantity). It must be 
            %               remembered that the CL used here is relative to
            %               a wing - body configuration (wing - fuselage CL
            %               without horizontal empennage)
            %   CL0     --> Lift coefficient when the aircraft is flying at
            %               a zero angle of attack
            %   CL_star --> Lift coefficient at which the lift curve stop
            %               to be linear (end of linearity for the lift
            %               curve)
            %   CLalpha --> Lift - curve slope coefficient; it must be
            %               expressed in [1/deg] if we want to define the
            %               angle of attack in degrees
            %   a       --> Coefficient a of the interpolation scheme
            %   b       --> Coefficient b of the interpolation scheme
            %   c       --> Coefficient c of the interpolation scheme
            %   OUTPUTS 
            %   ALPHA   --> Evaluated angle of attack at the prescribed
            %               flight condition
            %   To approximate the angle of attack, the following scheme
            %   is used:
            %   IF CL < CL_STAR                        
            %   LINEAR RELATIONSHIP BETWEEN LIFT COEFFICIENT AND ALPHA 
            %
            %            ( CL - CL0)
            %   ALPHA = ------------- 
            %              CLalpha
            %
            %   IF CL >= CL_STAR
            %   INTERPOLATION FORMULA 
            %   
            %          -b + sqrt(b^2 - 4*(a*(c - CL)))
            %  ALPHA = -------------------------------
            %                       2*a
            
            alphacalc = zeros(length(CL), 1);
            for i = 1:length(CL)
                if CL(i) < CL_star
                    x = (CL(i) - CL0)/(CLalpha);
                    alphacalc(i) = x;
                elseif CL(i) >= CL_star - 0.03
                    % p = polyfit(CL_supp, alpha_supp, 2);
                    x = polyval(p, CL(i));
                    % x = interp1(CL_supp, alpha_supp, CL(i), 'linear', 'extrap');
                    alphacalc(i) = x;
                end
            end
           % if CL < CL_star 
               % x         = (CL - CL0)/(CLalpha);
               % alphacalc = x;
           % elseif CL >= CL_star
               % SECOND DEGREE INTERPOLATION  
               % ++++++++++++++++++++++++++++++++++++ 
               % absolute   = abs(c - CL);          |
               % under_root = (b^2) - 4*a*absolute; |
               % x          = -b + sqrt(under_root);|
               % y          = 2*a;                  |
               % z          = x/y;                  |
               % alphacalc  = z;                    |
               % +++++++++++++++++++++++++++++++++++|
               
               % MATLAB INTERPOLATION 
               % +++++++++++++++++++++++++++++++++++
               
               % +++++++++++++++++++++++++++++++++++
            % end
            % Convert alpha in degrees
            % alphacalc = rad2deg(alphacalc);
            obj       = alphacalc;
        end
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        function obj = alpha_calc_lin(obj, CL, CL0, CLalpha)
            % METHOD 3
            % Angle of attack = ALPHA(obj, CL, CL0, CL_star, CLalpha, a, b, c)
            % 
            %   This function is able to evaluate the angle of attack of
            %   the aircraft in the prescribed flight condition. It must be
            %   noticed that 
            %   - to apply this function, certain aerodynamic and
            %     geometrical parameters of the aircraft must be known; and
            %   - the applied formula for the non linear part can be not so
            %     accurate and checks about the results must be performed.
            %   The function can be recursively applied to store AoA
            %   values inside a vector. For instance, if a vector which 
            %   contains CL values is defined, is possible to write
            %   
            %   alpha = zeros(length(CL), 1);
            %   for i = 1:length(CL) 
            %          alpha(i) = alpha_calc(obj, CL, CL0, CLalpha); 
            %   end
            %   
            %   The inputs must be checked to the correct application of
            %   the function, especially with respect to the normal -
            %   force/lift - force curve slope CLalpha. This must be
            %   defined in [1/deg] if the angle of attack have to be 
            %   expressed in degrees. 
            %
            %   INPUTS
            %   CL      --> Input lift coefficient, related to a particular
            %               flight condition which is completely defined in
            %               terms of airspeed V and load factor n (aircraft
            %               weight is a known quantity). It must be 
            %               remembered that the CL used here is relative to
            %               a wing - body configuration (wing - fuselage CL
            %               without horizontal empennage)
            %   CL0     --> Lift coefficient when the aircraft is flying at
            %               a zero angle of attack
            %   CLalpha --> Lift - curve slope coefficient; it must be
            %               expressed in [1/deg] if we want to define the
            %               angle of attack in degrees
            %   OUTPUTS 
            %   ALPHA   --> Evaluated angle of attack at the prescribed
            %               flight condition
            %   To approximate the angle of attack, the following scheme
            %   is used:
            %   IF CL < CL_STAR                        
            %   LINEAR RELATIONSHIP BETWEEN LIFT COEFFICIENT AND ALPHA 
            %
            %            ( CL - CL0)
            %   ALPHA = ------------- 
            %              CLalpha
            %
            
            alphacalc = 0.0;
            alphacalc = (CL - CL0)/(CLalpha);
            obj       = alphacalc;
        end
        % -----------------------------------------------------------------
        function obj = CL_calc(obj, n, W, g, V, rho, S)
            
            % METHOD 4
            % Lift coefficient = CL(obj, n, W, g, V, rho, S)
            % 
            %  This function is able to calculate the lift coefficient 
            %  related to a particular flight condition which is completely
            %  defined in terms of airspeed V and load factor n (aircraft
            %  weight and geometry are known quantities). It must be 
            %  noticed that to apply this function, certain aerodynamic and
            %  geometrical parameters of the aircraft must be known.
            %  This function can be used to have a convenient starting
            %  point to an iterative procedure aimed to evaluate the amount
            %  of lift produced by the horizontal empennage, which is
            %  often a negative lift. This results in an increase of the
            %  lift produced by the main lifting surface (i.e., the wing)
            %  and this fact must be taken into account. Also, the lift
            %  coefficient cannot be greater than CLMAX, evaluated via
            %  METHOD 1 contained inside this class file.
            %  
            %  INPUTS 
            %  n  --> Load factor relative to the chosen flight condition.
            %         The value could be relative to a critical point of the
            %         V - n diagram, such as the maximum dive speed point or
            %         the manoeuvring point (point D and point A on the 
            %         V - n diagram
            %  W  --> Max Takeoff Gross Weight of the chosen aircraft
            %         measured in kilograms
            %  g  --> Gravity acceleration in SI units [m/s^2]
            %  q  --> Dynamic pressure calculated via (0.5)*(V^2)*(rho),
            %         where
            %         V   --> Airspeed at the chosen flight condition
            %         rho --> Standard atmosphere density at selected 
            %                 flight condition
            %  S  --> Wing surface
            %
            %  OUTPUTS
            %  CL --> Lift coefficient evaluated at the prescribed flight
            %         condition and aircraft configuration
            %
            %  The function calculates the lift coefficient as follows: 
            %
            %        n * W * g
            %  CL = -----------
            %          q * S
            %  
            %  It must be remembered that the CL here defined is relative 
            %  to a wing - body configuration (wing - fuselage CL without 
            %  horizontal empennage.
            
            q    = 0.5*rho*(mpower(V, 2));
            MTOW = g*W;
            obj  = (n*MTOW)/(q*S);
        end
    end
%% Property 2
    properties
        pitch_equilibrium
        % Property defined to collect all the contributes related to 
        % aircraft pitching moment equilibrium. 
    end
    % ---------------------------------------------------------------------
    methods
        % -----------------------------------------------------------------
        function obj = CLWB_contrib(obj, CL, alpha, xac, xcg, bcg, MAC)
            
            % METHOD 1
            % Lift contribution to pitching mom. = CLWBcontrib(obj, CL, alpha, xac, xcg, bcg, MAC)
            % 
            %  This function is able to calculate the lift contribution
            %  related to the pitching moment equilibrium during a
            %  particular flight condition, defined in terms of airspeed V
            %  and load factor n. The equilibrium is taken about the
            %  aircraft center of gravity. This function take into account 
            %  all the terms that can be derived from simplified
            %  representation of the equilibrium condition; projection of 
            %  lift in a two dimensional fashion is performed, and all the 
            %  contributions are summed. 
            %  To be more specific, two contribution to the pitching moment
            %  are relate to the lift produced by the main lifting surface: 
            %  - a vertical component, which produce pitching moment if its
            %    arm with respect to the aircraft center of gravity is
            %    different from zero; 
            %  - an horizontal component, which produce pitching moment if
            %    its arm with respect to the aircraft center of gravity is
            %    different from zero. It must be remembered that, for a
            %    proper designed aircraft, the lift force horizontal
            %    component is relatively small.
            %  
            %  INPUTS 
            %  CL    --> Input lift coefficient, related to a particular
            %            flight condition which is completely defined in
            %            terms of airspeed V and load factor n (aircraft
            %            weight is a known quantity). It must be 
            %            remembered that the CL used here is relative to
            %            a wing - body configuration (wing - fuselage CL
            %            without horizontal empennage)
            %  alpha --> Input angle of attack, used in this function to
            %            correctly project along wind axes the lift force 
            %            components in its vertical and horizontal
            %            component
            %  xac   --> Position along the X - Axis of the aircraft 
            %            aerodynamic center; usually xac = 0.25 for a
            %            typical, subsonic airfoil when Xac is divided by
            %            the Mean Aerodynamic Chord.
            %  xcg   --> Position along the X - Axis of the aircraft center
            %            of gravity; for typical overall geom. arrangement
            %            and loading condition, the xcg could be anywhere
            %            between 0.2 and 0.5 when Xcg is divided by the
            %            Mean Aerodynamic chord
            %  bcg   --> Dimensional, vertical arm of the lift force
            %            horizontal component; this must be divided by the
            %            Mean Aerodynamic chord before to apply the
            %            equation for the pitching mom. equilibrium
            %  MAC   --> Mean Aerodynamic chord
            %
            %  OUTPUTS
            %  LIFT PITCHING MOMENT CONTRIBUTION
            %  It's the main lifting surface contribution to pitching
            %  moment and its lift.
            %
            %  The function calculates the lift coefficient contribution to 
            %  the pitching moment contribution is:
            %
            %                                                          bcg
            %  CLWB_contr = CL*cos(alpha)*(xac - xcg) - CL*sin(alpha)* ---
            %                                                          MAC  
            
            x    = bcg/MAC;
            y    = xac - xcg; 
            ter1 = cos(alpha);
            ter2 = sin(alpha);
            
            obj  = CL*ter1*y - CL*ter2*(x);
        end
        % -----------------------------------------------------------------
        function obj = CDWB_contrib(obj, CD, alpha, xac, xcg, bcg, MAC)
            
            % METHOD 2
            % Drag contribution to pitching mom. = CDWBcontrib(obj, CD, ...
            %                                    alpha, xac, xcg, bcg, MAC)
            % 
            %  This function is able to calculate the drag contribution
            %  related to the pitching moment equilibrium during a
            %  particular flight condition, defined in terms of airspeed V
            %  and load factor n. The equilibrium is taken about the
            %  aircraft center of gravity. This function take into account 
            %  all the terms that can be derived from simplified
            %  representation of the equilibrium condition; projection of 
            %  drag in a two dimensional fashion is performed, and all the 
            %  contributions are summed. 
            %  To be more specific, two contribution to the pitching moment
            %  are relate to the drag produced by the main lifting surface: 
            %  - a vertical component, which produce pitching moment if its
            %    arm with respect to the aircraft center of gravity is
            %    different from zero. It must be remembered that, for a
            %    typical design, the drag force vertical component is 
            %    relatively small;
            %  - an horizontal component, which produce pitching moment if
            %    its arm with respect to the aircraft center of gravity is
            %    different from zero. 
            %  
            %  INPUTS 
            %  CD    --> Input drag coefficient, related to a particular
            %            flight condition which is completely defined in
            %            terms of airspeed V and load factor n (aircraft
            %            weight is a known quantity). It must be 
            %            remembered that the CL used here is relative to
            %            a wing - body configuration (wing - fuselage CL
            %            without horizontal empennage)
            %  alpha --> Input angle of attack, used in this function to
            %            correctly project along wind axes the lift force 
            %            components in its vertical and horizontal
            %            component
            %  xac   --> Position along the X - Axis of the aircraft 
            %            aerodynamic center; usually xac = 0.25 for a
            %            typical, subsonic airfoil when Xac is divided by
            %            the Mean Aerodynamic Chord.
            %  xcg   --> Position along the X - Axis of the aircraft center
            %            of gravity; for typical overall geom. arrangement
            %            and loading condition, the xcg could be anywhere
            %            between 0.2 and 0.5 when Xcg is divided by the
            %            Mean Aerodynamic chord
            %  bcg   --> Dimensional, vertical arm of the lift force
            %            horizontal component; this must be divided by the
            %            Mean Aerodynamic chord before to apply the
            %            equation for the pitching mom. equilibrium
            %  MAC   --> Mean Aerodynamic chord
            %
            %  OUTPUTS
            %  DRAG PITCHING MOMENT CONTRIBUTION
            %  It's the main lifting surface contribution to pitching
            %  moment and its drag.
            %
            %  The function calculates the drag coefficient contribution to 
            %  the pitching moment contribution is:
            %
            %                              bcg
            %  CDWB_contr = CD*cos(alpha)* --- - CD*sin(alpha)*(xac - xcg)
            %                              MAC  
            
            obj  = CD*cos(alpha)*(bcg/MAC) - CD*sin(alpha)*(xac - xcg);
        end
        % -----------------------------------------------------------------
        function obj = CT_contr(obj, CD, h, MAC)
            % METHOD 3
            % Thrust contribution to pitching mom. = CTcontrib(obj, CD, h, MAC)
            % 
            %  This function is able to calculate the thrust contribution
            %  related to the pitching moment equilibrium during a
            %  particular flight condition, defined in terms of airspeed V
            %  and load factor n. The equilibrium is taken about the
            %  aircraft center of gravity. To model the thrust, we start 
            %  the obvious consideration that thrust must be at least equal
            %  to drag for sustained, levelled, steady flight; therefore,
            %  in this function the thrust coefficient is considered equal
            %  to the drag coefficient. The moment arm of thrust action 
            %  line with respect to the aircraft center of gravity is
            %  indicated with h.
            %  
            %  INPUTS 
            %  CD    --> Input drag coefficient, related to a particular
            %            flight condition which is completely defined in
            %            terms of airspeed V and load factor n (aircraft
            %            weight is a known quantity). It must be 
            %            remembered that the CD used here is relative to
            %            a wing - body configuration (wing - fuselage CD
            %            without horizontal empennage)
            %  h     --> Moment arm measured in meters of the thrust action
            %            line with respect to the center of gravity; this
            %            value must be divided by the Mean Aerodynamic
            %            chord 
            %  MAC   --> Mean Aerodynamic chord
            %
            %  OUTPUTS
            %  THRUST PITCHING MOMENT CONTRIBUTION
            %  It's the pitching moment contribution associated with thrust
            %  produced by the power train.
            %
            %  The function calculates the drag coefficient contribution to 
            %  the pitching moment contribution is:
            %
            %                  h
            %  CTcontr = CD * ---
            %                 MAC    
            
            obj = CD*(h/MAC);
        end
        % -----------------------------------------------------------------
        function obj = CM_aboutcg(obj, CM0, CM_landgear, CMCL, CL)
            % METHOD 4
            % Pitching mom. about c.g. = CM_aboutcg(obj, CM0, CM_landgear, CMCL, CL)
            % 
            %  This function is able to calculate the thrust contribution
            %  related to the pitching moment equilibrium during a
            %  particular flight condition, defined in terms of airspeed V
            %  and load factor n. The equilibrium is taken about the
            %  aircraft center of gravity. To model the thrust, we start 
            %  the obvious consideration that thrust must be at least equal
            %  to drag for sustained, levelled, steady flight; therefore,
            %  in this function the thrust coefficient is considered equal
            %  to the drag coefficient. The moment arm of thrust action 
            %  line with respect to the aircraft center of gravity is
            %  indicated with h.
            %  
            %  INPUTS 
            %  CM0         --> Zero - lift pitching moment of the wing -
            %                  body configuration 
            %  CM_landgear --> Pitching moment associated with the
            %                  extendend landing gear flight condition; it
            %                  is a constant value
            %  CMCL        --> Pitching moment curve slope, defined in
            %                  terms of lift coefficient variation
            %                                              d Cm
            %                  Pitching mom. curve slope = ----
            %                                              d CL
            %  CL          --> Input lift coefficient, related to a
            %                  particular flight condition which is 
            %                  completely defined in terms of airspeed V 
            %                  and load factor n (aircraft weight is a 
            %                  known quantity). It must be remembered that
            %                  the CL used here is relative to a 
            %                  wing - body configuration (wing - fuselage 
            %                  CL without horizontal empennage)
            %
            %  OUTPUTS
            %  PURE PITCHING MOMENT CONTRIBUTION
            %  Pure pitching moment associated with a typical flight 
            %  vehicle.
            %
            %  The function calculates the CM wing - body coefficient:
            %
            %                 d Cm
            %  CMcontr = CL * ---- + CM0 + (CM)_Landing_Gear
            %                 d CL
            
            obj = CM0 + CM_landgear + CL*(CMCL);
        end
        % -----------------------------------------------------------------
        function obj = CL_Tail(obj, CLwb_contr, CDwb_contr, CT_contr, CM, l_tail, MAC, xac, xcg, alpha)
            % METHOD 5
            % Horizontal tail lift coefficient = CL_Tail(obj, CLwb_contr, CDwb_contr, CT_contr, CM, l_tail, MAC, xac, xcg, alpha)
            % 
            %  This function is able to calculate the horizontal tail lift
            %  coefficient associated with the prescribed flight condition,
            %  defined in terms of airspeed V and load factor n. This lift
            %  coefficient is essential to correctly size the horizontal
            %  tail structure and to get a knowledge about the extra lift
            %  that the main lifting surface must produce to achieve
            %  steady, non - accelerated, levelled flight to compensate for
            %  the airloads acting on the horizontal tail. 
            %  
            %  INPUTS 
            %  CLwb_contr --> Pitching mom. contribution associated with
            %                 wing - body lift
            %  CDwb_contr --> Pitching mom. contribution associated with
            %                 wing - body drag
            %  CT_contr   --> Pitching mom. associated with non baricentral
            %                 thrust 
            %  CM         --> Pitching mom. of the wing body configuration 
            %  l_tail     --> Horizontal tail mom. arm 
            %  MAC        --> Mean Aerodynamic chord
            %  xac        --> Aerodynamic center position in unit of mean
            %                 aerodynamic chord
            %  xcg        --> Center of gravity position in unit of mean
            %                 aerodynamic chord
            %  alpha      --> Angle of attack at the prescribed flight
            %                 condition
            %
            %  OUTPUTS
            %  CL_Tail   ---> Tail lift coefficient
            %
            %  The function calculates the CL_Tail with the following
            %  equilibrium relationship:
            %
            %             (CLwb_contr + CDwb_contr + CT_contr + CM)
            %  CL_Tail =  -----------------------------------------
            %                l_tail    
            %               -------- + (xac - xcg) * cos(alpha)
            %                  MAC
            %
            
            obj = (CLwb_contr + CDwb_contr + CT_contr + CM)/((l_tail/MAC) + (xac - xcg)*cos(alpha));
        end
        % -----------------------------------------------------------------
        function obj = CL_Non_linear_model(obj, a, b, c)
            % METHOD 1
            % Maximum Lift Coefficient = CL_Non_linear_model(a, b, c)
            %   
            %   This method introduce a non-linear model for the lift
            %   coefficient calculation, able to take into account the fact
            %   that wings must compensate for the extra lift produced by the
            %   horizontal tail to achieve trimmed flight in any possible
            %   conditions. The function construct a second order
            %   polynomial to evaluate CL at high angle of attack.
            %   Applicability of this particular function is strictly
            %   limited to the case of the CS - VLA certification process
            %   for a MTOW = 100 kg drone examined as a test case.
            %   The equation used is the following:
            %   
            %   CLMax = a*(alpha)^2 + b*(alpha) + c
            %
            %   INPUTS 
            %   a     --> Coefficient of the second order term inside the used
            %             equation
            %   b     --> Coefficient of the linear term inside the used
            %             equation
            %   c     --> Constant coefficient inside the used equation 
            %   
            %   OUTPUTS 
            %   CL    --> Lift coefficients as function of angle of attack
            %   
            %   NOTE: It must be noticed that is not strictly necessary to
            %         define a wide range of AoA for the purpose of this
            %         function. This is done here for illustration
            %         purposes.
            
            AOA_aux = linspace(-20.0, 20.0, 1000);
            vect    = a.*(AOA_aux.^2) + b.*AOA_aux + c;
            obj     = vect;
        end
        % -----------------------------------------------------------------
        function obj = CL_WB_CompleteCurve(obj, CLzero, CLalpha, a, b, c, CLstar, AoA)
            % METHOD 1
            % Maximum Lift Coefficient = CL_Non_linear_model(a, b, c)
            %   
            %   This method introduce a non-linear model for the lift
            %   coefficient calculation, able to take into account the fact
            %   that wings must compensate for the extra lift produced by the
            %   horizontal tail to achieve trimmed flight in any possible
            %   conditions. The function construct a second order
            %   polynomial to evaluate CL at high angle of attack.
            %   Applicability of this particular function is strictly
            %   limited to the case of the CS - VLA certification process
            %   for a MTOW = 100 kg drone examined as a test case.
            %   The equation used is the following:
            %   
            %   CLMax = a*(alpha)^2 + b*(alpha) + c
            %
            %   INPUTS 
            %   a     --> Coefficient of the second order term inside the used
            %             equation
            %   b     --> Coefficient of the linear term inside the used
            %             equation
            %   c     --> Constant coefficient inside the used equation 
            %   
            %   OUTPUTS 
            %   CL    --> Lift coefficients as function of angle of attack
            %   
            %   NOTE: It must be noticed that is not strictly necessary to
            %         define a wide range of AoA for the purpose of this
            %         function. This is done here for illustration
            %         purposes.
            
            CL = CLzero + CLalpha*AoA;
            if length(CL) > 1
                for i = 1:length(CL)
                    x = CL(i) - CLstar; 
                    if abs(x) < 1.0*1e-2
                        j = i;
                        for z = j:length(CL)
                            CL(z) = a.*(AoA(z).^2) + b.*AoA(z) + c;
                        end
                    end
                end
            elseif length(CL) == 1
                x = CL - CLstar; 
                if abs(x) < 1.0*1e-2
                    CL = a*AoA^2 + b*AoA + c; 
                else
                    CL = CL;
                end
            end
            obj = CL;
        end
        % -----------------------------------------------------------------
    end
end

