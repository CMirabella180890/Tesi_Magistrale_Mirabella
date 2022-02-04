% ==== USEFUL FUNCTION DEFINED LOCALLY ====
% -------------------------------------------------------------------------
% CLMAX FUNCTION
CLmax_func = @(rho, S, V, WS, n) (2 / rho) * (1 / V^2) * (WS) * n;
% -------------------------------------------------------------------------
% GUST LOAD FACTOR - POSITIVE FLIGHT
nGust  = @(rho, V, a, kG, Ude, WS) 1 + (0.5 * rho * V * a * kG * Ude)/(WS); 
% -------------------------------------------------------------------------
% GUST LOAD FACTOR - INVERTED FLIGHT
nGust_inverted  = @(rho, V, a, kG, Ude, WS) 1 - (0.5 * rho * V * a * kG * Ude)/(WS); 
% -------------------------------------------------------------------------
% STALL SPEED FUNCTION
Vstall = @(WS, rho, CLmax, n) sqrt(WS * (2/rho) * (1/CLmax).*n); 
% -------------------------------------------------------------------------
%% FLIGHT ENVELOPE - FLAPS DEPLOYED 
%
%  In this script the flight envelope relative to a flaps-down
%  configuration of the aircraft is evaluated. Applicable airworthiness
%  rules applicable to this task are reported in the following lines. 
%
%% CS - VLA 345 High lift devices 
%  
%  (a) If flaps or similar high lift devices to be used for take-off,
%      approach or landing are installed, the aeroplane, with the flaps
%      fully deflected at VF, is assumed to be subjected to symmetrical
%      manoeuvres and gusts resulting in limit load factors within the
%      range determined by 
%      (1) Manoeuvring to a positive limit load factor of 2.0; and 
%      (2) Positive and negative gust of 7.62 m/s acting normale to the
%          flight path in level flight.
%
%  (b) VF must be assumed to be not less than 1.4*VS or 1.8*VSF, whichever
%      is greater, where 
%      VS --> Is the computed stalling speed with flaps retracted at the
%             design weight; and 
%     VSF --> Is the computed stalling speed with flaps fully extended at
%             the design weight.
%      However, if an automatic flap load limiting device is used, the
%      aeroplane may be designed for the critical combinations of airspeed
%      and flap position allowed by that device. 
% 
%  (c) In designing the flaps and supporting structures the following must
%      be accounted for:
%      (1) A head-on gust of 7.62 m/s (Equivalent airspeed).
%      (2) The slipstream effects specified in CS - VLA 457 (b). 
%
%  (d) In determining external loads on the aeroplane as a whole, thrust, 
%      slipstream and pitching acceleration may be assumed to be zero.
%
%  (e) The requirements of CS - VLA 457 and this paragraph may be complied
%      with separately or in combination. 

%% CS - VLA 457 Wing flaps 
%  (a) The wing flpas, their operating mechanisms and their supporting
%      structure must be designed for critical loads occurring in the
%      flaps-extended flight conditions with the flaps in any position.
%      However, if an automatic flap load limiting device is used, these
%      components may be designed for the critical combinations of airspeed
%      and flap position allowed by that device. 
%
%  (b) The effects of propeller slipstream, corresponding to take-off
%      power, must be taken into account not less than 1.4*VS, where VS is 
%      the computed stalling speed with flaps fully rectracted at the
%      design weight. For the investigation of slipstream effects, the load
%      factor may be assumed to be 1.0.

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% #########################################################################
% TAKEOFF
% #########################################################################
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% IMPLEMENTATION 
numb = 1e3;
nmax = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.Attributes.cs = " 345(a) ";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.wing_flaps_reg.Attributes.cs = " 347 ";

% CALCULATION OF THE LOAD FACTORS VECTOR 
n_flaps_vector = calcnflap(obj, nmax);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.n_flaps_vector.value = n_flaps_vector;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.n_flaps_vector.Attributes.unit = "g's";

% CALCULATION OF VS - CLEAN STALL SPEED
n1          = 1.0; 
rho0        = Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value;
rho         = Aircraft.Certification.ISA_Condition.Operative_ceiling.rho0.value;
Mass  
g           = Aircraft.Constants.g.value;
S           = Aircraft.Geometry.Wing.S.value;
b           = Aircraft.Geometry.Wing.b.value;
MGC         = S / b;
WS          = ( Mass * g ) / S;
CLmax_clean = Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value;
VS          = calcvs(obj, rho0, WS, CLmax_clean, n1);
% VS          = calcvs(obj, rho, WS, CLmax_clean, n1);

% CALCULATION OF VS1 - FLAPS DEPLOYED STALL SPEED
CLmax_takeoff = Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_takeoff.value;
VS1        = calcvs(obj, rho0, WS, CLmax_takeoff, n1);
% VS1        = calcvs(obj, rho, WS, CLmax_takeoff, n1);

% EVALUATION OF VF - FLAPS DEPLOYED AIRSPEED 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.VF.value = calcnVF(obj, VS, VS1);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.VF.Attributes.unit = "m/s";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.VF.Attributes.cs = " 345(b) ";

% STALLING SPEED VECTOR
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.VSpos_vec.value = calcvs(obj, rho0, WS, CLmax_takeoff, n_flaps_vector);
% Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.VSpos_vec.value = calcvs(obj, rho, WS, CLmax_takeoff, n_flaps_vector);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.VSpos_vec.Attributes.unit = "m/s";

% EVALUATION OF STALL AND FLAP MANOEUVRING POINT 
% POINT S 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointS.VS.value = VS1;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointS.VS.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointS.nS.value = 1.0;
nS = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointS.nS.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointS.nS.Attributes.unit = "g's"; 
% POINT A
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointA.VA.value = Vstall(WS, rho0, CLmax_takeoff, nmax);
% Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointA.VA.value = Vstall(WS, rho, CLmax_takeoff, nmax);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointA.VA.Attributes.unit = "m/s"; 
VA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointA.VA.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointA.nA.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
nA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointA.nA.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointA.nA.Attributes.unit = "g's"; 
% POINT F
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointF.VF.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.VF.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointF.VF.Attributes.unit = "m/s"; 
VF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointF.VF.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointF.nF.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
nF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointF.nF.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointF.nF.Attributes.unit = "g's"; 

% VECTOR OF STALL AIRSPEED TO PLOT THE ENVELOPE 
Reg            = Aircraft.Certification.Regulation.value;
Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;
VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointS.VS.value;
VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.VSpos_vec.value;
npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.n_flaps_vector.value;

n_from0toS     = linspace(0.0, nS, numb)';
V_from0toS     = VS*ones(numb, 1);

n_fromStoA     = linspace(nS, nA, numb)';
V_fromStoA     = Vstall(WS, rho0, CLmax_takeoff, n_fromStoA);
% V_fromStoA     = Vstall(WS, rho, CLmax_takeoff, n_fromStoA);

VS             = V_fromStoA(1);
V_from0toS     = VS*ones(numb, 1);

n_fromAtoF     = nmax*ones(numb, 1);
V_fromAtoF     = linspace(VA, VF, numb)';

n_fromFto0     = linspace(nF, 0.0, numb)';
V_fromFto0     = VF*ones(numb, 1);

flaps_envelope = figure;
hold on
grid on 
grid minor
ylim([-0.5 nmax+0.5])
xlim([0 VF+10])

plot(VSpos,      npos,       ':r', 'LineWidth', 0.2)
plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth',   1)
plot(V_fromAtoF, n_fromAtoF, '-r', 'LineWidth',   1)
plot(V_fromFto0, n_fromFto0, '-r', 'LineWidth',   1)
plot(V_from0toS, n_from0toS, '-r', 'LineWidth',   1)

plot(V_fromStoA(1), n_fromStoA(1), 'k.', 'MarkerSize', 10)
plot(V_fromAtoF(1), n_fromAtoF(1), 'k.', 'MarkerSize', 10)
plot(V_fromFto0(1), n_fromFto0(1), 'k.', 'MarkerSize', 10)

xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
% text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)
exportgraphics(flaps_envelope,'flapsenvelopediagramtakeoff.pdf','ContentType','vector')
exportgraphics(flaps_envelope,'flapsenvelopediagramtakeoff.png','ContentType','vector')

Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.Flaps_envelope.value = flaps_envelope;

% Saving figures inside correct folder
fprintf('Saving flapsenvelopediagramtakeoff.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile flapsenvelopediagramtakeoff.pdf Output
movefile flapsenvelopediagramtakeoff.png Output
% -----------------------------------------------------------------

% FLAPS DEPLOYED GUST ENVELOPE 
VF            = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointF.VF.value;
rho_operative = Aircraft.Certification.ISA_Condition.Operative_ceiling.rho0.value;
MAC           = Aircraft.Geometry.Wing.mac.value; 
CLalfa_rad    = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
g             = Aircraft.Constants.g.value;
Ude           = 7.62; % Gust magnitude

% CALCULATION OF THE MASS FACTOR
% mu_g = calcmug(obj, WS, MGC, CLalfa_rad, rho0, g);
mu_g = calcmug(obj, WS, MGC, CLalfa_rad, rho, g); 

% GUST ALLEVIATION FACTOR 
Kg   = calckg(obj, mu_g);

% CALCULATION OF THE GUST LOAD FACTOR AT V = VF 
nGUST_plus  = @(V) 1.0 + V .* ((0.5 * rho0 * CLalfa_rad * Kg * Ude) / ( WS ));
nGUST_minus = @(V) 1.0 - V .* ((0.5 * rho0 * CLalfa_rad * Kg * Ude) / ( WS ));
% nGUST_plus  = @(V) 1.0 + V .* ((0.5 * rho * CLalfa_rad * Kg * Ude) / ( WS ));
% nGUST_minus = @(V) 1.0 - V .* ((0.5 * rho * CLalfa_rad * Kg * Ude) / ( WS ));

V_gust      = linspace(0.0, VF, numb)'; 

n_gust_plus  = nGUST_plus(V_gust);
n_gust_minus = nGUST_minus(V_gust);

% STORE VALUES 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.Gust_envelope.mu_g.value                  = mu_g;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.Gust_envelope.mu_g.Attributes.unit        = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.Gust_envelope.Kg.value                    = Kg;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.Gust_envelope.Kg.Attributes.unit          = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.Gust_envelope.Vgust.value                 = V_gust;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.Gust_envelope.Vgust.Attributes.unit       = "m/s";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.Gust_envelope.ngust_plus.value            = n_gust_plus;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.Gust_envelope.ngust_plus.Attributes.unit  = "g's";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.Gust_envelope.ngust_minus.value           = n_gust_plus;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.Gust_envelope.ngust_minus.Attributes.unit = "g's";

% GUST ENVELOPE AND FLIGHT ENVELOPE SUPERPOSITION 
Reg            = Aircraft.Certification.Regulation.value;
Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;
VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointS.VS.value;
VF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointF.VF.value;
VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.VSpos_vec.value;
npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.n_flaps_vector.value;

flapsgust_envelope = figure;
hold on
grid on 
grid minor
ylim([-0.5 nmax+0.5])
xlim([0 VF+10])

plot(V_gust,     n_gust_plus,  ':k', 'LineWidth', 0.2)
plot(V_gust,     n_gust_minus, ':k', 'LineWidth', 0.2)
plot(VSpos,      npos,         ':r', 'LineWidth', 0.2)
plot(V_fromStoA, n_fromStoA,   '-r', 'LineWidth',   1)
plot(V_fromAtoF, n_fromAtoF,   '-b', 'LineWidth',   1)
plot(V_fromFto0, n_fromFto0,   '-b', 'LineWidth',   1)
plot(V_from0toS, n_from0toS,   '-b', 'LineWidth',   1)

plot(V_fromStoA(1), n_fromStoA(1), 'k.', 'MarkerSize', 10)
plot(V_fromAtoF(1), n_fromAtoF(1), 'k.', 'MarkerSize', 10)
plot(V_fromFto0(1), n_fromFto0(1), 'k.', 'MarkerSize', 10)

xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)
exportgraphics(flapsgust_envelope,'flaps_gust_envelopediagramtakeoff.pdf','ContentType','vector')
exportgraphics(flapsgust_envelope,'flaps_gust_envelopediagramtakeoff.png','ContentType','vector')

Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.Gust_envelope.diagram = flapsgust_envelope;
% Saving figures inside correct folder
fprintf('Saving flapsenvelopediagramtakeoff.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile flaps_gust_envelopediagramtakeoff.pdf Output
movefile flaps_gust_envelopediagramtakeoff.png Output        

% FINAL ENVELOPE WITH FLAPS DEPLOYED 
nmax           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
CLalfa_rad     = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.n_flaps_vector.value;
VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.VSpos_vec.value;
VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointS.VS.value;
V_gust         = linspace(VA, VF, 1e3*numb)'; 
n_g_plus       = nGUST_plus(V_gust);
n_g_minus      = nGUST_minus(V_gust);
nS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointS.nS.value;
VF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointF.VF.value;
nF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointF.nF.value;
Reg            = Aircraft.Certification.Regulation.value;
Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;

syms a b c V 
a        = (rho0 * CLmax_takeoff) / (2 * WS);
b        = (Kg * Ude * CLalfa_rad * rho0)/(2 * WS);
% a        = (rho * CLmax_takeoff) / (2 * WS);
% b        = (Kg * Ude * CLalfa_rad * rho)/(2 * WS);
c        = 1;
eqn      = a*V^2 - b*V - c;
Solution = vpasolve(eqn, V);

for i = 1:length(Solution)
    if Solution(i) > 0
        new_VA = cast(Solution(i), 'double');
        if abs(new_VA) > VA
            VA = abs(new_VA);
            nA = nGUST_plus(VA);

            % CASE STUDY FOR FLAPS FINAL ENVELOPE
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.case.value = 'Case 1';

            n_from0toS     = linspace(0.0, nS, numb)';
            V_from0toS     = VS * ones(numb, 1);
            nS             = n_from0toS(end);
            VS             = V_from0toS(end);

            n_fromStoA     = linspace(nS, nA, numb)';
            V_fromStoA     = Vstall(WS, rho0, CLmax_takeoff, n_fromStoA);
%             V_fromStoA     = Vstall(WS, rho, CLmax_takeoff, n_fromStoA);
            nA             = n_fromStoA(end);
            VA             = V_fromStoA(end);
            
            % DEFINITION OF NEW V FROM 0 TO S
            VS         = V_fromStoA(1);
            V_from0toS = VS * ones(numb, 1);

            V_fromAtoF     = linspace(VA, VF, numb)';
            n_fromAtoF     = linspace(nA, nmax, numb)';  
            nF             = n_fromAtoF(end);
            VF             = V_fromAtoF(end);      

            n_fromFto0     = linspace(nF, 0.0, numb)';
            V_fromFto0     = VF * ones(numb, 1);   
            % --------------------------------------------------------------------------------------------------------------------    
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointS.nS.value            = n_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointS.nS.Attributes.unit  = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointS.VS.value            = V_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointS.VS.Attributes.unit  = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_from0toS.value           = n_from0toS;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_from0toS.Attributes.unit = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_from0toS.value           = V_from0toS;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_from0toS.Attributes.unit = "m/s"; 
            % --------------------------------------------------------------------------------------------------------------------    
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointA.nA.value            = n_fromStoA(end);
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointA.nA.Attributes.unit  = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointA.VA.value            = V_fromStoA(end);
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointA.VA.Attributes.unit  = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromStoA.value           = n_fromStoA;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromStoA.Attributes.unit = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromStoA.value           = V_fromStoA;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromStoA.Attributes.unit = "m/s"; 
            % --------------------------------------------------------------------------------------------------------------------
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointF.nF.value            = n_fromAtoF(end);
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointF.nF.Attributes.unit  = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointF.VF.value            = V_fromAtoF(end);
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointF.VF.Attributes.unit  = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromAtoF.value           = n_fromAtoF;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromAtoF.Attributes.unit = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromAtoF.value           = V_fromAtoF;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromAtoF.Attributes.unit = "m/s";
            % --------------------------------------------------------------------------------------------------------------------
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromFto0.value           = n_fromFto0;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromFto0.Attributes.unit = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromFto0.value           = V_fromFto0;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromFto0.Attributes.unit = "m/s"; 
            
            final_envelope = figure;
            hold on
            grid on 
            grid minor
            ylim([-0.5 nmax+0.5])
            xlim([0 VF+10])

            plot(V_from0toS,  n_from0toS,   '-r', 'LineWidth',   1)
            plot(V_fromStoA,  n_fromStoA,   '-r', 'LineWidth',   1)
            plot(V_fromAtoF,  n_fromAtoF,   '-b', 'LineWidth',   1)
            plot(V_fromFto0,  n_fromFto0,   '-b', 'LineWidth',   1)

            plot(V_fromStoA(1),  n_fromStoA(1),  'k.', 'MarkerSize', 10)
            plot(V_fromAtoF(1),  n_fromAtoF(1),  'k.', 'MarkerSize', 10)
            plot(V_fromFto0(1),  n_fromFto0(1),  'k.', 'MarkerSize', 10)

            xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
            ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
            title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
            text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
            text(VS, nS,   '\fontname{Courier} S', 'FontSize', 6)
            text(VA, nA,   '\fontname{Courier} A', 'FontSize', 6)
            text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)   

            exportgraphics(final_envelope,'flaps_final_envelopediagramtakeoff.pdf','ContentType','vector')
            exportgraphics(final_envelope,'flaps_final_envelopediagramtakeoff.png','ContentType','vector')

            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.diagram = final_envelope;

            % Saving figures inside correct folder
            fprintf('Saving flaps_final_envelopediagramtakeoff.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile flaps_final_envelopediagramtakeoff.pdf Output
            movefile flaps_final_envelopediagramtakeoff.png Output   
            
        end


        if abs(new_VA) < VA

            % CASE STUDY FOR FLAPS FINAL ENVELOPE
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.case.value = 'Case 2';

            VA = VA; 
            nA = nA;

            if abs(n_g_plus) > nmax 

                n_from0toS     = linspace(0.0, nS, numb)';
                V_from0toS     = VS * ones(numb, 1);
                nS             = n_from0toS(end);
                VS             = V_from0toS(end);

                n_fromStoA     = linspace(nS, nA, numb)';
                V_fromStoA     = Vstall(WS, rho0, CLmax_takeoff, n_fromStoA);
%                 V_fromStoA     = Vstall(WS, rho, CLmax_takeoff, n_fromStoA);
                nA             = n_fromStoA(end);
                VA             = V_fromStoA(end);

                tol          = 1e-3;
                for i = 1:length(V_gust)
                    x = n_g_plus(i);
                    y = x - nmax;
                    if abs(y) < tol
                        row = i; 
                        VF1 = V_gust(row);
                    end
                end
                V_fromAtoF1 = linspace(VA, VF1, numb)';
                n_fromAtoF1 = nGUST_plus(V_fromAtoF1);
                nF1         = nGUST_plus(VF1); 

                V_fromF1toF = linspace(VF1, VF, numb)';
                n_fromF1toF = linspace(nF1, nF, numb)';

                V_fromFto0  = VF * ones(numb, 1);
                n_fromFto0  = linspace(nF, 0.0, numb)';
                % --------------------------------------------------------------------------------------------------------------------    
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointS.nS.value             = n_from0toS(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointS.nS.Attributes.unit   = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointS.VS.value             = V_from0toS(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointS.VS.Attributes.unit   = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_from0toS.value            = n_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_from0toS.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_from0toS.value            = V_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_from0toS.Attributes.unit  = "m/s"; 
                % -------------------------------------------------------------------------------------------------------------------- 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointA.nA.value             = n_fromStoA(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointA.nA.Attributes.unit   = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointA.VA.value             = V_fromStoA(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointA.VA.Attributes.unit   = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromStoA.value            = n_fromStoA;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromStoA.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromStoA.value            = V_fromStoA;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromStoA.Attributes.unit  = "m/s"; 
                % --------------------------------------------------------------------------------------------------------------------
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointF1.nF1.value           = n_fromAtoF1(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointF1.nF1.Attributes.unit = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointF1.VF1.value           = V_fromAtoF1(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointF1.VF1.Attributes.unit = "m/s";  
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromAtoF1.value           = n_fromAtoF1;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromAtoF1.Attributes.unit = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromAtoF1.value           = V_fromAtoF1;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromAtoF1.Attributes.unit = "m/s"; 
                % -------------------------------------------------------------------------------------------------------------------- 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointF.nF.value             = n_fromF1toF(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointF.nF.Attributes.unit   = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointF.VF.value             = V_fromF1toF(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointF.VF.Attributes.unit   = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromF1toF.value           = n_fromF1toF;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromF1toF.Attributes.unit = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromF1toF.value           = V_fromF1toF;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromF1toF.Attributes.unit = "m/s"; 
                % -------------------------------------------------------------------------------------------------------------------- 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromFto0.value            = n_fromFto0;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromFto0.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromFto0.value            = V_fromFto0;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromFto0.Attributes.unit  = "m/s"; 
                
                final_envelope = figure;
                hold on
                grid on 
                grid minor
                ylim([-0.5 nmax+0.5])
                xlim([0 VF+10])

                plot(V_from0toS,  n_from0toS,   '-r', 'LineWidth',   1)
                plot(V_fromStoA,  n_fromStoA,   '-r', 'LineWidth',   1)
                plot(V_fromAtoF1, n_fromAtoF1,  '-b', 'LineWidth',   1)
                plot(V_fromF1toF, n_fromF1toF,  '-b', 'LineWidth',   1)
                plot(V_fromFto0,  n_fromFto0,   '-b', 'LineWidth',   1)

                plot(V_fromStoA(1),  n_fromStoA(1),  'k.', 'MarkerSize', 10)
                plot(V_fromAtoF1(1), n_fromAtoF1(1), 'k.', 'MarkerSize', 10)
                plot(V_fromF1toF(1), n_fromF1toF(1), 'k.', 'MarkerSize', 10)
                plot(V_fromFto0(1),  n_fromFto0(1),  'k.', 'MarkerSize', 10)

                xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
                ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
                title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
                text(15, 1.8, Aircraft_name)                                      % Aircraft name inside the plot
                text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
                text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)   

                exportgraphics(final_envelope,'flaps_final_envelopediagramtakeoff.pdf','ContentType','vector')
                exportgraphics(final_envelope,'flaps_final_envelopediagramtakeoff.png','ContentType','vector')

                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.diagram = final_envelope;

                % Saving figures inside correct folder
                fprintf('Saving flaps_final_envelopediagram.pdf in: ');
                fprintf('\n'); 
                fprintf('%s\n', SaveFolder);
                % Moving file inside correct folder
                movefile flaps_final_envelopediagramtakeoff.pdf Output
                movefile flaps_final_envelopediagramtakeoff.png Output               

            elseif abs(n_g_plus) < nmax

                n_from0toS     = linspace(0.0, nS, numb)';
                V_from0toS     = VS * ones(numb, 1);
                nS             = n_from0toS(end);
                VS             = V_from0toS(end);

                n_fromStoA     = linspace(nS, nA, numb)';
                V_fromStoA     = Vstall(WS, rho0, CLmax_takeoff, n_fromStoA);
%                 V_fromStoA     = Vstall(WS, rho, CLmax_takeoff, n_fromStoA);
                nA             = n_fromStoA(end);
                VA             = V_fromStoA(end);

                V_fromAtoF = linspace(VA, VF, numb)';
                n_fromAtoF = nmax * ones(numb, 1);

                V_fromFto0  = VF * ones(numb, 1);
                n_fromFto0  = linspace(nF, 0.0, numb)';
                
                % --------------------------------------------------------------------------------------------------------------------    
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointS.nS.value             = n_from0toS(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointS.nS.Attributes.unit   = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointS.VS.value             = V_from0toS(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointS.VS.Attributes.unit   = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_from0toS.value            = n_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_from0toS.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_from0toS.value            = V_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_from0toS.Attributes.unit  = "m/s"; 
                % -------------------------------------------------------------------------------------------------------------------- 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointA.nA.value             = n_fromStoA(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointA.nA.Attributes.unit   = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointA.VA.value             = V_fromStoA(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointA.VA.Attributes.unit   = "m/s";    
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromStoA.value            = n_fromStoA;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromStoA.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromStoA.value            = V_fromStoA;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromStoA.Attributes.unit  = "m/s"; 
                % -------------------------------------------------------------------------------------------------------------------- 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointF.nF.value             = n_fromAtoF(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointF.nF.Attributes.unit   = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointF.VF.value             = V_fromAtoF(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.PointF.VF.Attributes.unit   = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromAtoF.value            = n_fromAtoF;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromAtoF.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromAtoF.value            = V_fromAtoF;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromAtoF.Attributes.unit  = "m/s"; 
                % -------------------------------------------------------------------------------------------------------------------- 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromFto0.value            = n_fromFto0;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromFto0.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromFto0.value            = V_fromFto0;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromFto0.Attributes.unit  = "m/s"; 

                final_envelope = figure;
                hold on
                grid on 
                grid minor
                ylim([-0.5 nmax+0.5])
                xlim([0 VF+10])

                plot(V_from0toS,  n_from0toS,   '-r', 'LineWidth',   1)
                plot(V_fromStoA,  n_fromStoA,   '-r', 'LineWidth',   1)
                plot(V_fromAtoF,  n_fromAtoF,   '-b', 'LineWidth',   1)
                plot(V_fromFto0,  n_fromFto0,   '-b', 'LineWidth',   1)

                plot(V_fromStoA(1),  n_fromStoA(1),  'k.', 'MarkerSize', 10)
                plot(V_fromAtoF(1),  n_fromAtoF(1), 'k.', 'MarkerSize', 10)
                plot(V_fromFto0(1),  n_fromFto0(1),  'k.', 'MarkerSize', 10)

                xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
                ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
                title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
                text(15, 1.8, Aircraft_name)                                      % Aircraft name inside the plot
                text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
                text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)   

                exportgraphics(final_envelope,'flaps_final_envelopediagramtakeoff.pdf','ContentType','vector')
                exportgraphics(final_envelope,'flaps_final_envelopediagramtakeoff.png','ContentType','vector')

                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.diagram = final_envelope;

                % Saving figures inside correct folder
                fprintf('Saving flaps_final_envelopediagramtakeoff.pdf in: ');
                fprintf('\n'); 
                fprintf('%s\n', SaveFolder);
                % Moving file inside correct folder
                movefile flaps_final_envelopediagramtakeoff.pdf Output
                movefile flaps_final_envelopediagramtakeoff.png Output               

            end
        end
    end
end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% #########################################################################
% #########################################################################
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% #########################################################################
% LANDING 
% #########################################################################
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% IMPLEMENTATION 
numb = 1e3;
nmax = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;

% CALCULATION OF THE LOAD FACTORS VECTOR 
n_flaps_vector = calcnflap(obj, nmax);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.n_flaps_vector.value = n_flaps_vector;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.n_flaps_vector.Attributes.unit = "g's";

% CALCULATION OF VS - CLEAN STALL SPEED
n1          = 1.0; 
rho0        = Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value;
rho         = Aircraft.Certification.ISA_Condition.Operative_ceiling.rho0.value;
Mass  
g           = Aircraft.Constants.g.value;
S           = Aircraft.Geometry.Wing.S.value;
b           = Aircraft.Geometry.Wing.b.value;
MGC         = S / b;
WS          = ( Mass * g ) / S;
CLmax_clean = Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value;
VS          = calcvs(obj, rho0, WS, CLmax_clean, n1);
% VS          = calcvs(obj, rho, WS, CLmax_clean, n1);

% CALCULATION OF VS1 - FLAPS DEPLOYED STALL SPEED
CLmax_landing = Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_landing.value;
VS1        = calcvs(obj, rho0, WS, CLmax_landing, n1);
% VS1        = calcvs(obj, rho, WS, CLmax_landing, n1);

% EVALUATION OF VF - FLAPS DEPLOYED AIRSPEED 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.VF.value = calcnVF(obj, VS, VS1);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.VF.Attributes.unit = "m/s";

% STALLING SPEED VECTOR
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.VSpos_vec.value = calcvs(obj, rho0, WS, CLmax_landing, n_flaps_vector);
% Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.VSpos_vec.value = calcvs(obj, rho, WS, CLmax_landing, n_flaps_vector);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.VSpos_vec.Attributes.unit = "m/s";

% EVALUATION OF STALL AND FLAP MANOEUVRING POINT 
% POINT S 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointS.VS.value = VS1;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointS.VS.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointS.nS.value = 1.0;
nS = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointS.nS.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointS.nS.Attributes.unit = "g's"; 
% POINT A
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointA.VA.value = Vstall(WS, rho0, CLmax_landing, nmax);
% Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointA.VA.value = Vstall(WS, rho, CLmax_landing, nmax);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointA.VA.Attributes.unit = "m/s"; 
VA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointA.VA.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointA.nA.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
nA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointA.nA.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointA.nA.Attributes.unit = "g's"; 
% POINT F
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointF.VF.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.VF.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointF.VF.Attributes.unit = "m/s"; 
VF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointF.VF.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointF.nF.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
nF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointF.nF.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointF.nF.Attributes.unit = "g's"; 

% VECTOR OF STALL AIRSPEED TO PLOT THE ENVELOPE 
Reg            = Aircraft.Certification.Regulation.value;
Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;
VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointS.VS.value;
VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.VSpos_vec.value;
npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.n_flaps_vector.value;

n_from0toS     = linspace(0.0, nS, numb)';
V_from0toS     = VS*ones(numb, 1);

n_fromStoA     = linspace(nS, nA, numb)';
V_fromStoA     = Vstall(WS, rho0, CLmax_landing, n_fromStoA);
% V_fromStoA     = Vstall(WS, rho, CLmax_landing, n_fromStoA);

n_fromAtoF     = nmax*ones(numb, 1);
V_fromAtoF     = linspace(VA, VF, numb)';

n_fromFto0     = linspace(nF, 0.0, numb)';
V_fromFto0     = VF*ones(numb, 1);

flaps_envelope = figure;
hold on
grid on 
grid minor
ylim([-0.5 nmax+0.5])
xlim([0 VF+10])

plot(VSpos,      npos,       ':r', 'LineWidth', 0.2)
plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth',   1)
plot(V_fromAtoF, n_fromAtoF, '-r', 'LineWidth',   1)
plot(V_fromFto0, n_fromFto0, '-r', 'LineWidth',   1)
plot(V_from0toS, n_from0toS, '-r', 'LineWidth',   1)

plot(V_fromStoA(1), n_fromStoA(1), 'k.', 'MarkerSize', 10)
plot(V_fromAtoF(1), n_fromAtoF(1), 'k.', 'MarkerSize', 10)
plot(V_fromFto0(1), n_fromFto0(1), 'k.', 'MarkerSize', 10)

xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
% text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)
exportgraphics(flaps_envelope,'flapsenvelopediagramlanding.pdf','ContentType','vector')
exportgraphics(flaps_envelope,'flapsenvelopediagramlanding.png','ContentType','vector')

Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.Flaps_envelope.value = flaps_envelope;

% Saving figures inside correct folder
fprintf('Saving flapsenvelopediagramlanding.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile flapsenvelopediagramlanding.pdf Output
movefile flapsenvelopediagramlanding.png Output
% -----------------------------------------------------------------

% FLAPS DEPLOYED GUST ENVELOPE 
VF            = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointF.VF.value;
rho0          = Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value;
rho_operative = Aircraft.Certification.ISA_Condition.Operative_ceiling.rho0.value;
MAC           = Aircraft.Geometry.Wing.mac.value; 
CLalfa_rad    = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
g             = Aircraft.Constants.g.value;
Ude           = 7.62; % Gust magnitude

% CALCULATION OF THE MASS FACTOR
% mu_g = calcmug(obj, WS, MGC, CLalfa_rad, rho0, g); 
mu_g = calcmug(obj, WS, MGC, CLalfa_rad, rho_operative, g); 

% GUST ALLEVIATION FACTOR 
Kg   = calckg(obj, mu_g);

% CALCULATION OF THE GUST LOAD FACTOR AT V = VF 
nGUST_plus  = @(V) 1.0 + V .* ((0.5 * rho0 * CLalfa_rad * Kg * Ude) / ( WS ));
nGUST_minus = @(V) 1.0 - V .* ((0.5 * rho0 * CLalfa_rad * Kg * Ude) / ( WS ));
% nGUST_plus  = @(V) 1.0 + V .* ((0.5 * rho_operative * CLalfa_rad * Kg * Ude) / ( WS ));
% nGUST_minus = @(V) 1.0 - V .* ((0.5 * rho_operative * CLalfa_rad * Kg * Ude) / ( WS ));

V_gust      = linspace(0.0, VF, numb)'; 

n_gust_plus  = nGUST_plus(V_gust);
n_gust_minus = nGUST_minus(V_gust);

% STORE VALUES 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.Gust_envelope.mu_g.value                  = mu_g;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.Gust_envelope.mu_g.Attributes.unit        = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.Gust_envelope.Kg.value                    = Kg;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.Gust_envelope.Kg.Attributes.unit          = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.Gust_envelope.Vgust.value                 = V_gust;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.Gust_envelope.Vgust.Attributes.unit       = "m/s";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.Gust_envelope.ngust_plus.value            = n_gust_plus;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.Gust_envelope.ngust_plus.Attributes.unit  = "g's";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.Gust_envelope.ngust_minus.value           = n_gust_plus;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.Gust_envelope.ngust_minus.Attributes.unit = "g's";

% GUST ENVELOPE AND FLIGHT ENVELOPE SUPERPOSITION 
Reg            = Aircraft.Certification.Regulation.value;
Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;
VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointS.VS.value;
VF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointF.VF.value;
VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.VSpos_vec.value;
npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.n_flaps_vector.value;

flapsgust_envelope = figure;
hold on
grid on 
grid minor
ylim([-0.5 nmax+0.5])
xlim([0 VF+10])

plot(V_gust,     n_gust_plus,  ':k', 'LineWidth', 0.2)
plot(V_gust,     n_gust_minus, ':k', 'LineWidth', 0.2)
plot(VSpos,      npos,         ':r', 'LineWidth', 0.2)
plot(V_fromStoA, n_fromStoA,   '-r', 'LineWidth',   1)
plot(V_fromAtoF, n_fromAtoF,   '-b', 'LineWidth',   1)
plot(V_fromFto0, n_fromFto0,   '-b', 'LineWidth',   1)
plot(V_from0toS, n_from0toS,   '-b', 'LineWidth',   1)

plot(V_fromStoA(1), n_fromStoA(1), 'k.', 'MarkerSize', 10)
plot(V_fromAtoF(1), n_fromAtoF(1), 'k.', 'MarkerSize', 10)
plot(V_fromFto0(1), n_fromFto0(1), 'k.', 'MarkerSize', 10)

xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)
exportgraphics(flapsgust_envelope,'flaps_gust_envelopediagramlanding.pdf','ContentType','vector')
exportgraphics(flapsgust_envelope,'flaps_gust_envelopediagramlanding.png','ContentType','vector')

Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.Gust_envelope.diagram = flapsgust_envelope;
% Saving figures inside correct folder
fprintf('Saving flapsenvelopediagramlanding.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile flaps_gust_envelopediagramlanding.pdf Output
movefile flaps_gust_envelopediagramlanding.png Output        

% FINAL ENVELOPE WITH FLAPS DEPLOYED 
nmax           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
CLalfa_rad     = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.n_flaps_vector.value;
VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.VSpos_vec.value;
VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointS.VS.value;
V_gust         = linspace(VA, VF, 1e3*numb)'; 
n_g_plus       = nGUST_plus(V_gust);
n_g_minus      = nGUST_minus(V_gust);
nS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointS.nS.value;
VF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointF.VF.value;
nF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointF.nF.value;
Reg            = Aircraft.Certification.Regulation.value;
Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;

syms a b c V 
a        = (rho0 * CLmax_landing) / (2 * WS);
b        = (Kg * Ude * CLalfa_rad * rho0)/(2 * WS);
% a        = (rho_operative * CLmax_landing) / (2 * WS);
% b        = (Kg * Ude * CLalfa_rad * rho_operative)/(2 * WS);
c        = 1;
eqn      = a*V^2 - b*V - c;
Solution = vpasolve(eqn, V);

for i = 1:length(Solution)
    if Solution(i) > 0
        new_VA = cast(Solution(i), 'double');
        if abs(new_VA) > VA
            VA = abs(new_VA);
            nA = nGUST_plus(VA);

            % CASE STUDY FOR FLAPS FINAL ENVELOPE
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.case.value = 'Case 1';

            n_from0toS     = linspace(0.0, nS, numb)';
            V_from0toS     = VS * ones(numb, 1);
            nS             = n_from0toS(end);
            VS             = V_from0toS(end);

            n_fromStoA     = linspace(nS, nA, numb)';
            V_fromStoA     = Vstall(WS, rho0, CLmax_landing, n_fromStoA);
%             V_fromStoA     = Vstall(WS, rho_operative, CLmax_landing, n_fromStoA);
            nA             = n_fromStoA(end);
            VA             = V_fromStoA(end);

            V_fromAtoF     = linspace(VA, VF, numb)';
            n_fromAtoF     = linspace(nA, nmax, numb)';  
            nF             = n_fromAtoF(end);
            VF             = V_fromAtoF(end);      

            n_fromFto0     = linspace(nF, 0.0, numb)';
            V_fromFto0     = VF * ones(numb, 1);        

            % --------------------------------------------------------------------------------------------------------------------    
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointS.nS.value             = n_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointS.nS.Attributes.unit   = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointS.VS.value             = V_from0toS(end);
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointS.VS.Attributes.unit   = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_from0toS.value           = n_from0toS;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_from0toS.Attributes.unit = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_from0toS.value           = V_from0toS;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_from0toS.Attributes.unit = "m/s"; 
            % -------------------------------------------------------------------------------------------------------------------- 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointA.nA.value             = n_fromStoA(end);
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointA.nA.Attributes.unit   = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointA.VA.value             = V_fromStoA(end);
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointA.VA.Attributes.unit   = "m/s";  
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromStoA.value           = n_fromStoA;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromStoA.Attributes.unit = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromStoA.value           = V_fromStoA;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromStoA.Attributes.unit = "m/s"; 
            % -------------------------------------------------------------------------------------------------------------------- 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointF.nF.value             = n_fromAtoF(end);
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointF.nF.Attributes.unit   = "g's";
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointF.VF.value             = V_fromAtoF(end);
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointF.VF.Attributes.unit   = "m/s";    
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromAtoF.value           = n_fromAtoF;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromAtoF.Attributes.unit = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromAtoF.value           = V_fromAtoF;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromAtoF.Attributes.unit = "m/s";
            % --------------------------------------------------------------------------------------------------------------------   
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromFto0.value           = n_fromFto0;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromFto0.Attributes.unit = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromFto0.value           = V_fromFto0;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromFto0.Attributes.unit = "m/s"; 
            
            final_envelope = figure;
            hold on
            grid on 
            grid minor
            ylim([-0.5 nmax+0.5])
            xlim([0 VF+10])

            plot(V_from0toS,  n_from0toS,   '-r', 'LineWidth',   1)
            plot(V_fromStoA,  n_fromStoA,   '-r', 'LineWidth',   1)
            plot(V_fromAtoF,  n_fromAtoF,   '-b', 'LineWidth',   1)
            plot(V_fromFto0,  n_fromFto0,   '-b', 'LineWidth',   1)

            plot(V_fromStoA(1),  n_fromStoA(1),  'k.', 'MarkerSize', 10)
            plot(V_fromAtoF(1),  n_fromAtoF(1),  'k.', 'MarkerSize', 10)
            plot(V_fromFto0(1),  n_fromFto0(1),  'k.', 'MarkerSize', 10)

            xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
            ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
            title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
            text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
            text(VS, nS,   '\fontname{Courier} S', 'FontSize', 6)
            text(VA, nA,   '\fontname{Courier} A', 'FontSize', 6)
            text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)   

            exportgraphics(final_envelope,'flaps_final_envelopediagramlanding.pdf','ContentType','vector')
            exportgraphics(final_envelope,'flaps_final_envelopediagramlanding.png','ContentType','vector')

            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.diagram = final_envelope;

            % Saving figures inside correct folder
            fprintf('Saving flaps_final_envelopediagramlanding.pdf in: ');
            fprintf('\n'); 
            fprintf('%s\n', SaveFolder);
            % Moving file inside correct folder
            movefile flaps_final_envelopediagramlanding.pdf Output
            movefile flaps_final_envelopediagramlanding.png Output     
        end


        if abs(new_VA) < VA

            % CASE STUDY FOR FLAPS FINAL ENVELOPE
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.case.value = 'Case 2';

            VA = VA; 
            nA = nA;

            if abs(n_g_plus) > nmax 

                n_from0toS     = linspace(0.0, nS, numb)';
                V_from0toS     = VS * ones(numb, 1);
                nS             = n_from0toS(end);
                VS             = V_from0toS(end);

                n_fromStoA     = linspace(nS, nA, numb)';
                V_fromStoA     = Vstall(WS, rho0, CLmax_landing, n_fromStoA);
%                 V_fromStoA     = Vstall(WS, rho_operative, CLmax_landing, n_fromStoA);
                nA             = n_fromStoA(end);
                VA             = V_fromStoA(end);

                tol          = 1e-3;
                for i = 1:length(V_gust)
                    x = n_g_plus(i);
                    y = x - nmax;
                    if abs(y) < tol
                        row = i; 
                        VF1 = V_gust(row);
                    end
                end
                V_fromAtoF1 = linspace(VA, VF1, numb)';
                n_fromAtoF1 = nGUST_plus(V_fromAtoF1);
                nF1         = nGUST_plus(VF1); 

                V_fromF1toF = linspace(VF1, VF, numb)';
                n_fromF1toF = linspace(nF1, nF, numb)';

                V_fromFto0  = VF * ones(numb, 1);
                n_fromFto0  = linspace(nF, 0.0, numb)';

                % --------------------------------------------------------------------------------------------------------------------    
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointS.nS.value             = n_from0toS(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointS.nS.Attributes.unit   = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointS.VS.value             = V_from0toS(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointS.VS.Attributes.unit   = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_from0toS.value            = n_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_from0toS.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_from0toS.value            = V_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_from0toS.Attributes.unit  = "m/s"; 
                % -------------------------------------------------------------------------------------------------------------------- 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointA.nA.value             = n_fromStoA(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointA.nA.Attributes.unit   = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointA.VA.value             = V_fromStoA(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointA.VA.Attributes.unit   = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromStoA.value            = n_fromStoA;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromStoA.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromStoA.value            = V_fromStoA;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromStoA.Attributes.unit  = "m/s"; 
                % -------------------------------------------------------------------------------------------------------------------- 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointF1.nF1.value             = n_fromAtoF1(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointF1.nF1.Attributes.unit   = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointF1.VF1.value             = V_fromAtoF1(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointF1.VF1.Attributes.unit   = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromAtoF1.value           = n_fromAtoF1;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromAtoF1.Attributes.unit = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromAtoF1.value           = V_fromAtoF1;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromAtoF1.Attributes.unit = "m/s"; 
                % -------------------------------------------------------------------------------------------------------------------- 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointF.nF.value             = n_fromF1toF(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointF.nF.Attributes.unit   = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointF.VF.value             = V_fromF1toF(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointF.VF.Attributes.unit   = "m/s";   
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromF1toF.value           = n_fromF1toF;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromF1toF.Attributes.unit = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromF1toF.value           = V_fromF1toF;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromF1toF.Attributes.unit = "m/s"; 
                % -------------------------------------------------------------------------------------------------------------------- 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromFto0.value            = n_fromFto0;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromFto0.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromFto0.value            = V_fromFto0;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromFto0.Attributes.unit  = "m/s"; 
                
                final_envelope = figure;
                hold on
                grid on 
                grid minor
                ylim([-0.5 nmax+0.5])
                xlim([0 VF+10])

                plot(V_from0toS,  n_from0toS,   '-r', 'LineWidth',   1)
                plot(V_fromStoA,  n_fromStoA,   '-r', 'LineWidth',   1)
                plot(V_fromAtoF1, n_fromAtoF1,  '-b', 'LineWidth',   1)
                plot(V_fromF1toF, n_fromF1toF,  '-b', 'LineWidth',   1)
                plot(V_fromFto0,  n_fromFto0,   '-b', 'LineWidth',   1)

                plot(V_fromStoA(1),  n_fromStoA(1),  'k.', 'MarkerSize', 10)
                plot(V_fromAtoF1(1), n_fromAtoF1(1), 'k.', 'MarkerSize', 10)
                plot(V_fromF1toF(1), n_fromF1toF(1), 'k.', 'MarkerSize', 10)
                plot(V_fromFto0(1),  n_fromFto0(1),  'k.', 'MarkerSize', 10)

                xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
                ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
                title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
                text(15, 1.8, Aircraft_name)                                      % Aircraft name inside the plot
                text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
                text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)   

                exportgraphics(final_envelope,'flaps_final_envelopediagramlanding.pdf','ContentType','vector')
                exportgraphics(final_envelope,'flaps_final_envelopediagramlanding.png','ContentType','vector')

                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.diagram = final_envelope;

                % Saving figures inside correct folder
                fprintf('Saving flaps_final_envelopediagramlanding.pdf in: ');
                fprintf('\n'); 
                fprintf('%s\n', SaveFolder);
                % Moving file inside correct folder
                movefile flaps_final_envelopediagramlanding.pdf Output
                movefile flaps_final_envelopediagramlanding.png Output               

            elseif abs(n_g_plus) < nmax

                n_from0toS     = linspace(0.0, nS, numb)';
                V_from0toS     = VS * ones(numb, 1);
                nS             = n_from0toS(end);
                VS             = V_from0toS(end);

                n_fromStoA     = linspace(nS, nA, numb)';
                V_fromStoA     = Vstall(WS, rho0, CLmax_landing, n_fromStoA);
%                 V_fromStoA     = Vstall(WS, rho_operative, CLmax_landing, n_fromStoA);
                nA             = n_fromStoA(end);
                VA             = V_fromStoA(end);

                V_fromAtoF = linspace(VA, VF, numb)';
                n_fromAtoF = nmax * ones(numb, 1);

                V_fromFto0  = VF * ones(numb, 1);
                n_fromFto0  = linspace(nF, 0.0, numb)';

                % --------------------------------------------------------------------------------------------------------------------    
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointS.nS.value             = n_from0toS(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointS.nS.Attributes.unit   = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointS.VS.value             = V_from0toS(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointS.VS.Attributes.unit   = "m/s";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_from0toS.value            = n_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_from0toS.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_from0toS.value            = V_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_from0toS.Attributes.unit  = "m/s"; 
                % --------------------------------------------------------------------------------------------------------------------  
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointA.nA.value             = n_fromStoA(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointA.nA.Attributes.unit   = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointA.VA.value             = V_fromStoA(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointA.VA.Attributes.unit   = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromStoA.value            = n_fromStoA;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromStoA.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromStoA.value            = V_fromStoA;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromStoA.Attributes.unit  = "m/s"; 
                % --------------------------------------------------------------------------------------------------------------------   
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointF.nF.value             = n_fromF1toF(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointF.nF.Attributes.unit   = "g's";
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointF.VF.value             = V_fromF1toF(end);
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.PointF.VF.Attributes.unit   = "m/s";   
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromAtoF.value            = n_fromAtoF;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromAtoF.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromAtoF.value            = V_fromAtoF;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromAtoF.Attributes.unit  = "m/s"; 
                % --------------------------------------------------------------------------------------------------------------------   
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromFto0.value            = n_fromFto0;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromFto0.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromFto0.value            = V_fromFto0;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromFto0.Attributes.unit  = "m/s"; 
                
                final_envelope = figure;
                hold on
                grid on 
                grid minor
                ylim([-0.5 nmax+0.5])
                xlim([0 VF+10])

                plot(V_from0toS,  n_from0toS,   '-r', 'LineWidth',   1)
                plot(V_fromStoA,  n_fromStoA,   '-r', 'LineWidth',   1)
                plot(V_fromAtoF,  n_fromAtoF,   '-b', 'LineWidth',   1)
                plot(V_fromFto0,  n_fromFto0,   '-b', 'LineWidth',   1)

                plot(V_fromStoA(1),  n_fromStoA(1),  'k.', 'MarkerSize', 10)
                plot(V_fromAtoF(1),  n_fromAtoF(1), 'k.', 'MarkerSize', 10)
                plot(V_fromFto0(1),  n_fromFto0(1),  'k.', 'MarkerSize', 10)

                xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
                ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
                title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
                text(15, 1.8, Aircraft_name)                                      % Aircraft name inside the plot
                text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
                text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)   

                exportgraphics(final_envelope,'flaps_final_envelopediagramlanding.pdf','ContentType','vector')
                exportgraphics(final_envelope,'flaps_final_envelopediagramlanding.png','ContentType','vector')

                Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.diagram = final_envelope;

                % Saving figures inside correct folder
                fprintf('Saving flaps_final_envelopediagramlanding.pdf in: ');
                fprintf('\n'); 
                fprintf('%s\n', SaveFolder);
                % Moving file inside correct folder
                movefile flaps_final_envelopediagramlanding.pdf Output
                movefile flaps_final_envelopediagramlanding.png Output               

            end
        end
    end
end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% #########################################################################
% #########################################################################
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% #########################################################################
% SUPERPOSITION
% #########################################################################
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% FLAPS - TAKEOFF
switch (Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.case.value)
    case 'Case 1'
        
        V_flap_takeoff_from0toS = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_from0toS.value;
        n_flap_takeoff_from0toS = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_from0toS.value;
        V_flap_takeoff_fromStoA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromStoA.value;
        n_flap_takeoff_fromStoA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromStoA.value;
        V_flap_takeoff_fromAtoF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromAtoF.value;
        n_flap_takeoff_fromAtoF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromAtoF.value;
        V_flap_takeoff_fromFto0 = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromFto0.value;
        n_flap_takeoff_fromFto0 = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromFto0.value;

        superposition = figure(46); 
        hold on
        grid on 
        grid minor
        
        plot(V_flap_takeoff_from0toS, n_flap_takeoff_from0toS, '-r', 'LineWidth', 1);
        plot(V_flap_takeoff_fromStoA, n_flap_takeoff_fromStoA, '-r', 'LineWidth', 1);
        plot(V_flap_takeoff_fromAtoF, n_flap_takeoff_fromAtoF, '-r', 'LineWidth', 1);
        plot(V_flap_takeoff_fromFto0, n_flap_takeoff_fromFto0, '-r', 'LineWidth', 1, 'DisplayName','Takeoff');
        % ----------------------------------------------------------------------------------------------------
        plot(V_flap_takeoff_from0toS(end), n_flap_takeoff_from0toS(end), '.k', 'MarkerSize', 10);
        text(V_flap_takeoff_from0toS(end), n_flap_takeoff_from0toS(end), 'S ', 'FontSize',    6);
        plot(V_flap_takeoff_fromStoA(end), n_flap_takeoff_fromStoA(end), '.k', 'MarkerSize', 10);
        text(V_flap_takeoff_fromStoA(end), n_flap_takeoff_fromStoA(end), 'A ', 'FontSize',    6);
        plot(V_flap_takeoff_fromAtoF(end), n_flap_takeoff_fromAtoF(end), '.k', 'MarkerSize', 10);
        text(V_flap_takeoff_fromAtoF(end), n_flap_takeoff_fromAtoF(end), 'F ', 'FontSize',    6);
        
        xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
        ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
        title("V~-~n diagram per ", Reg, "Interpreter", "latex")
        
    case 'Case 2'
        if abs(n_g_plus) > nmax
            
            V_flap_takeoff_from0toS  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_from0toS.value;
            n_flap_takeoff_from0toS  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_from0toS.value;
            V_flap_takeoff_fromStoA  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromStoA.value;
            n_flap_takeoff_fromStoA  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromStoA.value;
            V_flap_takeoff_fromAtoF1 = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromAtoF1.value;
            n_flap_takeoff_fromAtoF1 = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromAtoF1.value;
            V_flap_takeoff_fromF1toF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromF1toF.value;
            n_flap_takeoff_fromF1toF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromF1toF.value;
            V_flap_takeoff_fromFto0  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromFto0.value;
            n_flap_takeoff_fromFto0  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromFto0.value;
            
            figure(46);  
            hold on
            grid on 
            grid minor

            plot(V_flap_takeoff_from0toS, n_flap_takeoff_from0toS, '-r', 'LineWidth', 1);
            plot(V_flap_takeoff_fromStoA, n_flap_takeoff_fromStoA, '-r', 'LineWidth', 1);
            plot(V_flap_takeoff_fromAtoF1, n_flap_takeoff_fromAtoF1, '-r', 'LineWidth', 1);
            plot(V_flap_takeoff_fromF1toF, n_flap_takeoff_fromF1toF, '-r', 'LineWidth', 1);
            plot(V_flap_takeoff_fromFto0, n_flap_takeoff_fromFto0, '-r', 'LineWidth', 1, 'DisplayName','Takeoff');
            % ----------------------------------------------------------------------------------------------------
            plot(V_flap_takeoff_from0toS(end) , n_flap_takeoff_from0toS(end) , '.k', 'MarkerSize', 10);
            text(V_flap_takeoff_from0toS(end) , n_flap_takeoff_from0toS(end) , 'S ', 'FontSize',    6);
            plot(V_flap_takeoff_fromStoA(end) , n_flap_takeoff_fromStoA(end) , '.k', 'MarkerSize', 10);
            text(V_flap_takeoff_fromStoA(end) , n_flap_takeoff_fromStoA(end) , 'A ', 'FontSize',    6);
            plot(V_flap_takeoff_fromAtoF1(end), n_flap_takeoff_fromAtoF1(end), '.k', 'MarkerSize', 10);
            text(V_flap_takeoff_fromAtoF1(end), n_flap_takeoff_fromAtoF1(end), 'F1', 'FontSize',    6);
            plot(V_flap_takeoff_fromF1toF(end), n_flap_takeoff_fromF1toF(end), '.k', 'MarkerSize', 10);
            text(V_flap_takeoff_fromF1toF(end), n_flap_takeoff_fromF1toF(end), 'F ', 'FontSize',    6);
            
        
        elseif abs(n_g_plus) < nmax     
            
            V_flap_takeoff_from0toS  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_from0toS.value;
            n_flap_takeoff_from0toS  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_from0toS.value;
            V_flap_takeoff_fromStoA  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromStoA.value;
            n_flap_takeoff_fromStoA  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromStoA.value;
            V_flap_takeoff_fromAtoF  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromAtoF.value;
            n_flap_takeoff_fromAtoF  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromAtoF.value;
            V_flap_takeoff_fromFto0  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.V_fromFto0.value;
            n_flap_takeoff_fromFto0  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.final_envelope.n_fromFto0.value;
            
            figure(46);  
            hold on
            grid on 
            grid minor

            plot(V_flap_takeoff_from0toS, n_flap_takeoff_from0toS, '-r', 'LineWidth', 1);
            plot(V_flap_takeoff_fromStoA, n_flap_takeoff_fromStoA, '-r', 'LineWidth', 1);
            plot(V_flap_takeoff_fromAtoF, n_flap_takeoff_fromAtoF, '-r', 'LineWidth', 1);
            plot(V_flap_takeoff_fromFto0, n_flap_takeoff_fromFto0, '-r', 'LineWidth', 1, 'DisplayName','Takeoff');
            % ----------------------------------------------------------------------------------------------------
            plot(V_flap_takeoff_from0toS(end), n_flap_takeoff_from0toS(end), '.k', 'MarkerSize', 10);
            text(V_flap_takeoff_from0toS(end), n_flap_takeoff_from0toS(end), 'S ', 'FontSize',    6);
            plot(V_flap_takeoff_fromStoA(end), n_flap_takeoff_fromStoA(end), '.k', 'MarkerSize', 10);
            text(V_flap_takeoff_fromStoA(end), n_flap_takeoff_fromStoA(end), 'A ', 'FontSize',    6);
            plot(V_flap_takeoff_fromAtoF(end), n_flap_takeoff_fromAtoF(end), '.k', 'MarkerSize', 10);
            text(V_flap_takeoff_fromAtoF(end), n_flap_takeoff_fromAtoF(end), 'F ', 'FontSize',    6);
            
        end
end

% FLAPS - LANDING
switch (Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.case.value)
    case 'Case 1'
        
        V_flap_landing_from0toS = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_from0toS.value;
        n_flap_landing_from0toS = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_from0toS.value;
        V_flap_landing_fromStoA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromStoA.value;
        n_flap_landing_fromStoA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromStoA.value;
        V_flap_landing_fromAtoF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromAtoF.value;
        n_flap_landing_fromAtoF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromAtoF.value;
        V_flap_landing_fromFto0 = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromFto0.value;
        n_flap_landing_fromFto0 = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromFto0.value;
        
        figure(46);  
        hold on
        grid on 
        grid minor
        
        plot(V_flap_landing_from0toS, n_flap_landing_from0toS, '-b', 'LineWidth', 1);
        plot(V_flap_landing_fromStoA, n_flap_landing_fromStoA, '-b', 'LineWidth', 1);
        plot(V_flap_landing_fromAtoF, n_flap_landing_fromAtoF, '-b', 'LineWidth', 1);
        plot(V_flap_landing_fromFto0, n_flap_landing_fromFto0, '-b', 'LineWidth', 1, 'DisplayName','Landing');
        % ----------------------------------------------------------------------------------------------------
        plot(V_flap_landing_from0toS(end), n_flap_landing_from0toS(end), '.k', 'MarkerSize', 10);
        text(V_flap_landing_from0toS(end), n_flap_landing_from0toS(end), 'S ', 'FontSize',    6);
        plot(V_flap_landing_fromStoA(end), n_flap_landing_fromStoA(end), '.k', 'MarkerSize', 10);
        text(V_flap_landing_fromStoA(end), n_flap_landing_fromStoA(end), 'A ', 'FontSize',    6);
        plot(V_flap_landing_fromAtoF(end), n_flap_landing_fromAtoF(end), '.k', 'MarkerSize', 10);
        text(V_flap_landing_fromAtoF(end), n_flap_landing_fromAtoF(end), 'F ', 'FontSize',    6);
        
    case 'Case 2'
        if abs(n_g_plus) > nmax
            
            V_flap_landing_from0toS  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_from0toS.value;
            n_flap_landing_from0toS  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_from0toS.value;
            V_flap_landing_fromStoA  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromStoA.value;
            n_flap_landing_fromStoA  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromStoA.value;
            V_flap_landing_fromAtoF1 = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromAtoF1.value;
            n_flap_landing_fromAtoF1 = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromAtoF1.value;
            V_flap_landing_fromF1toF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromF1toF.value;
            n_flap_landing_fromF1toF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromF1toF.value;
            V_flap_landing_fromFto0  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromFto0.value;
            n_flap_landing_fromFto0  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromFto0.value;
            
            figure(46);  
            hold on
            grid on 
            grid minor

            plot(V_flap_landing_from0toS, n_flap_landing_from0toS, '-b', 'LineWidth', 1);
            plot(V_flap_landing_fromStoA, n_flap_landing_fromStoA, '-b', 'LineWidth', 1);
            plot(V_flap_landing_fromAtoF1, n_flap_landing_fromAtoF1, '-b', 'LineWidth', 1);
            plot(V_flap_landing_fromF1toF, n_flap_landing_fromF1toF, '-b', 'LineWidth', 1);
            plot(V_flap_landing_fromFto0, n_flap_landing_fromFto0, '-b', 'LineWidth', 1, 'DisplayName','Landing');
            % ----------------------------------------------------------------------------------------------------
            plot(V_flap_landing_from0toS(end) , n_flap_landing_from0toS(end) , '.k', 'MarkerSize', 10);
            text(V_flap_landing_from0toS(end) , n_flap_landing_from0toS(end) , 'S ', 'FontSize',    6);
            plot(V_flap_landing_fromStoA(end) , n_flap_landing_fromStoA(end) , '.k', 'MarkerSize', 10);
            text(V_flap_landing_fromStoA(end) , n_flap_landing_fromStoA(end) , 'A ', 'FontSize',    6);
            plot(V_flap_landing_fromAtoF1(end), n_flap_landing_fromAtoF1(end), '.k', 'MarkerSize', 10);
            text(V_flap_landing_fromAtoF1(end), n_flap_landing_fromAtoF1(end), 'F1', 'FontSize',    6);
            plot(V_flap_landing_fromF1toF(end), n_flap_landing_fromF1toF(end), '.k', 'MarkerSize', 10);
            text(V_flap_landing_fromF1toF(end), n_flap_landing_fromF1toF(end), 'F ', 'FontSize',    6);
            
        elseif abs(n_g_plus) < nmax     
            
            V_flap_landing_from0toS  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_from0toS.value;
            n_flap_landing_from0toS  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_from0toS.value;
            V_flap_landing_fromStoA  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromStoA.value;
            n_flap_landing_fromStoA  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromStoA.value;
            V_flap_landing_fromAtoF  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromAtoF.value;
            n_flap_landing_fromAtoF  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromAtoF.value;
            V_flap_landing_fromFto0  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.V_fromFto0.value;
            n_flap_landing_fromFto0  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.final_envelope.n_fromFto0.value;
            
            figure(46);  
            hold on
            grid on 
            grid minor

            plot(V_flap_landing_from0toS, n_flap_landing_from0toS, '-b', 'LineWidth', 1);
            plot(V_flap_landing_fromStoA, n_flap_landing_fromStoA, '-b', 'LineWidth', 1);
            plot(V_flap_landing_fromAtoF, n_flap_landing_fromAtoF, '-b', 'LineWidth', 1);
            plot(V_flap_landing_fromFto0, n_flap_landing_fromFto0, '-b', 'LineWidth', 1, 'DisplayName','Landing');
            % ----------------------------------------------------------------------------------------------------
            plot(V_flap_landing_from0toS(end), n_flap_landing_from0toS(end), '.k', 'MarkerSize', 10);
            text(V_flap_landing_from0toS(end), n_flap_landing_from0toS(end), 'S ', 'FontSize',    6);
            plot(V_flap_landing_fromStoA(end), n_flap_landing_fromStoA(end), '.k', 'MarkerSize', 10);
            text(V_flap_landing_fromStoA(end), n_flap_landing_fromStoA(end), 'A ', 'FontSize',    6);
            plot(V_flap_landing_fromAtoF(end), n_flap_landing_fromAtoF(end), '.k', 'MarkerSize', 10);
            text(V_flap_landing_fromAtoF(end), n_flap_landing_fromAtoF(end), 'F ', 'FontSize',    6);
            
        end
end

% CLEAN
switch (Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Straight_flight.value)
    case 'Case 1'
        
        if max(n_gust_cruise_plus) > nmax
            
            V_flap_clean_from0toS                   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_from0toS.value;
            n_flap_clean_from0toS                   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_from0toS.value;
            V_flap_clean_Positive_stall_speed       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value;
            n_flap_clean_Positive_stall_load_factor = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value;
            V_flap_clean_fromA1toC1                 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromA1toC1.value;
            n_flap_clean_fromA1toC1                 = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromA1toC1.value;
            V_flap_clean_fromC1toC                  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromC1toC.value;
            n_flap_clean_fromC1toC                  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromC1toC.value;
            V_flap_clean_fromCtoC2                  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromCtoC2.value;
            n_flap_clean_fromCtoC2                  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromCtoC2.value;
            V_flap_clean_fromC2toD                  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromC2toD.value;
            n_flap_clean_fromC2toD                  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromC2toD.value;
            V_flap_clean_fromDto0                   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromDto0.value;
            n_flap_clean_fromDto0                   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromDto0.value;
            
            figure(46);  
            hold on
            grid on 
            grid minor

            plot(V_flap_clean_from0toS, n_flap_clean_from0toS, '-k', 'LineWidth', 1);
            plot(V_flap_clean_Positive_stall_speed, n_flap_clean_Positive_stall_load_factor, '-k', 'LineWidth', 1);
            plot(V_flap_clean_fromA1toC1, n_flap_clean_fromA1toC1, '-k', 'LineWidth', 1);
            plot(V_flap_clean_fromC1toC, n_flap_clean_fromC1toC, '-k', 'LineWidth', 1);
            plot(V_flap_clean_fromCtoC2, n_flap_clean_fromCtoC2, '-k', 'LineWidth', 1); 
            plot(V_flap_clean_fromC2toD, n_flap_clean_fromC2toD, '-k', 'LineWidth', 1);   
            plot(V_flap_clean_fromDto0, n_flap_clean_fromDto0, '-k', 'LineWidth', 1, 'DisplayName','Clean'); 
            
        elseif max(n_gust_cruise_plus) < nmax
            
            V_flap_clean_from0toS                   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_from0toS.value;
            n_flap_clean_from0toS                   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_from0toS.value;
            V_flap_clean_Positive_stall_speed       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value;
            n_flap_clean_Positive_stall_load_factor = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value;
            V_flap_clean_fromA1toC                  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromA1toC.value;
            n_flap_clean_fromA1toC                  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromA1toC.value;
            V_flap_clean_fromCtoD                   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromCtoD.value;
            n_flap_clean_fromCtoD                   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromCtoD.value;
            V_flap_clean_fromDto0                   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromDto0.value;
            n_flap_clean_fromDto0                   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromDto0.value;
            
            figure(46);  
            hold on
            grid on 
            grid minor

            plot(V_flap_clean_from0toS, n_flap_clean_from0toS, '-k', 'LineWidth', 1);
            plot(V_flap_clean_Positive_stall_speed, n_flap_clean_Positive_stall_load_factor, '-k', 'LineWidth', 1);
            plot(V_flap_clean_fromA1toC, n_flap_clean_fromA1toC, '-k', 'LineWidth', 1);
            plot(V_flap_clean_fromCtoD, n_flap_clean_fromCtoD, '-k', 'LineWidth', 1);   
            plot(V_flap_clean_fromDto0, n_flap_clean_fromDto0, '-k', 'LineWidth', 1, 'DisplayName','Clean');   
            
        end
        
    case 'Case 2'
        
        V_flap_clean_from0toS                   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_from0toS.value;
        n_flap_clean_from0toS                   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_from0toS.value;
        V_flap_clean_Positive_stall_speed       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value;
        n_flap_clean_Positive_stall_load_factor = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value;
        V_flap_clean_fromCtoA2                  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromCtoA2.value;
        n_flap_clean_fromCtoA2                  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromCtoA2.value;
        V_flap_clean_fromA2toD                  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromA2toD.value;
        n_flap_clean_fromA2toD                  = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromA2toD.value;
        V_flap_clean_fromDto0                   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.V_fromDto0.value;
        n_flap_clean_fromDto0                   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.n_fromDto0.value;
        
        figure(46);  
        hold on
        grid on 
        grid minor

        plot(V_flap_clean_from0toS, n_flap_clean_from0toS, '-k', 'LineWidth', 1);
        plot(V_flap_clean_Positive_stall_speed, n_flap_clean_Positive_stall_load_factor, '-k', 'LineWidth', 1);
        plot(V_flap_clean_fromCtoA2, n_flap_clean_fromCtoA2, '-k', 'LineWidth', 1);
        plot(V_flap_clean_fromA2toD, n_flap_clean_fromA2toD, '-k', 'LineWidth', 1);   
        plot(V_flap_clean_fromDto0, n_flap_clean_fromDto0, '-k', 'LineWidth', 1, 'DisplayName','Clean');   
         
end

% EXPORT FIGURE
exportgraphics(superposition, 'Superposition.pdf', 'ContentType', 'vector')
exportgraphics(superposition, 'Superposition.png', 'ContentType', 'vector')

% Saving figures inside correct folder
fprintf('Saving Superposition.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile Superposition.pdf Output
movefile Superposition.png Output 
% -----------------------------------------------------------------


