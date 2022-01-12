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

% CALCULATION OF THE LOAD FACTORS VECTOR 
n_flaps_vector = calcnflap(obj, nmax);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value = n_flaps_vector;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.Attributes.unit = "g's";

% CALCULATION OF VS - CLEAN STALL SPEED
n1          = 1.0; 
rho0        = Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value;
Mass  
g           = Aircraft.Constants.g.value;
S           = Aircraft.Geometry.Wing.S.value;
WS          = ( Mass * g ) / S;
CLmax_clean = Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value;
VS          = calcvs(obj, rho0, WS, CLmax_clean, n1);

% CALCULATION OF VS1 - FLAPS DEPLOYED STALL SPEED
CLmax_takeoff = Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_takeoff.value;
VS1        = calcvs(obj, rho0, WS, CLmax_takeoff, n1);

% EVALUATION OF VF - FLAPS DEPLOYED AIRSPEED 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.VF.value = calcnVF(obj, VS, VS1);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.VF.Attributes.unit = "m/s";

% STALLING SPEED VECTOR
Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value = calcvs(obj, rho0, WS, CLmax_takeoff, n_flaps_vector);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.Attributes.unit = "m/s";

% EVALUATION OF STALL AND FLAP MANOEUVRING POINT 
% POINT S 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value = VS1;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.value = 1.0;
nS = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.Attributes.unit = "g's"; 
% POINT A
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.VA.value = Vstall(WS, rho0, CLmax_takeoff, nmax);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.VA.Attributes.unit = "m/s"; 
VA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.VA.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.nA.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
nA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.nA.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.nA.Attributes.unit = "g's"; 
% POINT F
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VF.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.Attributes.unit = "m/s"; 
VF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
nF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.Attributes.unit = "g's"; 

% VECTOR OF STALL AIRSPEED TO PLOT THE ENVELOPE 
Reg            = Aircraft.Certification.Regulation.value;
Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;
VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value;
VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value;
npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value;

n_from0toS     = linspace(0.0, nS, numb)';
V_from0toS     = VS*ones(numb, 1);

n_fromStoA     = linspace(nS, nA, numb)';
V_fromStoA     = Vstall(WS, rho0, CLmax_takeoff, n_fromStoA);

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

Aircraft.Certification.Regulation.SubpartC.Flapsloads.Flaps_envelope.value = flaps_envelope;

% Saving figures inside correct folder
fprintf('Saving flapsenvelopediagramtakeoff.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile flapsenvelopediagramtakeoff.pdf Output
movefile flapsenvelopediagramtakeoff.png Output
% -----------------------------------------------------------------

% FLAPS DEPLOYED GUST ENVELOPE 
VF            = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
rho_operative = Aircraft.Certification.ISA_Condition.Operative_ceiling.rho0.value;
MAC           = Aircraft.Geometry.Wing.mac.value; 
CLalfa_rad    = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
g             = Aircraft.Constants.g.value;
Ude           = 7.62; % Gust magnitude

% CALCULATION OF THE MASS FACTOR
mu_g = calcmug(obj, WS, MAC, CLalfa_rad, rho0, g); 

% GUST ALLEVIATION FACTOR 
Kg   = calckg(obj, mu_g);

% CALCULATION OF THE GUST LOAD FACTOR AT V = VF 
nGUST_plus  = @(V) 1.0 + V .* ((0.5 * rho0 * CLalfa_rad * Kg * Ude) / ( WS ));
nGUST_minus = @(V) 1.0 - V .* ((0.5 * rho0 * CLalfa_rad * Kg * Ude) / ( WS ));

V_gust      = linspace(0.0, VF, numb)'; 

n_gust_plus  = nGUST_plus(V_gust);
n_gust_minus = nGUST_minus(V_gust);

% STORE VALUES 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.mu_g.value                  = mu_g;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.mu_g.Attributes.unit        = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.Kg.value                    = Kg;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.Kg.Attributes.unit          = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.Vgust.value                 = V_gust;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.ngust_plus.value            = n_gust_plus;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.ngust_plus.Attributes.unit  = "g's";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.ngust_minus.value           = n_gust_plus;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.ngust_minus.Attributes.unit = "g's";

% GUST ENVELOPE AND FLIGHT ENVELOPE SUPERPOSITION 
Reg            = Aircraft.Certification.Regulation.value;
Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;
VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value;
VF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value;
npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value;

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

Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.diagram = flapsgust_envelope;
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
npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value;
VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value;
VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value;
V_gust         = linspace(VA, VF, 1e3*numb)'; 
n_g_plus       = nGUST_plus(V_gust);
n_g_minus      = nGUST_minus(V_gust);
nS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.value;
VF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
nF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.value;
Reg            = Aircraft.Certification.Regulation.value;
Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;

syms a b c V 
a        = (rho0 * CLmax_takeoff) / (2 * WS);
b        = (Kg * Ude * CLalfa_rad * rho0)/(2 * WS);
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
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.case.value = 'Case 1';

            n_from0toS     = linspace(0.0, nS, numb)';
            V_from0toS     = VS * ones(numb, 1);
            nS             = n_from0toS(end);
            VS             = V_from0toS(end);

            n_fromStoA     = linspace(nS, nA, numb)';
            V_fromStoA     = Vstall(WS, rho0, CLmax_takeoff, n_fromStoA);
            nA             = n_fromStoA(end);
            VA             = V_fromStoA(end);

            V_fromAtoF     = linspace(VA, VF, numb)';
            n_fromAtoF     = linspace(nA, nmax, numb)';  
            nF             = n_fromAtoF(end);
            VF             = V_fromAtoF(end);      

            n_fromFto0     = linspace(nF, 0.0, numb)';
            V_fromFto0     = VF * ones(numb, 1);        

            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_from0toS.value           = n_from0toS;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_from0toS.Attributes.unit = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_from0toS.value           = V_from0toS;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_from0toS.Attributes.unit = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromStoA.value           = n_fromStoA;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromStoA.Attributes.unit = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromStoA.value           = V_fromStoA;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromStoA.Attributes.unit = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromAtoF.value           = n_fromAtoF;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromAtoF.Attributes.unit = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromAtoF.value           = V_fromAtoF;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromAtoF.Attributes.unit = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromFto0.value           = n_fromFto0;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromFto0.Attributes.unit = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromFto0.value           = V_fromFto0;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromFto0.Attributes.unit = "m/s"; 
            
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

            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.diagram = final_envelope;

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
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.case.value = 'Case 2';

            VA = VA; 
            nA = nA;

            if abs(n_g_plus) > nmax 

                n_from0toS     = linspace(0.0, nS, numb)';
                V_from0toS     = VS * ones(numb, 1);
                nS             = n_from0toS(end);
                VS             = V_from0toS(end);

                n_fromStoA     = linspace(nS, nA, numb)';
                V_fromStoA     = Vstall(WS, rho0, CLmax_takeoff, n_fromStoA);
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

                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_from0toS.value            = n_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_from0toS.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_from0toS.value            = V_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_from0toS.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromStoA.value            = n_fromStoA;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromStoA.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromStoA.value            = V_fromStoA;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromStoA.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromAtoF1.value           = n_fromAtoF1;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromAtoF1.Attributes.unit = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromAtoF1.value           = V_fromAtoF1;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromAtoF1.Attributes.unit = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromF1toF.value           = n_fromF1toF;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromF1toF.Attributes.unit = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromF1toF.value           = V_fromF1toF;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromF1toF.Attributes.unit = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromFto0.value            = n_fromFto0;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromFto0.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromFto0.value            = V_fromFto0;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromFto0.Attributes.unit  = "m/s"; 
                
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

                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.diagram = final_envelope;

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
                nA             = n_fromStoA(end);
                VA             = V_fromStoA(end);

                V_fromAtoF = linspace(VA, VF, numb)';
                n_fromAtoF = nmax * ones(numb, 1);

                V_fromFto0  = VF * ones(numb, 1);
                n_fromFto0  = linspace(nF, 0.0, numb)';
                
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_from0toS.value            = n_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_from0toS.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_from0toS.value            = V_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_from0toS.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromStoA.value            = n_fromStoA;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromStoA.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromStoA.value            = V_fromStoA;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromStoA.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromAtoF.value            = n_fromAtoF;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromAtoF.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromAtoF.value            = V_fromAtoF;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromAtoF.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromFto0.value            = n_fromFto0;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromFto0.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromFto0.value            = V_fromFto0;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromFto0.Attributes.unit  = "m/s"; 

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

                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.diagram = final_envelope;

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
Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value = n_flaps_vector;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.Attributes.unit = "g's";

% CALCULATION OF VS - CLEAN STALL SPEED
n1          = 1.0; 
rho0        = Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value;
Mass  
g           = Aircraft.Constants.g.value;
S           = Aircraft.Geometry.Wing.S.value;
WS          = ( Mass * g ) / S;
CLmax_clean = Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value;
VS          = calcvs(obj, rho0, WS, CLmax_clean, n1);

% CALCULATION OF VS1 - FLAPS DEPLOYED STALL SPEED
CLmax_landing = Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_landing.value;
VS1        = calcvs(obj, rho0, WS, CLmax_landing, n1);

% EVALUATION OF VF - FLAPS DEPLOYED AIRSPEED 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.VF.value = calcnVF(obj, VS, VS1);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.VF.Attributes.unit = "m/s";

% STALLING SPEED VECTOR
Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value = calcvs(obj, rho0, WS, CLmax_landing, n_flaps_vector);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.Attributes.unit = "m/s";

% EVALUATION OF STALL AND FLAP MANOEUVRING POINT 
% POINT S 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value = VS1;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.value = 1.0;
nS = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.Attributes.unit = "g's"; 
% POINT A
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.VA.value = Vstall(WS, rho0, CLmax_landing, nmax);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.VA.Attributes.unit = "m/s"; 
VA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.VA.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.nA.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
nA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.nA.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.nA.Attributes.unit = "g's"; 
% POINT F
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VF.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.Attributes.unit = "m/s"; 
VF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
nF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.Attributes.unit = "g's"; 

% VECTOR OF STALL AIRSPEED TO PLOT THE ENVELOPE 
Reg            = Aircraft.Certification.Regulation.value;
Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;
VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value;
VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value;
npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value;

n_from0toS     = linspace(0.0, nS, numb)';
V_from0toS     = VS*ones(numb, 1);

n_fromStoA     = linspace(nS, nA, numb)';
V_fromStoA     = Vstall(WS, rho0, CLmax_landing, n_fromStoA);

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

Aircraft.Certification.Regulation.SubpartC.Flapsloads.Flaps_envelope.value = flaps_envelope;

% Saving figures inside correct folder
fprintf('Saving flapsenvelopediagramlanding.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile flapsenvelopediagramlanding.pdf Output
movefile flapsenvelopediagramlanding.png Output
% -----------------------------------------------------------------

% FLAPS DEPLOYED GUST ENVELOPE 
VF            = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
rho_operative = Aircraft.Certification.ISA_Condition.Operative_ceiling.rho0.value;
MAC           = Aircraft.Geometry.Wing.mac.value; 
CLalfa_rad    = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
g             = Aircraft.Constants.g.value;
Ude           = 7.62; % Gust magnitude

% CALCULATION OF THE MASS FACTOR
mu_g = calcmug(obj, WS, MAC, CLalfa_rad, rho0, g); 

% GUST ALLEVIATION FACTOR 
Kg   = calckg(obj, mu_g);

% CALCULATION OF THE GUST LOAD FACTOR AT V = VF 
nGUST_plus  = @(V) 1.0 + V .* ((0.5 * rho0 * CLalfa_rad * Kg * Ude) / ( WS ));
nGUST_minus = @(V) 1.0 - V .* ((0.5 * rho0 * CLalfa_rad * Kg * Ude) / ( WS ));

V_gust      = linspace(0.0, VF, numb)'; 

n_gust_plus  = nGUST_plus(V_gust);
n_gust_minus = nGUST_minus(V_gust);

% STORE VALUES 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.mu_g.value                  = mu_g;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.mu_g.Attributes.unit        = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.Kg.value                    = Kg;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.Kg.Attributes.unit          = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.Vgust.value                 = V_gust;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.ngust_plus.value            = n_gust_plus;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.ngust_plus.Attributes.unit  = "g's";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.ngust_minus.value           = n_gust_plus;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.ngust_minus.Attributes.unit = "g's";

% GUST ENVELOPE AND FLIGHT ENVELOPE SUPERPOSITION 
Reg            = Aircraft.Certification.Regulation.value;
Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;
VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value;
VF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value;
npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value;

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

Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.diagram = flapsgust_envelope;
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
npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value;
VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value;
VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value;
V_gust         = linspace(VA, VF, 1e3*numb)'; 
n_g_plus       = nGUST_plus(V_gust);
n_g_minus      = nGUST_minus(V_gust);
nS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.value;
VF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
nF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.value;
Reg            = Aircraft.Certification.Regulation.value;
Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;

syms a b c V 
a        = (rho0 * CLmax_landing) / (2 * WS);
b        = (Kg * Ude * CLalfa_rad * rho0)/(2 * WS);
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
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.case.value = 'Case 1';

            n_from0toS     = linspace(0.0, nS, numb)';
            V_from0toS     = VS * ones(numb, 1);
            nS             = n_from0toS(end);
            VS             = V_from0toS(end);

            n_fromStoA     = linspace(nS, nA, numb)';
            V_fromStoA     = Vstall(WS, rho0, CLmax_landing, n_fromStoA);
            nA             = n_fromStoA(end);
            VA             = V_fromStoA(end);

            V_fromAtoF     = linspace(VA, VF, numb)';
            n_fromAtoF     = linspace(nA, nmax, numb)';  
            nF             = n_fromAtoF(end);
            VF             = V_fromAtoF(end);      

            n_fromFto0     = linspace(nF, 0.0, numb)';
            V_fromFto0     = VF * ones(numb, 1);        

            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_from0toS.value           = n_from0toS;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_from0toS.Attributes.unit = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_from0toS.value           = V_from0toS;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_from0toS.Attributes.unit = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromStoA.value           = n_fromStoA;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromStoA.Attributes.unit = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromStoA.value           = V_fromStoA;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromStoA.Attributes.unit = "m/s"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromAtoF.value           = n_fromAtoF;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromAtoF.Attributes.unit = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromAtoF.value           = V_fromAtoF;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromAtoF.Attributes.unit = "m/s";
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromFto0.value           = n_fromFto0;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromFto0.Attributes.unit = "g's"; 
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromFto0.value           = V_fromFto0;
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromFto0.Attributes.unit = "m/s"; 
            
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

            Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.diagram = final_envelope;

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
            Aircraft.Certification.Regulation.SubpartC.Flapsloads.case.value = 'Case 2';

            VA = VA; 
            nA = nA;

            if abs(n_g_plus) > nmax 

                n_from0toS     = linspace(0.0, nS, numb)';
                V_from0toS     = VS * ones(numb, 1);
                nS             = n_from0toS(end);
                VS             = V_from0toS(end);

                n_fromStoA     = linspace(nS, nA, numb)';
                V_fromStoA     = Vstall(WS, rho0, CLmax_landing, n_fromStoA);
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

                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_from0toS.value            = n_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_from0toS.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_from0toS.value            = V_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_from0toS.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromStoA.value            = n_fromStoA;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromStoA.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromStoA.value            = V_fromStoA;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromStoA.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromAtoF1.value           = n_fromAtoF1;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromAtoF1.Attributes.unit = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromAtoF1.value           = V_fromAtoF1;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromAtoF1.Attributes.unit = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromF1toF.value           = n_fromF1toF;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromF1toF.Attributes.unit = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromF1toF.value           = V_fromF1toF;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromF1toF.Attributes.unit = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromFto0.value            = n_fromFto0;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromFto0.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromFto0.value            = V_fromFto0;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromFto0.Attributes.unit  = "m/s"; 
                
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

                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.diagram = final_envelope;

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
                nA             = n_fromStoA(end);
                VA             = V_fromStoA(end);

                V_fromAtoF = linspace(VA, VF, numb)';
                n_fromAtoF = nmax * ones(numb, 1);

                V_fromFto0  = VF * ones(numb, 1);
                n_fromFto0  = linspace(nF, 0.0, numb)';

                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_from0toS.value            = n_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_from0toS.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_from0toS.value            = V_from0toS;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_from0toS.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromStoA.value            = n_fromStoA;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromStoA.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromStoA.value            = V_fromStoA;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromStoA.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromAtoF.value            = n_fromAtoF;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromAtoF.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromAtoF.value            = V_fromAtoF;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromAtoF.Attributes.unit  = "m/s"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromFto0.value            = n_fromFto0;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromFto0.Attributes.unit  = "g's"; 
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromFto0.value            = V_fromFto0;
                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromFto0.Attributes.unit  = "m/s"; 
                
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

                Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.diagram = final_envelope;

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
switch (Aircraft.Certification.Regulation.SubpartC.Flapsloads.case.value)
    case 'Case 1'
        
        V_flap_takeoff_from0toS = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_from0toS.value;
        n_flap_takeoff_from0toS = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_from0toS.value;
        V_flap_takeoff_fromStoA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromStoA.value;
        n_flap_takeoff_fromStoA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromStoA.value;
        V_flap_takeoff_fromAtoF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromAtoF.value;
        n_flap_takeoff_fromAtoF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromAtoF.value;
        V_flap_takeoff_fromFto0 = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromFto0.value;
        n_flap_takeoff_fromFto0 = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromFto0.value;

        superposition = figure(46); 
        hold on
        grid on 
        grid minor
        
        plot(V_flap_takeoff_from0toS, n_flap_takeoff_from0toS, '-r', 'LineWidth', 1);
        plot(V_flap_takeoff_fromStoA, n_flap_takeoff_fromStoA, '-r', 'LineWidth', 1);
        plot(V_flap_takeoff_fromAtoF, n_flap_takeoff_fromAtoF, '-r', 'LineWidth', 1);
        plot(V_flap_takeoff_fromFto0, n_flap_takeoff_fromFto0, '-r', 'LineWidth', 1);
        xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
        ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
        title("V~-~n diagram per ", Reg, "Interpreter", "latex")
        
    case 'Case 2'
        if abs(n_g_plus) > nmax
            
            V_flap_takeoff_from0toS  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_from0toS.value;
            n_flap_takeoff_from0toS  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_from0toS.value;
            V_flap_takeoff_fromStoA  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromStoA.value;
            n_flap_takeoff_fromStoA  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromStoA.value;
            V_flap_takeoff_fromAtoF1 = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromAtoF1.value;
            n_flap_takeoff_fromAtoF1 = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromAtoF1.value;
            V_flap_takeoff_fromF1toF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromF1toF.value;
            n_flap_takeoff_fromF1toF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromF1toF.value;
            V_flap_takeoff_fromFto0  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromFto0.value;
            n_flap_takeoff_fromFto0  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromFto0.value;
            
            figure(46);  
            hold on
            grid on 
            grid minor

            plot(V_flap_takeoff_from0toS, n_flap_takeoff_from0toS, '-r', 'LineWidth', 1);
            plot(V_flap_takeoff_fromStoA, n_flap_takeoff_fromStoA, '-r', 'LineWidth', 1);
            plot(V_flap_takeoff_fromAtoF1, n_flap_takeoff_fromAtoF1, '-r', 'LineWidth', 1);
            plot(V_flap_takeoff_fromF1toF, n_flap_takeoff_fromF1toF, '-r', 'LineWidth', 1);
            plot(V_flap_takeoff_fromFto0, n_flap_takeoff_fromFto0, '-r', 'LineWidth', 1);
        
        elseif abs(n_g_plus) < nmax     
            
            V_flap_takeoff_from0toS  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_from0toS.value;
            n_flap_takeoff_from0toS  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_from0toS.value;
            V_flap_takeoff_fromStoA  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromStoA.value;
            n_flap_takeoff_fromStoA  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromStoA.value;
            V_flap_takeoff_fromAtoF  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromAtoF.value;
            n_flap_takeoff_fromAtoF  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromAtoF.value;
            V_flap_takeoff_fromFto0  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.V_fromFto0.value;
            n_flap_takeoff_fromFto0  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.takeoff.n_fromFto0.value;
            
            figure(46);  
            hold on
            grid on 
            grid minor

            plot(V_flap_takeoff_from0toS, n_flap_takeoff_from0toS, '-r', 'LineWidth', 1);
            plot(V_flap_takeoff_fromStoA, n_flap_takeoff_fromStoA, '-r', 'LineWidth', 1);
            plot(V_flap_takeoff_fromAtoF, n_flap_takeoff_fromAtoF, '-r', 'LineWidth', 1);
            plot(V_flap_takeoff_fromFto0, n_flap_takeoff_fromFto0, '-r', 'LineWidth', 1);
            
        end
end

% FLAPS - LANDING
switch (Aircraft.Certification.Regulation.SubpartC.Flapsloads.case.value)
    case 'Case 1'
        
        V_flap_landing_from0toS = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_from0toS.value;
        n_flap_landing_from0toS = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_from0toS.value;
        V_flap_landing_fromStoA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromStoA.value;
        n_flap_landing_fromStoA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromStoA.value;
        V_flap_landing_fromAtoF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromAtoF.value;
        n_flap_landing_fromAtoF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromAtoF.value;
        V_flap_landing_fromFto0 = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromFto0.value;
        n_flap_landing_fromFto0 = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromFto0.value;
        
        figure(46);  
        hold on
        grid on 
        grid minor
        
        plot(V_flap_landing_from0toS, n_flap_landing_from0toS, '-b', 'LineWidth', 1);
        plot(V_flap_landing_fromStoA, n_flap_landing_fromStoA, '-b', 'LineWidth', 1);
        plot(V_flap_landing_fromAtoF, n_flap_landing_fromAtoF, '-b', 'LineWidth', 1);
        plot(V_flap_landing_fromFto0, n_flap_landing_fromFto0, '-b', 'LineWidth', 1);
        
    case 'Case 2'
        if abs(n_g_plus) > nmax
            
            V_flap_landing_from0toS  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_from0toS.value;
            n_flap_landing_from0toS  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_from0toS.value;
            V_flap_landing_fromStoA  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromStoA.value;
            n_flap_landing_fromStoA  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromStoA.value;
            V_flap_landing_fromAtoF1 = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromAtoF1.value;
            n_flap_landing_fromAtoF1 = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromAtoF1.value;
            V_flap_landing_fromF1toF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromF1toF.value;
            n_flap_landing_fromF1toF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromF1toF.value;
            V_flap_landing_fromFto0  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromFto0.value;
            n_flap_landing_fromFto0  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromFto0.value;
            
            figure(46);  
            hold on
            grid on 
            grid minor

            plot(V_flap_landing_from0toS, n_flap_landing_from0toS, '-b', 'LineWidth', 1);
            plot(V_flap_landing_fromStoA, n_flap_landing_fromStoA, '-b', 'LineWidth', 1);
            plot(V_flap_landing_fromAtoF1, n_flap_landing_fromAtoF1, '-b', 'LineWidth', 1);
            plot(V_flap_landing_fromF1toF, n_flap_landing_fromF1toF, '-b', 'LineWidth', 1);
            plot(V_flap_landing_fromFto0, n_flap_landing_fromFto0, '-b', 'LineWidth', 1);
            
        elseif abs(n_g_plus) < nmax     
            
            V_flap_landing_from0toS  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_from0toS.value;
            n_flap_landing_from0toS  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_from0toS.value;
            V_flap_landing_fromStoA  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromStoA.value;
            n_flap_landing_fromStoA  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromStoA.value;
            V_flap_landing_fromAtoF  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromAtoF.value;
            n_flap_landing_fromAtoF  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromAtoF.value;
            V_flap_landing_fromFto0  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.V_fromFto0.value;
            n_flap_landing_fromFto0  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.landing.n_fromFto0.value;
            
            figure(46);  
            hold on
            grid on 
            grid minor

            plot(V_flap_landing_from0toS, n_flap_landing_from0toS, '-b', 'LineWidth', 1);
            plot(V_flap_landing_fromStoA, n_flap_landing_fromStoA, '-b', 'LineWidth', 1);
            plot(V_flap_landing_fromAtoF, n_flap_landing_fromAtoF, '-b', 'LineWidth', 1);
            plot(V_flap_landing_fromFto0, n_flap_landing_fromFto0, '-b', 'LineWidth', 1);
            
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
            plot(V_flap_clean_fromDto0, n_flap_clean_fromDto0, '-k', 'LineWidth', 1);      
            
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
            plot(V_flap_clean_fromDto0, n_flap_clean_fromDto0, '-k', 'LineWidth', 1);   
            
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
        plot(V_flap_clean_fromDto0, n_flap_clean_fromDto0, '-k', 'LineWidth', 1);   
         
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


% CLEAN

% switch (Straight_flight_Case)
%     % CASE 1: VA greater than the intercept
%     case 'Case 1'
%         for i = 1:length(Solution)
%             new_VA = cast(Solution(i), 'double');
% %             if abs(new_VA) > VA
% %                 VA = abs(new_VA);
% %                 nA = nGUST_plus(VA);
% %             elseif abs(new_VA) < VA
% %                 VA = VA; 
% %                 nA = nA;
% %             end
% %         end
%             if max(n_g_plus) > nmax 
%                 
%                 VA = VA;
%                 nA = nA;
% 
%                 % FROM 0 TO S
%                 n_from0toS = linspace(0.0, n1, numb)';
%                 V_from0toS = VS * ones(numb, 1);
% 
%                 % FROM S TO A 
%                 n_fromStoA     = linspace(nS, nA, numb)';
%                 V_fromStoA     = Vstall(WS, rho0, CLmax_takeoff, n_fromStoA);
% 
%                 % FROM A TO F1 
%                 tol          = 1e-3;
%                 for i = 1:length(V_g)
%                     x = n_g_plus(i);
%                     y = x - nmax;
%                     if abs(y) < tol
%                         row = i; 
%                         VF1 = V_g(row);
%                     end
%                 end
%                 V_fromAtoF1 = linspace(VA, VF1, numb)';
%                 n_fromAtoF1 = nmax * ones(length(V_fromAtoF1),1);
%                 nF1          = nmax;
% 
%                 % FROM F1 TO F
%                 V_fromF1toF = linspace(VF1, VF, numb)'; 
%                 n_fromF1toF = nGUST_plus(V_fromF1toF);
%                 nF          = n_fromF1toF(end);
%                 VF          = V_fromF1toF(end);
% 
%                 % FROM F TO 0
%                 V_fromFto0 = VF * ones(numb, 1);
%                 n_fromFto0 = linspace(nF, 0.0, numb)';
% 
%                 final_envelope = figure;
%                 hold on
%                 grid on 
%                 grid minor
%                 ylim([-0.5 nmax+0.5])
%                 xlim([0 VF+10])
% 
%     %             plot(V_gust,     n_gust_plus,  ':k', 'LineWidth', 0.2)
%     %             plot(V_gust,     n_gust_minus, ':k', 'LineWidth', 0.2)
%     %             plot(VSpos,      npos,         ':r', 'LineWidth', 0.2)
%                 plot(V_from0toS,  n_from0toS,   '-r', 'LineWidth',   1)
%                 plot(V_fromStoA,  n_fromStoA,   '-r', 'LineWidth',   1)
%                 plot(V_fromAtoF1, n_fromAtoF1,  '-b', 'LineWidth',   1)
%                 plot(V_fromF1toF, n_fromF1toF,  '-b', 'LineWidth',   1)
%                 plot(V_fromFto0,  n_fromFto0,   '-b', 'LineWidth',   1)
% 
%                 plot(V_fromStoA(1),  n_fromStoA(1),  'k.', 'MarkerSize', 10)
%                 plot(V_fromAtoF1(1), n_fromAtoF1(1), 'k.', 'MarkerSize', 10)
%                 plot(V_fromF1toF(1), n_fromF1toF(1), 'k.', 'MarkerSize', 10)
%                 plot(V_fromFto0(1),  n_fromFto0(1),  'k.', 'MarkerSize', 10)
% 
%                 xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
%                 ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
%                 title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
%                 text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
%                 text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
%                 text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)   
% 
%                 exportgraphics(final_envelope,'flaps_final_envelopediagram.pdf','ContentType','vector')
%                 exportgraphics(final_envelope,'flaps_final_envelopediagram.png','ContentType','vector')
% 
%                 Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.diagram = final_envelope;
% 
%                 % Saving figures inside correct folder
%                 fprintf('Saving flaps_final_envelopediagram.pdf in: ');
%                 fprintf('\n'); 
%                 fprintf('%s\n', SaveFolder);
%                 % Moving file inside correct folder
%                 movefile flaps_final_envelopediagram.pdf Output
%                 movefile flaps_final_envelopediagram.png Output
%             
%         elseif (max(n_g_plus) < nmax)  
% 
%             VA = new_VA;
%             nA = nGUST_plus(VA);
%             
%             % FROM 0 TO S
%             n_from0toS = linspace(0.0, n1, numb)';
%             V_from0toS = VS * ones(numb, 1);
%             
%             % FROM S TO A 
%             n_fromStoA     = linspace(nS, nA, numb)';
%             V_fromStoA     = Vstall(WS, rho0, CLmax_takeoff, n_fromStoA);
%             
%             % FROM A TO F
%             V_fromAtoF = linspace(VA, VF, numb)'; 
%             n_fromAtoF = nmax * ones(numb, 1);
%             nF          = n_fromAtoF(end);
%             VF          = V_fromAtoF(end);
%             
%             % FROM F TO 0
%             V_fromFto0 = VF * ones(numb, 1);
%             n_fromFto0 = linspace(nF, 0.0, numb)';
%             
%             final_envelope = figure;
%             hold on
%             grid on 
%             grid minor
%             ylim([-0.5 nmax+0.5])
%             xlim([0 VF+10])
% 
% %             plot(V_gust,     n_gust_plus,  ':k', 'LineWidth', 0.2)
% %             plot(V_gust,     n_gust_minus, ':k', 'LineWidth', 0.2)
% %             plot(VSpos,      npos,         ':r', 'LineWidth', 0.2)
%             plot(V_from0toS,  n_from0toS,   '-r', 'LineWidth',   1)
%             plot(V_fromStoA,  n_fromStoA,   '-r', 'LineWidth',   1)
%             plot(V_fromAtoF,  n_fromAtoF,   '-b', 'LineWidth',   1)
%             plot(V_fromFto0,  n_fromFto0,   '-b', 'LineWidth',   1)
% 
%             plot(V_fromStoA(1),  n_fromStoA(1),  'k.', 'MarkerSize', 10)
%             plot(V_fromAtoF(1),  n_fromAtoF(1),  'k.', 'MarkerSize', 10)
%             plot(V_fromFto0(1),  n_fromFto0(1),  'k.', 'MarkerSize', 10)
% 
%             xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
%             ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
%             title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
%             text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
%             text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
%             text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)   
%             
%             exportgraphics(final_envelope,'flaps_final_envelopediagram.pdf','ContentType','vector')
%             exportgraphics(final_envelope,'flaps_final_envelopediagram.png','ContentType','vector')
% 
%             Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.diagram = final_envelope;
% 
%             % Saving figures inside correct folder
%             fprintf('Saving flaps_final_envelopediagram.pdf in: ');
%             fprintf('\n'); 
%             fprintf('%s\n', SaveFolder);
%             % Moving file inside correct folder
%             movefile flaps_final_envelopediagram.pdf Output
%             movefile flaps_final_envelopediagram.png Output
%             end
%         end
%         
%     % CASE 1: VA lower than the intercept    
%     case 'Case 2'
%         syms a b c V 
%         a        = (rho0 * CLmax_takeoff) / (2 * WS);
%         b        = (Kg * Ude * CLalfa_rad * rho0)/(2 * WS);
%         c        = 1;
%         eqn      = a*V^2 - b*V - c;
%         Solution = vpasolve(eqn, V);
% 
%         for i = 1:length(Solution)
%            new_VA = cast(Solution(i), 'double');
%            if abs(new_VA) > VA
%                VA = abs(new_VA);
%                nA = nGUST_plus(VA);
%            elseif abs(new_VA) < VA
%                VA = VA; 
%                nA = nA;
%            end
%         end
% 
%         % POINT A
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.VA.value           = Vstall(WS, rho0, CLmax_takeoff, nA);
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.VA.Attributes.unit = "m/s"; 
%         VA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.VA.value;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.nA.value           = nA;
%         nA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.nA.value;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.nA.Attributes.unit = "g's"; 
% 
%         n_from0toS     = linspace(0.0, nS, numb);
%         V_from0toS     = VS*ones(numb, 1);
% 
%         n_fromStoA     = linspace(nS, nA, numb);
%         V_fromStoA     = Vstall(WS, rho0, CLmax_takeoff, n_fromStoA);
% 
%         n_fromAtoF     = [nA nF];
%         V_fromAtoF     = [VA VF];
% 
%         n_fromFto0     = [nF 0.0];
%         V_fromFto0     = [VF VF];
% 
%         final_envelope = figure;
%         hold on
%         grid on 
%         grid minor
%         ylim([-0.5 nmax+0.5])
%         xlim([0 VF+10])
%         plot(V_gust, n_gust_plus,    ':k', 'LineWidth', 0.2)
%         plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth',   1)
%         plot(V_fromAtoF, n_fromAtoF, '-r', 'LineWidth',   1)
%         plot(V_fromFto0, n_fromFto0, '-r', 'LineWidth',   1)
%         plot(V_from0toS, n_from0toS, '-r', 'LineWidth',   1) 
%  
%         plot(V_fromStoA(1),  n_fromStoA(1),  'k.', 'MarkerSize', 10)
%         plot(V_fromAtoF(1),  n_fromAtoF(1),  'k.', 'MarkerSize', 10)
%         plot(V_fromFto0(1),  n_fromFto0(1),  'k.', 'MarkerSize', 10)
%        
%         xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
%         ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
%         title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
%         text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
%         text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
%         text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)
%         exportgraphics(flapsgust_envelope,'flaps_final_envelopediagram.pdf','ContentType','vector')
%         exportgraphics(flapsgust_envelope,'flaps_final_envelopediagram.png','ContentType','vector')
% 
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.diagram = final_envelope;
% 
%         % Saving figures inside correct folder
%         fprintf('Saving flaps_final_envelopediagram.pdf in: ');
%         fprintf('\n'); 
%         fprintf('%s\n', SaveFolder);
%         % Moving file inside correct folder
%         movefile flaps_final_envelopediagram.pdf Output
%         movefile flaps_final_envelopediagram.png Output        
% end


% switch (Straight_flight_Case)
%     % CASE 1: VA greater than the intercept
%     case 'Case 1'
%         % IMPLEMENTATION 
%         numb = 1e3;
%         nmax = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
% 
%         % CALCULATION OF THE LOAD FACTORS VECTOR 
%         n_flaps_vector = calcnflap(obj, nmax);
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value = n_flaps_vector;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.Attributes.unit = "g's";
% 
%         % CALCULATION OF VS - CLEAN STALL SPEED
%         n1    = 1.0; 
%         rho0  = Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value;
%         x     = 0;
%         Mass  = Aircraft.Weight.I_Level.W_maxTakeOff.value + x;
%         g     = Aircraft.Constants.g.value;
%         S     = Aircraft.Geometry.Wing.S.value;
%         WS    = ( Mass * g ) / S;
%         CLmax = Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value;
%         VS    = calcvs(obj, rho0, WS, CLmax, n1);
% 
%         % CALCULATION OF VS1 - FLAPS DEPLOYED STALL SPEED
%         CLmax_takeoff = Aircraft.Certification.Aerodynamic_data.Flaps.CLmax_takeoffs.value;
%         VS1        = calcvs(obj, rho0, WS, CLmax_takeoff, n1);
% 
%         % EVALUATION OF VF - FLAPS DEPLOYED AIRSPEED 
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.VF.value = calcnVF(obj, VS, VS1);
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.VF.Attributes.unit = "m/s";
% 
%         % STALLING SPEED VECTOR
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value = calcvs(obj, rho0, WS, CLmax_takeoff, n_flaps_vector);
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.Attributes.unit = "m/s";
% 
%         % EVALUATION OF STALL AND FLAP MANOEUVRING POINT 
%         % POINT S 
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value = VS1;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.Attributes.unit = "m/s"; 
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.value = 1.0;
%         nS = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.value;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.Attributes.unit = "g's"; 
%         % POINT A
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.VA.value = Vstall(WS, rho0, CLmax_takeoff, nmax);
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.VA.Attributes.unit = "m/s"; 
%         VA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.VA.value;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.nA.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
%         nA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.nA.value;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.nA.Attributes.unit = "g's"; 
%         % POINT F
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VF.value;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.Attributes.unit = "m/s"; 
%         VF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
%         nF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.value;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.Attributes.unit = "g's"; 
%         
%         % VECTOR OF STALL AIRSPEED TO PLOT THE ENVELOPE 
%         Reg            = Aircraft.Certification.Regulation.value;
%         Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;
%         VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value;
%         VF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
%         VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value;
%         npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value;
% 
%         n_from0toS     = linspace(0.0, nS, numb);
%         V_from0toS     = VS*ones(numb, 1);
% 
%         n_fromStoA     = linspace(nS, nA, numb);
%         V_fromStoA     = Vstall(WS, rho0, CLmax_takeoff, n_fromStoA);
% 
%         n_fromAtoF     = nmax*ones(numb, 1);
%         V_fromAtoF     = linspace(VA, VF, numb);
% 
%         n_fromFto0     = linspace(nF, 0.0, numb);
%         V_fromFto0     = VF*ones(numb, 1);
% 
%         flaps_envelope = figure;
%         hold on
%         grid on 
%         grid minor
%         ylim([-0.5 nmax+0.5])
%         xlim([0 VF+10])
%         plot(VSpos, npos, ':r', 'LineWidth',0.2)
%         plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth',1)
%         plot(V_fromAtoF, n_fromAtoF, '-b', 'LineWidth',1)
%         plot(V_fromFto0, n_fromFto0, '-b', 'LineWidth',1)
%         plot(V_from0toS, n_from0toS, '-b', 'LineWidth',1)
%         xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
%         ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
%         title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
%         % text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
%         text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
%         text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)
%         exportgraphics(flaps_envelope,'flapsenvelopediagram.pdf','ContentType','vector')
%         exportgraphics(flaps_envelope,'flapsenvelopediagram.png','ContentType','vector')
% 
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Flaps_envelope.value = flaps_envelope;
% 
%         % Saving figures inside correct folder
%         fprintf('Saving flapsenvelopediagram.pdf in: ');
%         fprintf('\n'); 
%         fprintf('%s\n', SaveFolder);
%         % Moving file inside correct folder
%         movefile flapsenvelopediagram.pdf Output
%         movefile flapsenvelopediagram.png Output
%         % -----------------------------------------------------------------
%         
%         % FLAPS DEPLOYED GUST ENVELOPE 
%         VF            = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
%         rho_operative = Aircraft.Certification.ISA_Condition.Operative_ceiling.rho0.value;
%         MAC           = Aircraft.Geometry.Wing.mac.value; 
%         CLalfa_rad    = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
%         g             = Aircraft.Constants.g.value;
%         Ude           = 7.62; % Gust magnitude
% 
%         % CALCULATION OF THE MASS FACTOR
%         mu_g = calcmug(obj, WS, MAC, CLalfa_rad, rho0, g); 
% 
%         % GUST ALLEVIATION FACTOR 
%         Kg   = calckg(obj, mu_g);
% 
%         % CALCULATION OF THE GUST LOAD FACTOR AT V = VF 
%         nGUST  = @(V) 1.0 + V * ((0.5 * rho0 * CLalfa_rad * Kg * Ude) / ( WS ));
%         V_gust = linspace(0.0, VF, numb); 
%         n_gust = nGUST(V_gust);
% 
%         % STORE VALUES 
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.mu_g.value           = mu_g;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.mu_g.Attributes.unit = "Non dimensional";
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.Kg.value             = Kg;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.Kg.Attributes.unit   = "Non dimensional";
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.Vgust.value          = V_gust;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.ngust.value          = n_gust;
%             
%         % GUST ENVELOPE AND FLIGHT ENVELOPE SUPERPOSITION 
%         Reg            = Aircraft.Certification.Regulation.value;
%         Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;
%         VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value;
%         VF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
%         VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value;
%         npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value;
% 
%         flapsgust_envelope = figure;
%         hold on
%         grid on 
%         grid minor
%         ylim([-0.5 nmax+0.5])
%         xlim([0 VF+10])
%         plot(V_gust, n_gust, ':k', 'LineWidth', 0.2)
%         plot(VSpos, npos, ':r', 'LineWidth', 0.2)
%         plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth',1)
%         plot(V_fromAtoF, n_fromAtoF, '-b', 'LineWidth',1)
%         plot(V_fromFto0, n_fromFto0, '-b', 'LineWidth',1)
%         plot(V_from0toS, n_from0toS, '-b', 'LineWidth',1)
%         xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
%         ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
%         title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
%         text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
%         text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
%         text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)
%         exportgraphics(flapsgust_envelope,'flaps_gust_envelopediagram.pdf','ContentType','vector')
%         exportgraphics(flapsgust_envelope,'flaps_gust_envelopediagram.png','ContentType','vector')
% 
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.diagram = flapsgust_envelope;
%         % Saving figures inside correct folder
%         fprintf('Saving flapsenvelopediagram.pdf in: ');
%         fprintf('\n'); 
%         fprintf('%s\n', SaveFolder);
%         % Moving file inside correct folder
%         movefile flaps_gust_envelopediagram.pdf Output
%         movefile flaps_gust_envelopediagram.png Output        
%             
%         % FINAL ENVELOPE WITH FLAPS DEPLOYED 
%         nmax           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
%         CLalfa_rad     = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
%         npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value;
%         VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value;
%         VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value;
%         V_g            = linspace(VS, VF, numb); 
%         n_g            = nGUST(V_g);
%         nS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.value;
%         VF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
%         nF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.value;
%         Reg            = Aircraft.Certification.Regulation.value;
%         Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;
% 
%         syms a b c V 
%         a        = (rho0 * CLmax_takeoff) / (2 * WS);
%         b        = (Kg * Ude * CLalfa_rad * rho0)/(2 * WS);
%         c        = 1;
%         eqn      = a*V^2 - b*V - c;
%         Solution = vpasolve(eqn, V);
% 
%         for i = 1:length(Solution)
%            new_VA = cast(Solution(i), 'double');
%            if abs(new_VA) > VA
%                VA1 = abs(new_VA);
%                nA1 = nGUST(VA1);
%            elseif abs(new_VA) < VA
%                VA1 = VA; 
%                nA1 = nA;
%            end
%         end
%             
%             if (max(n_gust) > nA) & (VA < VA1)
% 
%             % POINT A
%             Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA1.VA1.value           = Vstall(WS, rho0, CLmax_takeoff, nA1);
%             Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA1.VA1.Attributes.unit = "m/s"; 
%             VA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA1.VA1.value;
%             Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA1.nA1.value           = nA1;
%             nA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA1.nA1.value;
%             Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA1.nA1.Attributes.unit = "g's"; 
% 
%             n_from0toS     = linspace(0.0, nS, numb);
%             V_from0toS     = VS*ones(numb, 1);
% 
%             n_fromStoA     = linspace(nS, nA, numb);
%             V_fromStoA     = Vstall(WS, rho0, CLmax_takeoff, n_fromStoA);
%             
%             n_fromAtoA1    = nmax * ones(length(numb), 1); 
%             V_fromAtoA1    = linspace(VA, VA1, numb)'; 
%             
%             
%             
%             n_fromAtoF     = [nA nF];
%             V_fromAtoF     = [VA VF];
% 
%             n_fromFto0     = [nF 0.0];
%             V_fromFto0     = [VF VF];
% 
%             final_envelope = figure;
%             hold on
%             grid on 
%             grid minor
%             ylim([-0.5 nmax+0.5])
%             xlim([0 VF+10])
%             plot(V_gust, n_gust, ':k', 'LineWidth', 0.2)
%             plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth',1)
%             plot(V_fromAtoF, n_fromAtoF, '-b', 'LineWidth',1)
%             plot(V_fromFto0, n_fromFto0, '-b', 'LineWidth',1)
%             plot(V_from0toS, n_from0toS, '-b', 'LineWidth',1)
%             xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
%             ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
%             title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
%             text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
%             text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
%             text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)
%             exportgraphics(flapsgust_envelope,'flaps_final_envelopediagram.pdf','ContentType','vector')
%             exportgraphics(flapsgust_envelope,'flaps_final_envelopediagram.png','ContentType','vector')
% 
%             Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.diagram = final_envelope;
% 
%             % Saving figures inside correct folder
%             fprintf('Saving flaps_final_envelopediagram.pdf in: ');
%             fprintf('\n'); 
%             fprintf('%s\n', SaveFolder);
%             % Moving file inside correct folder
%             movefile flaps_final_envelopediagram.pdf Output
%             movefile flaps_final_envelopediagram.png Output
%             
%         elseif max(n_gust) < nA
%             
%             syms a b c V 
%             a        = (rho0 * CLmax_takeoff) / (2 * WS);
%             b        = (Kg * Ude * CLalfa_rad * rho0)/(2 * WS);
%             c        = 1;
%             eqn      = a*V^2 - b*V - c;
%             Solution = vpasolve(eqn, V);
% 
%             for i = 1:length(Solution)
%                new_VA = cast(Solution(i), 'double');
%                if abs(new_VA) > VA
%                    VA = abs(new_VA);
%                    nA = nGUST(VA);
%                elseif abs(new_VA) < VA
%                    VA = VA; 
%                    nA = nA;
%                end
%             end
% 
%             % POINT A
%             Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.VA.value           = Vstall(WS, rho0, CLmax_takeoff, nA);
%             Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.VA.Attributes.unit = "m/s"; 
%             VA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.VA.value;
%             Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.nA.value           = nA;
%             nA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.nA.value;
%             Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.nA.Attributes.unit = "g's"; 
% 
%             n_from0toS     = linspace(0.0, nS, numb);
%             V_from0toS     = VS*ones(numb, 1);
% 
%             n_fromStoA     = linspace(nS, nA, numb);
%             V_fromStoA     = Vstall(WS, rho0, CLmax_takeoff, n_fromStoA);
% 
%             n_fromAtoF     = [nA nF];
%             V_fromAtoF     = [VA VF];
% 
%             n_fromFto0     = [nF 0.0];
%             V_fromFto0     = [VF VF];
% 
%             final_envelope = figure;
%             hold on
%             grid on 
%             grid minor
%             ylim([-0.5 nmax+0.5])
%             xlim([0 VF+10])
%             plot(V_gust, n_gust, ':k', 'LineWidth', 0.2)
%             plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth',1)
%             plot(V_fromAtoF, n_fromAtoF, '-b', 'LineWidth',1)
%             plot(V_fromFto0, n_fromFto0, '-b', 'LineWidth',1)
%             plot(V_from0toS, n_from0toS, '-b', 'LineWidth',1)
%             xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
%             ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
%             title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
%             text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
%             text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
%             text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)
%             exportgraphics(flapsgust_envelope,'flaps_final_envelopediagram.pdf','ContentType','vector')
%             exportgraphics(flapsgust_envelope,'flaps_final_envelopediagram.png','ContentType','vector')
% 
%             Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.diagram = final_envelope;
% 
%             % Saving figures inside correct folder
%             fprintf('Saving flaps_final_envelopediagram.pdf in: ');
%             fprintf('\n'); 
%             fprintf('%s\n', SaveFolder);
%             % Moving file inside correct folder
%             movefile flaps_final_envelopediagram.pdf Output
%             movefile flaps_final_envelopediagram.png Output
%             
%         end
%         % -----------------------------------------------------------------
%     % CASE 1: VA lower than the intercept
%     case 'Case 2'
%         % IMPLEMENTATION 
%         numb = 1e3;
%         nmax = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
% 
%         % CALCULATION OF THE LOAD FACTORS VECTOR 
%         n_flaps_vector = calcnflap(obj, nmax);
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value = n_flaps_vector;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.Attributes.unit = "g's";
% 
%         % CALCULATION OF VS - CLEAN STALL SPEED
%         n1          = 1.0; 
%         rho0        = Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value;
%         x           = 0;
%         Mass        = Aircraft.Weight.I_Level.W_maxTakeOff.value + x;
%         g           = Aircraft.Constants.g.value;
%         S           = Aircraft.Geometry.Wing.S.value;
%         WS          = ( Mass * g ) / S;
%         CLmax_clean = Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value;
%         VS          = calcvs(obj, rho0, WS, CLmax_clean, n1);
% 
%         % CALCULATION OF VS1 - FLAPS DEPLOYED STALL SPEED
%         CLmax_takeoff = Aircraft.Certification.Aerodynamic_data.Flaps.CLmax_takeoffs.value;
%         VS1        = calcvs(obj, rho0, WS, CLmax_takeoff, n1);
% 
%         % EVALUATION OF VF - FLAPS DEPLOYED AIRSPEED 
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.VF.value = calcnVF(obj, VS, VS1);
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.VF.Attributes.unit = "m/s";
% 
%         % STALLING SPEED VECTOR
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value = calcvs(obj, rho0, WS, CLmax_takeoff, n_flaps_vector);
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.Attributes.unit = "m/s";
% 
%         % EVALUATION OF STALL AND FLAP MANOEUVRING POINT 
%         % POINT S 
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value = VS1;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.Attributes.unit = "m/s"; 
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.value = 1.0;
%         nS = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.value;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.Attributes.unit = "g's"; 
%         % POINT A
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.VA.value = Vstall(WS, rho0, CLmax_takeoff, nmax);
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.VA.Attributes.unit = "m/s"; 
%         VA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.VA.value;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.nA.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
%         nA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.nA.value;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.nA.Attributes.unit = "g's"; 
%         % POINT F
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VF.value;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.Attributes.unit = "m/s"; 
%         VF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
%         nF = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.value;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.Attributes.unit = "g's"; 
% 
%         % VECTOR OF STALL AIRSPEED TO PLOT THE ENVELOPE 
%         Reg            = Aircraft.Certification.Regulation.value;
%         Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;
%         VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value;
%         VF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
%         VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value;
%         npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value;
% 
%         n_from0toS     = linspace(0.0, nS, numb);
%         V_from0toS     = VS*ones(numb, 1);
% 
%         n_fromStoA     = linspace(nS, nA, numb);
%         V_fromStoA     = Vstall(WS, rho0, CLmax_takeoff, n_fromStoA);
% 
%         n_fromAtoF     = nmax*ones(numb, 1);
%         V_fromAtoF     = linspace(VA, VF, numb);
% 
%         n_fromFto0     = linspace(nF, 0.0, numb);
%         V_fromFto0     = VF*ones(numb,1);
% 
%         flaps_envelope = figure;
%         hold on
%         grid on 
%         grid minor
%         ylim([-0.5 nmax+0.5])
%         xlim([0 VF+10])
%         plot(VSpos, npos, ':r', 'LineWidth',0.2)
%         plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth',1)
%         plot(V_fromAtoF, n_fromAtoF, '-b', 'LineWidth',1)
%         plot(V_fromFto0, n_fromFto0, '-b', 'LineWidth',1)
%         plot(V_from0toS, n_from0toS, '-b', 'LineWidth',1)
%         xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
%         ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
%         title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
%         % text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
%         text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
%         text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)
%         exportgraphics(flaps_envelope,'flapsenvelopediagram.pdf','ContentType','vector')
%         exportgraphics(flaps_envelope,'flapsenvelopediagram.png','ContentType','vector')
% 
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Flaps_envelope.value = flaps_envelope;
% 
%         % Aircraft.Certification.Regulation.SubpartC.Flapsloads.Flaps_envelope.value = flapsenvelope_diagram(obj, npos, nmax, VSpos, VS, VF, Reg, Aircraft_name);
% 
%         % Saving figures inside correct folder
%         fprintf('Saving flapsenvelopediagram.pdf in: ');
%         fprintf('\n'); 
%         fprintf('%s\n', SaveFolder);
%         % Moving file inside correct folder
%         movefile flapsenvelopediagram.pdf Output
%         movefile flapsenvelopediagram.png Output
% 
%         % FLAPS DEPLOYED GUST ENVELOPE 
%         VF            = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
%         rho_operative = Aircraft.Certification.ISA_Condition.Operative_ceiling.rho0.value;
%         MAC           = Aircraft.Geometry.Wing.mac.value; 
%         CLalfa_rad    = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
%         g             = Aircraft.Constants.g.value;
%         Ude           = 7.62; % Gust magnitude
% 
%         % CALCULATION OF THE MASS FACTOR
%         mu_g = calcmug(obj, WS, MAC, CLalfa_rad, rho0, g); 
% 
%         % GUST ALLEVIATION FACTOR 
%         Kg   = calckg(obj, mu_g);
% 
%         % CALCULATION OF THE GUST LOAD FACTOR AT V = VF 
%         nGUST  = @(V) 1.0 + V * (( 0.5 * rho0 * CLalfa_rad * Kg * Ude) / (WS));
%         V_gust = linspace(0.0, VF, numb); 
%         n_gust = nGUST(V_gust);
% 
%         % STORE VALUES 
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.mu_g.value = mu_g;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.mu_g.Attributes.unit = "Non dimensional";
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.Kg.value = Kg;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.Kg.Attributes.unit = "Non dimensional";
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.Vgust.value = V_gust;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.ngust.value = n_gust;
% 
%         % GUST ENVELOPE AND FLIGHT ENVELOPE SUPERPOSITION 
%         Reg            = Aircraft.Certification.Regulation.value;
%         Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;
%         VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value;
%         VF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
%         VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value;
%         npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value;
% 
%         flapsgust_envelope = figure;
%         hold on
%         grid on 
%         grid minor
%         ylim([-0.5 nmax+0.5])
%         xlim([0 VF+10])
%         plot(V_gust, n_gust, ':k', 'LineWidth', 0.2)
%         plot(VSpos, npos, ':r', 'LineWidth', 0.2)
%         plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth',1)
%         plot(V_fromAtoF, n_fromAtoF, '-b', 'LineWidth',1)
%         plot(V_fromFto0, n_fromFto0, '-b', 'LineWidth',1)
%         plot(V_from0toS, n_from0toS, '-b', 'LineWidth',1)
%         xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
%         ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
%         title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
%         text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
%         text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
%         text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)
%         exportgraphics(flapsgust_envelope,'flaps_gust_envelopediagram.pdf','ContentType','vector')
%         exportgraphics(flapsgust_envelope,'flaps_gust_envelopediagram.png','ContentType','vector')
% 
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.diagram = flapsgust_envelope;
%         % Saving figures inside correct folder
%         fprintf('Saving flapsenvelopediagram.pdf in: ');
%         fprintf('\n'); 
%         fprintf('%s\n', SaveFolder);
%         % Moving file inside correct folder
%         movefile flaps_gust_envelopediagram.pdf Output
%         movefile flaps_gust_envelopediagram.png Output
% 
%         %% FINAL ENVELOPE WITH FLAPS DEPLOYED 
%         nmax           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
%         CLalfa_rad     = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
%         npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value;
%         VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value;
%         VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value;
%         V_g            = linspace(VS, VF, numb); 
%         n_g            = nGUST(V_g);
%         nS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.value;
%         VF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
%         nF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.value;
%         Reg            = Aircraft.Certification.Regulation.value;
%         Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;
% 
%         syms a b c V 
%         a        = (rho0 * CLmax_takeoff) / (2 * WS);
%         b        = (Kg * Ude * CLalfa_rad * rho0)/(2 * WS);
%         c        = 1;
%         eqn      = a*V^2 - b*V - c;
%         Solution = vpasolve(eqn, V);
% 
%         for i = 1:length(Solution)
%            new_VA = cast(Solution(i), 'double');
%            if abs(new_VA) > VA
%                VA = abs(new_VA);
%                nA = nGUST(VA);
%            elseif abs(new_VA) < VA
%                VA = VA; 
%                nA = nA;
%            end
%         end
% 
%         % POINT A
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.VA.value           = Vstall(WS, rho0, CLmax_takeoff, nA);
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.VA.Attributes.unit = "m/s"; 
%         VA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.VA.value;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.nA.value           = nA;
%         nA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.nA.value;
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.nA.Attributes.unit = "g's"; 
% 
%         n_from0toS     = linspace(0.0, nS, numb);
%         V_from0toS     = VS*ones(numb, 1);
% 
%         n_fromStoA     = linspace(nS, nA, numb);
%         V_fromStoA     = Vstall(WS, rho0, CLmax_takeoff, n_fromStoA);
% 
%         n_fromAtoF     = [nA nF];
%         V_fromAtoF     = [VA VF];
% 
%         n_fromFto0     = [nF 0.0];
%         V_fromFto0     = [VF VF];
% 
%         final_envelope = figure;
%         hold on
%         grid on 
%         grid minor
%         ylim([-0.5 nmax+0.5])
%         xlim([0 VF+10])
%         plot(V_gust, n_gust, ':k', 'LineWidth', 0.2)
%         plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth',1)
%         plot(V_fromAtoF, n_fromAtoF, '-b', 'LineWidth',1)
%         plot(V_fromFto0, n_fromFto0, '-b', 'LineWidth',1)
%         plot(V_from0toS, n_from0toS, '-b', 'LineWidth',1)
%         xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
%         ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
%         title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
%         text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
%         text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
%         text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)
%         exportgraphics(flapsgust_envelope,'flaps_final_envelopediagram.pdf','ContentType','vector')
%         exportgraphics(flapsgust_envelope,'flaps_final_envelopediagram.png','ContentType','vector')
% 
%         Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.diagram = final_envelope;
% 
%         % Saving figures inside correct folder
%         fprintf('Saving flaps_final_envelopediagram.pdf in: ');
%         fprintf('\n'); 
%         fprintf('%s\n', SaveFolder);
%         % Moving file inside correct folder
%         movefile flaps_final_envelopediagram.pdf Output
%         movefile flaps_final_envelopediagram.png Output
% 
% end
