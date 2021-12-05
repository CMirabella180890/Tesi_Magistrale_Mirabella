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

%% IMPLEMENTATION 
numb = 1e3;
nmax = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;

% CALCULATION OF THE LOAD FACTORS VECTOR 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value = calcnflap(obj, nmax);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.Attributes.unit = "g's";

% CALCULATION OF VS - CLEAN STALL SPEED
n     = 1.0; 
rho   = Aircraft.Certification.ISA_Condition.rho0.value;
WS    = Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value;
CLmax = Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value;
VS    = calcvs(obj, rho, WS, CLmax, n);

% CALCULATION OF VS1 - FLAPS DEPLOYED STALL SPEED
CLmax = Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_flaps.value;
VS1   = calcvs(obj, rho, WS, CLmax, n);

% EVALUATION OF VF - FLAPS DEPLOYED AIRSPEED 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.VF.value = calcnVF(obj, VS, VS1);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.VF.Attributes.unit = "m/s";

% STALLING SPEED VECTOR
Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value = calcvs(obj, rho, WS, CLmax, Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.Attributes.unit = "m/s";

% EVALUATION OF STALL AND FLAP MANOEUVRING POINT 
% POINT S 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value = VS1;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.value = 1.0;
nS = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.Attributes.unit = "g's"; 
% POINT A
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointA.VA.value = Vstall(WS, rho, CLmax, nmax);
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
VF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value;
npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value;

n_from0toS     = linspace(0.0, nS, numb);
V_from0toS     = VS*ones(numb, 1);

n_fromStoA     = linspace(nS, nA, numb);
V_fromStoA     = Vstall(WS, rho, CLmax, n_fromStoA);

n_fromAtoF     = nmax*ones(numb, 1);
V_fromAtoF     = linspace(VA, VF, numb);

n_fromFto0     = linspace(nF, 0.0, numb);
V_fromFto0     = VF*ones(numb,1);

flaps_envelope = figure;
hold on
grid on 
grid minor
ylim([-0.5 nmax+0.5])
xlim([0 VF+10])
plot(VSpos, npos, ':r', 'LineWidth',0.2)
plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth',1)
plot(V_fromAtoF, n_fromAtoF, '-b', 'LineWidth',1)
plot(V_fromFto0, n_fromFto0, '-b', 'LineWidth',1)
plot(V_from0toS, n_from0toS, '-b', 'LineWidth',1)
xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
% text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)
exportgraphics(flaps_envelope,'flapsenvelopediagram.pdf','ContentType','vector')
exportgraphics(flaps_envelope,'flapsenvelopediagram.png','ContentType','vector')

Aircraft.Certification.Regulation.SubpartC.Flapsloads.Flaps_envelope.value = flaps_envelope;

% Aircraft.Certification.Regulation.SubpartC.Flapsloads.Flaps_envelope.value = flapsenvelope_diagram(obj, npos, nmax, VSpos, VS, VF, Reg, Aircraft_name);

% Saving figures inside correct folder
fprintf('Saving flapsenvelopediagram.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile flapsenvelopediagram.pdf Output
movefile flapsenvelopediagram.png Output

%% FLAPS DEPLOYED GUST ENVELOPE 
VF  = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
WS  = Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value;
rho = Aircraft.Certification.ISA_Condition.rho0.value;
MAC = Aircraft.Geometry.Wing.mac.value; 
a   = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
g   = Aircraft.Constants.g.value;
Ude = 7.62; % Gust magnitude

% CALCULATION OF THE MASS FACTOR
mu_g = calcmug(obj, WS, MAC, a, rho, g); 

% GUST ALLEVIATION FACTOR 
Kg   = calckg(obj, mu_g);

% CALCULATION OF THE GUST LOAD FACTOR AT V = VF 
nGUST  = @(V) 1.0 + V*((0.5*rho*a*Kg*Ude)/(WS));
V_gust = linspace(0.0, VF, numb); 
n_gust = nGUST(V_gust);

% STORE VALUES 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.mu_g.value = mu_g;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.mu_g.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.Kg.value = Kg;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.Kg.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.Vgust.value = V_gust;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.ngust.value = n_gust;

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
plot(V_gust, n_gust, ':k', 'LineWidth', 0.2)
plot(VSpos, npos, ':r', 'LineWidth', 0.2)
plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth',1)
plot(V_fromAtoF, n_fromAtoF, '-b', 'LineWidth',1)
plot(V_fromFto0, n_fromFto0, '-b', 'LineWidth',1)
plot(V_from0toS, n_from0toS, '-b', 'LineWidth',1)
xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)
exportgraphics(flapsgust_envelope,'flaps_gust_envelopediagram.pdf','ContentType','vector')
exportgraphics(flapsgust_envelope,'flaps_gust_envelopediagram.png','ContentType','vector')

Aircraft.Certification.Regulation.SubpartC.Flapsloads.Gust_envelope.diagram = flapsgust_envelope;
% Saving figures inside correct folder
fprintf('Saving flapsenvelopediagram.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile flaps_gust_envelopediagram.pdf Output
movefile flaps_gust_envelopediagram.png Output

%% FINAL ENVELOPE WITH FLAPS DEPLOYED 
nmax           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
CLalfa         = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value;
VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value;
VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value;
V_g            = linspace(VS, VF, numb); 
n_g            = nGUST(V_g);
nS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.value;
VF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
nF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.value;
Reg            = Aircraft.Certification.Regulation.value;
Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;

syms a b c V 
a        = (rho*CLmax)/(2*WS);
b        = (Kg*Ude*CLalfa*rho)/(2*WS);
c        = 1;
eqn      = a*V^2 - b*V - c;
Solution = vpasolve(eqn, V);

for i = 1:length(Solution)
   new_VA = cast(Solution(i), 'double');
   if abs(new_VA) > VA
       VA = abs(new_VA);
       nA = nGUST(VA);
   elseif abs(new_VA) < VA
       VA = VA; 
       nA = nA;
   end
end

% POINT A
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.VA.value           = Vstall(WS, rho, CLmax, nA);
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.VA.Attributes.unit = "m/s"; 
VA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.VA.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.nA.value           = nA;
nA = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.nA.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.Final_envelope.PointA.nA.Attributes.unit = "g's"; 

n_from0toS     = linspace(0.0, nS, numb);
V_from0toS     = VS*ones(numb, 1);

n_fromStoA     = linspace(nS, nA, numb);
V_fromStoA     = Vstall(WS, rho, CLmax, n_fromStoA);

n_fromAtoF     = [nA nF];
V_fromAtoF     = [VA VF];

n_fromFto0     = [nF 0.0];
V_fromFto0     = [VF VF];

final_envelope = figure;
hold on
grid on 
grid minor
ylim([-0.5 nmax+0.5])
xlim([0 VF+10])
plot(V_gust, n_gust, ':k', 'LineWidth', 0.2)
plot(V_fromStoA, n_fromStoA, '-r', 'LineWidth',1)
plot(V_fromAtoF, n_fromAtoF, '-b', 'LineWidth',1)
plot(V_fromFto0, n_fromFto0, '-b', 'LineWidth',1)
plot(V_from0toS, n_from0toS, '-b', 'LineWidth',1)
xlabel("Airspeed - $V$ (m/s)", "Interpreter", "latex")
ylabel("Load factor - $n$ (g's)", "Interpreter", "latex")
title("Flaps envelope diagram per ", Reg, "Interpreter", "latex") % Applied regulation from 'Aircraft' struct
text(15, 1.8, Aircraft_name)                                % Aircraft name inside the plot
text(VS, nS, '\fontname{Courier} S', 'FontSize', 6)
text(VF, nmax, '\fontname{Courier} F', 'FontSize', 6)
exportgraphics(flapsgust_envelope,'flaps_final_envelopediagram.pdf','ContentType','vector')
exportgraphics(flapsgust_envelope,'flaps_final_envelopediagram.png','ContentType','vector')

Aircraft.Certification.Regulation.SubpartC.Flapsloads.final_envelope.diagram = final_envelope;

% Saving figures inside correct folder
fprintf('Saving flaps_final_envelopediagram.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile flaps_final_envelopediagram.pdf Output
movefile flaps_final_envelopediagram.png Output
