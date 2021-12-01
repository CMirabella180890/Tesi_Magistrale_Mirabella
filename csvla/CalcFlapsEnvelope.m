
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
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.nS.Attributes.unit = "g's"; 
% POINT F
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VF.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.value = Aircraft.Certification.Regulation.SubpartC.Flapsloads.nmax.value;
Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.nF.Attributes.unit = "g's"; 

% VECTOR OF STALL AIRSPEED TO PLOT THE ENVELOPE 
Reg            = Aircraft.Certification.Regulation.value;
Aircraft_name  = Aircraft.Certification.Aircraft_Name.value;
VS             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value;
VF             = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointF.VF.value;
VSpos          = Aircraft.Certification.Regulation.SubpartC.Flapsloads.VSpos_vec.value;
npos           = Aircraft.Certification.Regulation.SubpartC.Flapsloads.n_flaps_vector.value;

Aircraft.Certification.Regulation.SubpartC.Flapsloads.Flaps_envelope.value = flapsenvelope_diagram(obj, npos, nmax, VSpos, VS, VF, Reg, Aircraft_name);

% Saving figures inside correct folder
fprintf('Saving flapsenvelopediagram.pdf in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile flapsenvelopediagram.pdf Output
movefile flapsenvelopediagram.png Output



