%% VERTICAL TAIL SURFACES 
%  Digitalization of the airworthiness rules for the vertical tail loads
%  calculation and structural sizing. 

%% CS - VLA 441 MANOEUVRING LOADS 
%
%  (a) At speeds up to VA, the vertical tail surfaces must be designed to
%      withistand the following conditions.In computing the tail loads, the
%      yawing velocity may be assumed to be zero
%      (1) whit the aeroplane in unaccelerated flight at zero yaw, it is
%          assumed that the rudder control is suddenly displaced to the
%          maximum deflection, as limited by the control stops or by limit
%          pilot forces; 
%      (2) with the rudder deflected as specified in sub-paragraph (a)(1)
%          of this paragraph, it is assumed that the aeroplane yaws to the
%          resulting sideslip angle. In lieu of a rational analysis, an
%          overswing angle equal to 1.3 times the static sideslip angle of
%          sub-paragraph (a)(3) of this paragraph may be assumed; 
%      (3) A yaw angle of 15.0 degrees with the rudder control maintained
%          in the neutral position, except as limited by pilot strength. 
%
%  (b) The average loading of Appendix B, B11 and figure B1 of Appendix B 
%      and the distribution of figures B6, B7, B8 of Appendix B may be used
%      instead of requirements of sub-paragraphs (a)(2), (a)(1) and (a)(3)
%      of this paragraph, respectively.
%
%  (c) The yaw angles specified in sub-paragraph (a)(3) of this paragraph
%      may be reduced if the yaw angle chosen for a particular speed 
%      cannot be exceeded in 
%      (1) steady sideslip conditions; and 
%      (2) uncoordinated rolls from steep banks.

%% AMC VLA 441 MANOEUVRING LOADS 
%  -----------------------------
%  INTERPRETATIVE MATERIAL AND ACCEPTABLE MEANS OF COMPLIANCE 
% 
%  For aeroplanes where the horizontal tail is supported by the vertical
%  tail, the tail surfaces and their supporting structure, including the
%  rear portion of the fuselage, should be designed to whitstand the
%  prescribed loading on the vertica ltail and the roll-moments induced by
%  the horizontal tail acting in the same direction.

%  For T - Tails, in the abscence of a more rational analysis, the rolling
%  moment induced by deflection of the vertical rudder may be computed as
%  follows: 
%                            rho0
%  M_rudder = (0.3) * S_ht * ---- * beta * V^2 * b_ht 
%                              2
%  where 
%
%  M_rudder = induced roll - moment at horizontal tail (N * m) 
%  b_ht     = span of horizonta tail (m) 
%  V        = airspeed (m/s)
%  S_ht     = area of horizontal tail (m^2)
%  rho0     = air densidity at sea level (kg/m^3)
%  beta     = angle of zero - lift line due to rudder deflection; this
%             angle can be defined as follow         
%                         dL
%                       ------- * eta * f_n
%                        d eta 
%             where 
%               dL
%             ------- = change of zero - lift angle of eta * f_ = 1
%              d eta         
%       
%             eta     = rudder deflection 
%             f_n     = effectivity factor in accordance with angle of
%                       rudder deflection

% APPLICABLE LOAD FACTOR 
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.load_factor_n.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value;
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.load_factor_n.Attributes.unit = "g's";
% WING LOADING IN PASCALS
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.value = Aircraft.Certification.Performance.I_Level.Wing_loading_SI.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.load_factor_n.value;
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.Attributes.unit = "Pa"; 
% AREA OF THE VERTICAL TAIL 
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value = 0.10;
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.Attributes.unit = "m^2";
% PRODUCT OF WING LOADING TIMES THE VERTICAL TAIL SURFACES 
%
% P = WS * S * (1/g) --> [P] = [kg]
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Reference_mass.value = (1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.S_vertical_tail.value;
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Reference_mass.Attributes.unit = "kg";
% VERTICAL TAIL SPAN b_vt
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.b_vt.value = 0.4375;
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.b_vt.Attributes.unit = "m";
% CALCULATING THE RATIO OF WING LOADING IN KILOGRAMS PER UNIT SPAN 
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Spanwise_Load_ratio.value = (1/(Aircraft.Constants.g.value))*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Manoeuvring_Wing_Loading.value*Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.b_vt.value;
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Spanwise_Load_ratio.Attributes.unit = "kg/m";
% VERTICAL TAIL ROOT CHORD 
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.value = 0.3136;
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.Attributes.unit = "m";
% QUARTER CHORD DISTANCE: C_ROOT_VT / 4
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_quarter.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.value/4;
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_quarter.Attributes.unit = "m";
% CHORDWISE KILOGRAMS OF LOAD DISTRIBUTION 
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Chordwise_Load_ratio.value = Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Reference_mass.value/Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.c_root_vt.value;
Aircraft.Certification.Regulation.SubpartC.VerticalTailLoads.Chordwise_Load_ratio.Attributes.unit = "kg/m";



%% CS - VLA 443 GUST LOADS 
%  
%  (a) Vertical tail surfaces must be designed to whithstand, in
%      unaccelerated flight at speed V = VC, lateral gusts of the values
%      prescribed for VC in CS - VLA 333 (c).
%
%  (b) In the absence of a more rational analysis, the gust load must be
%      computed as follows 
%
%             K_gt * U_de * V * a_vt * S_vt
%      L_vt = -----------------------------
%                          16.3
%
%      where 
%
%      U_de  = derived gust velocities (m/s);
% 
%      L_vt  = vertical tail loads (daN);
%     
%               (0.88) * mu_gt 
%      K_gt  = --------------- = gust alleviation factor;
%                5.3 + mu_gt
%
%                          2 * M                   K^2
%      mu_gt = ------------------------------- * ------- = lat. mass ratio;
%              rho * c_bar_t * g * a_vt * S_vt   (l_t^2)
%              
%              where 
%                
%              M    = aeroplane mass (kg); 
%              rho  = air density (kg/m^3); 
%              S_vt = area of vertical tail (m^2);
%              l_t  = distance from aeroplane c.g. to lift centre of
%                     vertical surface (m); 
%              a_vt = lift curve slope of vertical tail (1/rad); 
%              V    = aeroplane equivalent speed (m/s); 
%              K    = radius of gyration in yaw (m);
%              g    = acceleration due to gravity (m/s^2).
%
%  (c) The average loading in figure B5 and the distribution in figure B8
%      of Appendix B may be used. 

%% AMC VLA 443 GUST LOADS 
%  -----------------------------
%  INTERPRETATIVE MATERIAL AND ACCEPTABLE MEANS OF COMPLIANCE 
%
%  For aeroplanes where the horizontal tail is supported by the vertical
%  tail, the tail surfaces and their supporting structure including the
%  rear portion of the fuselage should be designed to whithstand the
%  prescribed loadings on the vertical tail and the roll - moments induced
%  by the horizontal tail acting in the same direction. 
%
%  For T - Tails in the abscence of a more rational analysis, the rolling
%  moment induced by gust load may be computed as follow 
% 
%                            rho0
%  M_rudder = (0.3) * S_ht * ---- * V * U * b_ht * K_gf
%    
%  where 
%
%  M_rudder = induced roll - moment at horizontal tail (N * m);
%  K_gf     = gust factor = 1.2; 
%  b_ht     = span of horizontal tail (m); 
%  S_ht     = area of horizontal tail (m^2); 
%  rho0     = density of air at sea level (kg/m^3);
%  V        = airspeed of flight (m/s);
%  U        = gust speed (m/s).

%% CS - VLA 445 OUTBOARD FINS 
%
%  (a) If outboard fins are on the horizonta tail surface, the tail
%      surfaces must be designed for the maximum horizontal surface load in
%      combination with the corresponding loads induced on the vertical
%      surfaces by endplates effects. These induced effects need not be
%      combined with other vertical surface loads. 
%
%  (b) If outboard fins extend above and below the horizontal surface, the 
%      critical vertical surface loading, the load per unit area as
%      determined under CS - VLA 441 and 443, must be applied to 
%      (1) the part of the vertical surfaces above the horizontal surface
%          with 80% of that loading applied to the part below the 
%          horizontal surface; and 
%      (2) the part of the vertical surfaces below the horizontal surface
%          with 80% of that loading applied to the part above the 
%          horizontal surface. 
%
%  (c) The endplate effects of outboard fins must be taken into account in
%      applying the yawing conditions of CS - VLA 441 and 443 to the
%      vertical surfaces in sub - paragraph (b) of this paragraph.