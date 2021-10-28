%% Script to evaluate shear and bending moment distr. along the main wing span

%   DESCRIPTION
%    In this script, the lift and drag distribution along the span will be 
%    used to evaluate shear and bending moment along the main wing span.
%    First, it is necessary to store the critical point, which can be
%    exctracted from the V-N diagram. 

%% STORE CRITICAL POINTS OF THE FINAL ENVELOPE
% It is crucial to store inside the struct variable 'Aircraft' all the
% critical points coming from the Final envelope diagram; this points are
% associated with the most critical, in-flight structural loads for the 
% main wing and horizontal tail structures. These loads are used to size
% structural elements (i.e. main wing spar, main wing ribs...) of the
% aircraft. The points are (V, n) couples.

% For a more in depth description, look ahead in this file
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value = cl_unit_lift(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CL.value(3), ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CL.value(4), ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value(3,:)', ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value(4,:)');
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.Attributes.unit = "Non dimensional";  

% Half-span stations vector 
half_span = linspace(0, ...
    Aircraft.Geometry.Wing.b.value*0.5, ...
    length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value));

% Point S of the final envelope
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.y_half_span.value = half_span;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.y_half_span.Attributes.unit = "m";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value = 'Point S';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value(end);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.value = 0.5*Aircraft.Certification.ISA_Condition.rho0.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.value^2);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.Attributes.unit = "N/m^2";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.nS.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value(1); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.nS.Attributes.unit = "g's"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_posstall_new.value(1); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.Attributes.unit = "Non dimensional"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CD_S.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value(1);  
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CD_S.Attributes.unit = "Non dimensional"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alpha_S.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_posstall_new.value(1); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alpha_S.Attributes.unit = "Degrees";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_posstall_new.value(1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LS.Attributes.unit = "daN"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTail_S.value = (0.5)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.value^2)* ...
        (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.rho.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_positivestall.value(1))*(1E-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.LHTail_S.Attributes.unit = "daN";

% Point A of the final envelope
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.y_half_span.value = half_span;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.y_half_span.Attributes.unit = "m";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.point_name.value = 'Point A';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value = 0.5*Aircraft.Certification.ISA_Condition.rho0.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value^2);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.Attributes.unit = "N/m^2";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.Attributes.unit = "g's"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_posstall_new.value(1); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.Attributes.unit = "Non dimensional"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CD_A.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value(1);  
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CD_A.Attributes.unit = "Non dimensional"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alpha_A.value = alpha_calc(obj1, ...
                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alpha_A.Attributes.unit = "Degrees";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA.value = (0.5)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.value^2)* ...
        (Aircraft.Certification.ISA_Condition.rho.value)*(Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_posstall_new.value(1))*(1E-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LA.Attributes.unit = "daN"; 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTail_A.value = (0.5)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Design_manoeuvring_speed_VA.value^2)* ...
%         (Aircraft.Geometry.Wing.S.value)*(Aircraft.Certification.ISA_Condition.rho.value)*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_positivestall.value(1))*(1E-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTail_A.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_positivestall.value(end);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.LHTail_A.Attributes.unit = "daN";

% Point C of the final envelope
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.y_half_span.value = half_span;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.y_half_span.Attributes.unit = "m";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value = 'Point C';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_speed.value(end);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value = 0.5*Aircraft.Certification.ISA_Condition.rho0.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value^2);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.Attributes.unit = "N/m^2"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Positive_stall_load_factor.value(end); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.Attributes.unit = "g's"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_posstall_new.value(end); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.Attributes.unit = "Non dimensional"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CD_C.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_positivestall.value(end);  
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CD_C.Attributes.unit = "Non dimensional"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alpha_C.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_posstall_new.value(end); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alpha_C.Attributes.unit = "Degrees";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_posstall_new.value(end);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LC.Attributes.unit = "daN"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTail_C.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_positivestall.value(end);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.LHTail_C.Attributes.unit = "daN";

% Point D of the final envelope 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.y_half_span.value = half_span;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.y_half_span.Attributes.unit = "m";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value = 'Point D';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_dive.value(end);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value = 0.5*Aircraft.Certification.ISA_Condition.rho0.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value^2);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.Attributes.unit = "N/m^2";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive.value(end);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.Attributes.unit = "g";  
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromDtoE_new.value(1); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_D.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromDtoE.value(1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_D.Attributes.unit = "Non dimensional";  
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alpha_D.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_new.value(1); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alpha_D.Attributes.unit = "Degrees";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromDtoE_new.value(1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LD.Attributes.unit = "daN"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTail_D.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromDtoE.value(1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.LHTail_D.Attributes.unit = "daN";

% Point F of the final envelope
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.y_half_span.value = half_span;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.y_half_span.Attributes.unit = "m";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value = 'Point F';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Airspeed_dive.value(end);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value = 0.5*Aircraft.Certification.ISA_Condition.rho0.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.value^2);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.Attributes.unit = "N/m^2";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_dive.value(end);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.Attributes.unit = "g"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value(1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value(1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.Attributes.unit = "Non dimensional";  
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alpha_F.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromFtoE_new.value(1); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alpha_F.Attributes.unit = "Radians"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE_new.value(1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LF.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTail_F.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromFtoE.value(1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.LHTail_F.Attributes.unit = "daN";

% Point G of the final envelope
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.y_half_span.value = half_span;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.y_half_span.Attributes.unit = "m";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.point_name.value = 'Point G';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Negative_Design_manoeuvring_speed_VG.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG.Attributes.unit = "m/s";   
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value = 0.5*Aircraft.Certification.ISA_Condition.rho0.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG.value^2);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.Attributes.unit = "N/m^2";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.nG.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.nG.Attributes.unit = "g";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value = -(1/Aircraft.Geometry.Wing.S.value)*((Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value)/Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.value = Aircraft.Certification.Aerodynamic_data.CD0.value + (Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value^2)/(pi*Aircraft.Certification.Aerodynamic_data.e.value*Aircraft.Geometry.Wing.AR.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.Attributes.unit = "Non dimensional";  
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alpha_G.value = ((Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value - Aircraft.Certification.Aerodynamic_data.CL0.value)/(Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value))-Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alpha_G.Attributes.unit = "degrees"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value*Aircraft.Geometry.Wing.S.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value*(1e-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LG.Attributes.unit = "daN";   
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTail_G.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_Tail_negativestall.value(end)*Aircraft.Geometry.Wing.S.value*(1e-1); % Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromFtoE.value(end);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.LHTail_G.Attributes.unit = "daN"; 

% Point E of the final envelope
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.y_half_span.value = half_span;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.y_half_span.Attributes.unit = "m";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.point_name.value = 'Point E';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.VE.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_speed.value(end);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.VE.Attributes.unit = "m/s"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value = 0.5*Aircraft.Certification.ISA_Condition.rho0.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.VE.value^2);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.Attributes.unit = "N/m^2";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.nE.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Negative_stall_load_factor.value(end);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.nE.Attributes.unit = "g"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CL_fromFtoE_new.value(end);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CD_E.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.CD_fromFtoE.value(end);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CD_E.Attributes.unit = "Non dimensional";  
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alpha_E.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.alpha_fromDtoE_new.value(end); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alpha_E.Attributes.unit = "Degrees"; 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.WBLift_fromFtoE_new.value(end);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LE.Attributes.unit = "daN";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LHTail_E.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.HTail_Lift_fromFtoE.value(end);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.LTHail_E.Attributes.unit = "daN"; 

%% LIFT CURVE AND LIFT COEFFICIENT EVALUATED AT THE FINAL ENVELOPE POINTS

Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Diagram_lift_coefficient_comparison.value = Lift_coefficients_Points(Aircraft.Certification.Aerodynamic_data.AOA_aux.value, ...
                           Aircraft.Certification.Aerodynamic_data.CL_Full_model.value, ...
                           Aircraft.Certification.Aerodynamic_data.CL_Full_model_invertedflight.value, ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alpha_S.value, ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alpha_A.value, ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alpha_C.value, ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alpha_D.value + Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value, ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alpha_F.value + Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value, ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alpha_G.value - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value, ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alpha_E.value + Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value, ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value, ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value, ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value, ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value, ...
                           abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value), ...
                           abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value), ...
                           abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.value), ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value, ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.point_name.value, ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value, ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value, ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value, ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.point_name.value, ...
                           Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.point_name.value);
pause(1);
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Diagram_lift_coefficient_comparison.value, 'LiftComparison.pdf', 'ContentType', 'vector')
% Saving figures inside correct folder
fprintf('Saving in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile LiftComparison.pdf Output

%% LOAD A CLASS OF FUNCTIONS USEFUL TO EVALUATE ALL THE REQUIRED DATA
obj2 = ShearBendingTorsion; 
pause(5/1000);
% Test speed
% To check the work, a test case issue is provided. OpenVSP is called with
% the following value of the airspeed. Also, the chosen flight condition is
% relative to a lift coefficient equal to one (CL = 1). All the geometrical
% parameters are inserted inside the OpenVSP input file. 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Test_speed.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VA.value; 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Test_speed.Attributes.unit = "m/s";

%% DEFINE A CHORD DISTRIBUTION 
% Lift coefficient distribution at a global CL equal to one
%
% cl_at_CL1 = cl_unit_lift(CL1, CL2, cl1, cl2)
% A simple function to evaluate the lift coefficient distribution along the
% span cl = cl(y) when the associated global lift coefficient of the whole
% wing is equal to 1.0; the function use a method similar to that suggested
% by Abbott in Theory of Wing Section. See the complete documentation
% inside the cl_unit_lift.m file
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value = cl_unit_lift(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CL.value(3), ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CL.value(4), ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value(3,:)', ...
                                                                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value(4,:)');
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.Attributes.unit = "Non dimensional";        

% A vector which contains all the stations along the main wing span.
Aircraft.Geometry.Wing.y.value = linspace(0, ...
                                          Aircraft.Geometry.Wing.b.value*0.5, ...
                                          length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value));
Aircraft.Geometry.Wing.y.Attributes.unit = "m";
                                    
% Main wing taper ratio
Aircraft.Geometry.Wing.taper.value = Aircraft.Geometry.Wing.ctip.value/Aircraft.Geometry.Wing.croot.value;
Aircraft.Geometry.Wing.taper.Attributes.unit = "Non dimensional";

% Calculation of a chord distribution with a convenient, simple function.
% 
% c(y) = calc_chord(Swing, taper_ratio, span, y)
% A complete documentation of this function is included inside the class
% ShearBendingTorsion.m
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value = calc_chord(obj2, Aircraft.Geometry.Wing.S.value, ...
                                                                                         Aircraft.Geometry.Wing.taper.value, ...
                                                                                         Aircraft.Geometry.Wing.b.value, ...
                                                                                         Aircraft.Geometry.Wing.y.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.Attributes.unit = "m";

%% POINT S CALCULATIONS                 
% Lift coefficient distribution along the span at the Point S
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cl_S.value = abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.CL_S.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cl_S.Attributes.unit = "Non dimensional";

% Drag coefficient ditribution along the span at the Point S (close to
% stall)
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cd_S.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:)';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cd_S.Attributes.unit = "Non dimensional";

% Pitching moment coefficient distribution along the span at Point S
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cm_S.value =  Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value(7,:)';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cm_S.Attributes.unit = "Non dimensional";

%% PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
% In this section of the code two vectors are defined to store the product 
% c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCl_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cl_S.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCl_distr.Attributes.unit = "m";

Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCd_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cd_S.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCd_distr.Attributes.unit = "m";

%% ALPHA VECTOR ALONG THE MAIN WING SPAN 
cd .. 
cd utilities
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

% Calling the function alpha_calc
% Angle of attack = ALPHA(obj, CL, CL0, CL_star, CLalpha, a, b, c)
% This function is included in the class file aero_model.m; a complete
% description of this function is included inside that file.
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_alongspan_deg.value = alpha_calc(obj1, ...
                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cl_S.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p);
%                          alpha_calc_lin(obj1, ...
%                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cl_S.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);

Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_alongspan_deg.Attributes.unit = "degrees"; 

% Convert to radians 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_alongspan_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_alongspan_deg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_alongspan_rad.Attributes.unit = "radians";

% Changing folder 
cd .. 
cd csvla
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

%% TOTAL ANGLE OF ATTACK DATA STORAGE

% AoA_Tot = AoA + Twist_angle of the main wing
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_deg.value = Aircraft.Geometry.Wing.twist_angle.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.alpha_S.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_deg.Attributes.unit = "Degrees"; 

% Convert in radians 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_deg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_rad.Attributes.unit = "Radians";

%% PROJECTION OF FORCES ALONG WING AXES - NORMAL AND AXIAL FORCES 

% Calculation of the normal force coefficient
% N = calc_normal_force(AoA_Tot, cCl, cCd)
% This function will be used to evaluate the normal force coefficients
% distribution along the span; it is possible to fin a complete
% documentation inside the class file ShearBendingTorsion.m 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCz.value = calc_normal_force(obj2, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_deg.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCl_distr.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCd_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCz.Attributes.unit = "m"; 

% Calculation of the axial force coefficient 
% A = calc_normal_force(AoA_Tot, cCl, cCd)
% This function will be used to evaluate the axial force coefficients
% distribution along the span; it is possible to fin a complete
% documentation inside the class file ShearBendingTorsion.m 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCa.value = calc_axial_force(obj2, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.AoA_Tot_deg.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCl_distr.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCd_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCa.Attributes.unit = "m"; 

% Normal force 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Normal_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCz.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Normal_force.Attributes.unit = "N/m";

% Axial force 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Axial_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cCa.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Axial_force.Attributes.unit = "N/m";

%% SHEAR FORCE CALCULATION 
% A = calc_shear_force(AoA_Tot, y, cCZ)
% A complete description of this function is available inside the class
% file ShearBendingTorsion.m 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_distr.value = calc_shear_force(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Normal_force.value)*(1e-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_distr.Attributes.unit = "daN";

%% BENDING MOMENT CALCULATION 
% BM = calc_bend_mom(y, S)
% A complete description of this function is included inside the class file
% ShearBendingTorsion.m
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Bend_mom_distr.value = calc_bend_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Bend_mom_distr.Attributes.unit = "daN*m";

%% PLANS FOR THE STRUCTURAL DIMENSIONING
% To correctly size the aerostructures of the main lifting surface 
% it is necessary to apply the procedure just developed to the 
% critical points coming from the V-N diagram. Those point represents 
% the most demanding flight conditions that our aircraft could survive. 
% Those points are all stored inside: 
% --> Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope
% Retrieve the values and apply formulas to them. 

% Pitching moment per unit length
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.m_distr.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cm_S.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cm_S.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.m_distr.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.cm_S.value(i)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qS.value*((Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value(i))^2);
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.m_distr.Attributes.unit = "N";

% Torque applied
% T = calc_tors_mom(obj, y, m)
% A complete distribution of this function is included inside the class
% file ShearBendingTorsion.m
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Tors_mom_distr.value = calc_tors_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.m_distr.value)*(1e-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Tors_mom_distr.Attributes.unit = "daN*m";

% Subplots with Shear, Bending moment and Torsion
%  fig1 = Bending_Shear_diag(y, Shear, Bend_mom, Torsion, Point)
% A complete description of this function is included inside the class file
% ShearBendingTorsion.m
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_BendMom_diagram.value = Shear_Bending_Torsion_diag(obj2, flip(Aircraft.Geometry.Wing.y.value), ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Bend_mom_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Tors_mom_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value);
pause(1);
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointS.pdf', 'ContentType', 'vector')
% Saving figures inside correct folder
fprintf('Saving in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile ShearBendingTorsionDiagramPointS.pdf Output

%% POINT A CALCULATIONS                 
% Lift coefficient distribution along the span at the Point A
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cl_A.value = abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.CL_A.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cl_A.Attributes.unit = "Non dimensional";

% Drag coefficient ditribution along the span at the Point A (close to
% stall)
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cd_A.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:)';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cd_A.Attributes.unit = "Non dimensional";

% Pitching moment coefficient distribution along the span at Point A
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cm_A.value =  Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value(7,:)';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cm_A.Attributes.unit = "Non dimensional";

%% PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
% In this section of the code two vectors are defined to store the product 
% c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCl_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cl_A.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCl_distr.Attributes.unit = "m";

Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCd_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cd_A.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCd_distr.Attributes.unit = "m";

%% ALPHA VECTOR ALONG THE MAIN WING SPAN 
cd .. 
cd utilities
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

% Calling the function alpha_calc
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_alongspan_deg.value = alpha_calc(obj1, ...
                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cl_A.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p);
%                          alpha_calc_lin(obj1, ...
%                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cl_A.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_alongspan_deg.Attributes.unit = "degrees"; 

% Convert to radians 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_alongspan_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_alongspan_deg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_alongspan_rad.Attributes.unit = "radians";

% Changing folder 
cd .. 
cd csvla
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

%% TOTAL ANGLE OF ATTACK DATA STORAGE

% AoA_Tot = AoA + Twist_angle of the main wing
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_deg.value = Aircraft.Geometry.Wing.twist_angle.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.alpha_A.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_deg.Attributes.unit = "Degrees"; 

% Convert in radians 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_deg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_rad.Attributes.unit = "Radians";

%% PROJECTION OF FORCES ALONG WING AXES - NORMAL AND AXIAL FORCES 

% Calculation of the normal force coefficient
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCz.value = calc_normal_force(obj2, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_deg.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCl_distr.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCd_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCz.Attributes.unit = "m"; 

% Calculation of the axial force coefficient 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCa.value = calc_axial_force(obj2, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.AoA_Tot_deg.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCl_distr.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCd_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCa.Attributes.unit = "m"; 

% Normal force 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Normal_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCz.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Normal_force.Attributes.unit = "N/m";

% Axial force 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Axial_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cCa.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Axial_force.Attributes.unit = "N/m";

%% SHEAR FORCE CALCULATION 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_distr.value = calc_shear_force(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Normal_force.value)*(1e-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_distr.Attributes.unit = "daN";

%% BENDING MOMENT CALCULATION 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Bend_mom_distr.value = calc_bend_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Bend_mom_distr.Attributes.unit = "daN*m";

%% SHEAR, BENDING MOMENT AND TORSION DIAGRAM

% Pitching moment per unit length
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.m_distr.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cm_A.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cm_A.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.m_distr.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cm_A.value(i)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value*((Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value(i))^2);
end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.m_distr.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.cm_A.value.*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.qA.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value.^2);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.m_distr.Attributes.unit = "N";

Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Tors_mom_distr.value = calc_tors_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.m_distr.value)*(1e-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Tors_mom_distr.Attributes.unit = "daN*m";

% Subplots with Shear, Bending moment and Torsion
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_BendMom_diagram.value = Shear_Bending_Torsion_diag(obj2, flip(Aircraft.Geometry.Wing.y.value), ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Bend_mom_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Tors_mom_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.point_name.value);
pause(1);
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointA.pdf', 'ContentType', 'vector')
% Saving figures inside correct folder
fprintf('Saving in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile ShearBendingTorsionDiagramPointA.pdf Output

%% POINT C CALCULATIONS                 
% Lift coefficient distribution along the span at the Point C
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cl_C.value = abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.CL_C.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cl_C.Attributes.unit = "Non dimensional";

% Drag coefficient ditribution along the span at the Point C (close to
% stall)
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cd_C.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:)';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cd_C.Attributes.unit = "Non dimensional";

% Pitching moment coefficient distribution along the span at Point C
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.value =  Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value(7,:)';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.Attributes.unit = "Non dimensional";

%% IMPORT DATA FROM THE CORRECT SOURCE

%% PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
% In this section of the code two vectors are defined to store the product 
% c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCl_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cl_C.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCl_distr.Attributes.unit = "m";

Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCd_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cd_C.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCd_distr.Attributes.unit = "m";

%% ALPHA VECTOR ALONG THE MAIN WING SPAN 
cd .. 
cd utilities
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

% Calling the function alpha_calc
% Angle of attack = ALPHA(obj, CL, CL0, CL_star, CLalpha, a, b, c)
% 
%   This function is able to evaluate the angle of attack of
%   the aircraft in the prescribed flight condition. It must be
%   noticed that. Search the complete description inside the class
%   file aero_model.m
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_alongspan_deg.value = alpha_calc(obj1, ...
                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cl_C.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p);
%                          alpha_calc_lin(obj1, ...
%                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cl_C.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_alongspan_deg.Attributes.unit = "degrees"; 

% Convert to radians 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_alongspan_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_alongspan_deg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_alongspan_rad.Attributes.unit = "radians";

% Changing folder 
cd .. 
cd csvla
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

%% TOTAL ANGLE OF ATTACK DATA STORAGE

% AoA_Tot = AoA + Twist_angle of the main wing
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_deg.value = Aircraft.Geometry.Wing.twist_angle.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.alpha_C.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_deg.Attributes.unit = "Degrees"; 

% Convert in radians 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_deg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_rad.Attributes.unit = "Radians";

%% PROJECTION OF FORCES ALONG WING AXES - NORMAL AND AXIAL FORCES 

% Calculation of the normal force coefficient
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCz.value = calc_normal_force(obj2, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_deg.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCl_distr.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCd_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCz.Attributes.unit = "m"; 

% Calculation of the axial force coefficient 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCa.value = calc_axial_force(obj2, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.AoA_Tot_deg.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCl_distr.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCd_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCa.Attributes.unit = "m"; 

% Normal force 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Normal_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCz.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Normal_force.Attributes.unit = "N/m";

% Axial force 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Axial_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cCa.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Axial_force.Attributes.unit = "N/m";

%% SHEAR FORCE CALCULATION 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_distr.value = calc_shear_force(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Normal_force.value)*(1e-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_distr.Attributes.unit = "daN";

%% BENDING MOMENT CALCULATION 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Bend_mom_distr.value = calc_bend_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Bend_mom_distr.Attributes.unit = "daN*m";

%% SHEAR, BENDING MOMENT AND TORSION DIAGRAM

% Pitching moment per unit length
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.m_distr.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.m_distr.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.value(i)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value*((Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value(i))^2);
end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.m_distr.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.cm_C.value.*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.qC.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value.^2);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.m_distr.Attributes.unit = "N";

Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Tors_mom_distr.value = calc_tors_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.m_distr.value)*(1e-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Tors_mom_distr.Attributes.unit = "daN*m";

% Subplots with shear, bending moment and torsion diagram
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_BendMom_diagram.value = Shear_Bending_Torsion_diag(obj2, flip(Aircraft.Geometry.Wing.y.value), ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Bend_mom_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Tors_mom_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value);
pause(1);
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointC.pdf', 'ContentType', 'vector')
% Saving figures inside correct folder
fprintf('Saving in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile ShearBendingTorsionDiagramPointC.pdf Output

%% POINT D CALCULATIONS                 
% Lift coefficient distribution along the span at the Point D
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cl_D.value = abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CL_D.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cl_D.Attributes.unit = "Non dimensional";

% Interpolation to obtain Cd at the desired global CL
x   = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:));
y   = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(:,1));
xi  = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(7,:));
yi  = 1:0.1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value(:,1));
[XX,YY] = meshgrid(x,y);
[XI,YI] = meshgrid(xi,yi);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cd.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value, ...
                                                                        XI, YI, 'spline');
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cd.Attributes.unit = 'Non dimensional';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_Interpolation_Graph.value = cd_interpolation_graph(x, y, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value, ... 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cd.value, XI, YI);            
            
% Pitching moment interpolation along the span
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cm.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value, ...
                                                                        XI, YI, 'spline'); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cm.Attributes.unit = 'Non dimensional';
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CM_Interpolation_Graph.value = cm_interpolation_graph(x, y, ...
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value, ... 
                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cm.value, XI, YI);
            
% Export the interpolation carpet plot
% Cd interpolation and saving diagram
pause(1);
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_Interpolation_Graph.value, 'CdInterpolationdiag.pdf', 'ContentType', 'vector')
% Saving figures inside correct folder
fprintf('Saving in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile CdInterpolationdiag.pdf Output

% Cd interpolation and saving diagram
pause(1);
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CM_Interpolation_Graph.value, 'CmInterpolationdiag.pdf', 'ContentType', 'vector')
% Saving figures inside correct folder
fprintf('Saving in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile CmInterpolationdiag.pdf Output
            
% Integrate along semi-span and assign drag and pitching moment along span
Aircraft.Geometry.Wing.half_span_y.value = linspace(0, ...
                                          Aircraft.Geometry.Wing.b.value*0.5, ...
                                          length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value));
Aircraft.Geometry.Wing.half_span_y.Attributes.unit = 'm';

% Selection of the interpolated distribution of CD and CM
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Global_CD.value = zeros(length(yi), 1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Global_CD.Attributes.unit = "Non dimensional";
for i = 1:length(yi)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Global_CD.value(i) = trapz(Aircraft.Geometry.Wing.half_span_y.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cd.value(i,:));
    if abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Global_CD.value(i) - Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.CD_D.value) < 1e-3
       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cd_D.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cd.value(i,:)';
       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Interpolated_Cm.value(i,:)';
    end
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cd_D.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.Attributes.unit = "Non dimensional";

%% IMPORT DATA FROM THE CORRECT SOURCE

%% PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
% In this section of the code two vectors are defined to store the product 
% c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCl_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cl_D.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCl_distr.Attributes.unit = "m";

Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCd_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cd_D.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCd_distr.Attributes.unit = "m";

%% ALPHA VECTOR ALONG THE MAIN WING SPAN 
cd .. 
cd utilities
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

% Calling the function alpha_calc
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_alongspan_deg.value = alpha_calc(obj1, ...
                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cl_D.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p);
%                          alpha_calc_lin(obj1, ...
%                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cl_D.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_alongspan_deg.Attributes.unit = "degrees"; 

% Convert to radians 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_alongspan_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_alongspan_deg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_alongspan_rad.Attributes.unit = "radians";

% Changing folder 
cd .. 
cd csvla
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

%% TOTAL ANGLE OF ATTACK DATA STORAGE

% AoA_Tot = AoA + Twist_angle of the main wing
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_deg.value = Aircraft.Geometry.Wing.twist_angle.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.alpha_D.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_deg.Attributes.unit = "Degrees"; 

% Convert in radians 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_deg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_rad.Attributes.unit = "Radians";

%% PROJECTION OF FORCES ALONG WING AXES - NORMAL AND AXIAL FORCES 

% Calculation of the normal force coefficient
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCz.value = calc_normal_force(obj2, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_deg.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCl_distr.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCd_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCz.Attributes.unit = "m"; 

% Calculation of the axial force coefficient 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCa.value = calc_axial_force(obj2, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.AoA_Tot_deg.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCl_distr.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCd_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCa.Attributes.unit = "m"; 

% Normal force 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Normal_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCz.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Normal_force.Attributes.unit = "N/m";

% Axial force 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Axial_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cCa.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Axial_force.Attributes.unit = "N/m";

%% SHEAR FORCE CALCULATION 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_distr.value = calc_shear_force(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Normal_force.value)*(1e-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_distr.Attributes.unit = "daN";

%% BENDING MOMENT CALCULATION 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Bend_mom_distr.value = calc_bend_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Bend_mom_distr.Attributes.unit = "daN*m";

%% SHEAR, BENDING MOMENT AND TORSION DIAGRAM

% Pitching moment per unit length
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.m_distr.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.m_distr.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.value(i)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value*((Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value(i))^2);
end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.m_distr.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.cm_D.value.*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.qD.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value.^2);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.m_distr.Attributes.unit = "N";

Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Tors_mom_distr.value = calc_tors_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.m_distr.value)*(1e-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Tors_mom_distr.Attributes.unit = "daN*m";

% Subplots with shear, bending moment and torsion 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_BendMom_diagram.value = Shear_Bending_Torsion_diag(obj2, flip(Aircraft.Geometry.Wing.y.value), ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Bend_mom_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Tors_mom_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value);
pause(1);
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointD.pdf', 'ContentType', 'vector')
% Saving figures inside correct folder
fprintf('Saving in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile ShearBendingTorsionDiagramPointD.pdf Output

%% POINT F CALCULATIONS                 
% Lift coefficient distribution along the span at the Point F
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cl_F.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CL_F.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cl_F.Attributes.unit = "Non dimensional";

% Interpolation to obtain Cd at the desired global CL
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Cd.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value, ...
                                                                        XI, YI, 'spline');
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Cd.Attributes.unit = 'Non dimensional';

% Pitching moment interpolation along the span
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Cm.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value, ...
                                                                        XI, YI, 'spline'); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Cm.Attributes.unit = 'Non dimensional';
          
% Integrate along semi-span and assign drag and pitching moment along span
Aircraft.Geometry.Wing.half_span_y.value = linspace(0, ...
                                          Aircraft.Geometry.Wing.b.value*0.5, ...
                                          length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value));
Aircraft.Geometry.Wing.half_span_y.Attributes.unit = 'm';

% Selection of CD and CM distribution for Point F
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Global_CD.value = zeros(length(yi), 1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Global_CD.Attributes.unit = "Non dimensional";
for i = 1:length(yi)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Global_CD.value(i) = trapz(Aircraft.Geometry.Wing.half_span_y.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Cd.value(i,:));
    if abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Global_CD.value(i) - Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.CD_F.value) < 1e-3
       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cd_F.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Cd.value(i,:)';
       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cm_F.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Interpolated_Cm.value(i,:)';
    end
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cd_F.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cm_F.Attributes.unit = "Non dimensional";

%% IMPORT DATA FROM THE CORRECT SOURCE

%% PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
% In this section of the code two vectors are defined to store the product 
% c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCl_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cl_F.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCl_distr.Attributes.unit = "m";

Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCd_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cd_F.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCd_distr.Attributes.unit = "m";

%% ALPHA VECTOR ALONG THE MAIN WING SPAN 
cd .. 
cd utilities
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

% Calling the function alpha_calc
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_alongspan_deg.value = alpha_calc(obj1, ...
                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cl_F.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p) - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
% alpha_calc_lin(obj1, ...
%                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cl_F.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value) - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_alongspan_deg.Attributes.unit = "degrees"; 

% Convert to radians 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_alongspan_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_alongspan_deg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_alongspan_rad.Attributes.unit = "radians";

% Changing folder 
cd .. 
cd csvla
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

%% TOTAL ANGLE OF ATTACK DATA STORAGE

% AoA_Tot = AoA + Twist_angle of the main wing
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_Tot_deg.value = Aircraft.Geometry.Wing.twist_angle.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.alpha_F.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_Tot_deg.Attributes.unit = "Degrees"; 

% Convert in radians 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_Tot_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_Tot_deg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_Tot_rad.Attributes.unit = "Radians";

%% PROJECTION OF FORCES ALONG WING AXES - NORMAL AND AXIAL FORCES 

% Calculation of the normal force coefficient
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCz.value = calc_normal_force(obj2, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_Tot_deg.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCl_distr.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCd_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCz.Attributes.unit = "m"; 

% Calculation of the axial force coefficient 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCa.value = calc_axial_force(obj2, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.AoA_Tot_deg.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCl_distr.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCd_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCa.Attributes.unit = "m"; 

% Normal force 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Normal_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCz.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Normal_force.Attributes.unit = "N/m";

% Axial force 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Axial_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cCa.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Axial_force.Attributes.unit = "N/m";

%% SHEAR FORCE CALCULATION 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Shear_distr.value = calc_shear_force(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Normal_force.value)*(1e-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Shear_distr.Attributes.unit = "daN";

%% BENDING MOMENT CALCULATION 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Bend_mom_distr.value = calc_bend_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Shear_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Bend_mom_distr.Attributes.unit = "daN*m";

%% SHEAR, BENDING MOMENT AND TORSION DIAGRAM

% Pitching moment per unit length
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.m_distr.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cm_F.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cm_F.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.m_distr.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cm_F.value(i)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value*((Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value(i))^2);
end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.m_distr.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.cm_F.value.*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.qF.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value.^2);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.m_distr.Attributes.unit = "N";

Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Tors_mom_distr.value = calc_tors_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.m_distr.value)*(1e-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Tors_mom_distr.Attributes.unit = "daN*m";

% Subplots with shear, bending moment and torsion
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Shear_BendMom_diagram.value = Shear_Bending_Torsion_diag(obj2, flip(Aircraft.Geometry.Wing.y.value), ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Shear_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Bend_mom_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Tors_mom_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value);
pause(1);
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointF.pdf', 'ContentType', 'vector')
% Saving figures inside correct folder
fprintf('Saving in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile ShearBendingTorsionDiagramPointF.pdf Output

%% POINT G CALCULATIONS                 
% Lift coefficient distribution along the span at the Point G
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cl_G.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CL_G.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cl_G.Attributes.unit = "Non dimensional";

% Interpolation to obtain Cd at the desired global CL
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Cd.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value, ...
                                                                        XI, YI, 'spline');
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Cd.Attributes.unit = 'Non dimensional';

% Pitching moment interpolation along the span
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Cm.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value, ...
                                                                        XI, YI, 'spline'); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Cm.Attributes.unit = 'Non dimensional';
           
% Integrate along semi-span and assign drag and pitching moment along span
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Global_CD.value = zeros(length(yi), 1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Global_CD.Attributes.unit = "Non dimensional";
for i = 1:length(yi)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Global_CD.value(i) = trapz(Aircraft.Geometry.Wing.half_span_y.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Cd.value(i,:));
    if abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Global_CD.value(i) - Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.CD_G.value) < 1e-2
       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cd_G.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Cd.value(i,:)';
       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cm_G.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Interpolated_Cm.value(i,:)';
    end
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cd_G.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cm_G.Attributes.unit = "Non dimensional";

%% IMPORT DATA FROM THE CORRECT SOURCE

%% PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
% In this section of the code two vectors are defined to store the product 
% c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCl_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cl_G.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCl_distr.Attributes.unit = "m";

Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCd_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cd_G.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCd_distr.Attributes.unit = "m";

%% ALPHA VECTOR ALONG THE MAIN WING SPAN 
cd .. 
cd utilities
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

% Calling the function alpha_calc
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_alongspan_deg.value = alpha_calc(obj1, ...
                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cl_G.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p) - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
%                          alpha_calc_lin(obj1, ...
%                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cl_G.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value) - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_alongspan_deg.Attributes.unit = "degrees"; 

% Convert to radians 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_alongspan_rad.value = rad2deg(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_alongspan_deg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_alongspan_rad.Attributes.unit = "Degrees";

% Changing folder 
cd .. 
cd csvla
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

%% TOTAL ANGLE OF ATTACK DATA STORAGE

% AoA_Tot = AoA + Twist_angle of the main wing
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_Tot_deg.value = Aircraft.Geometry.Wing.twist_angle.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.alpha_G.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_Tot_deg.Attributes.unit = "Degrees"; 

% Convert in radians 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_Tot_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_Tot_deg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_Tot_rad.Attributes.unit = "Radians";

%% PROJECTION OF FORCES ALONG WING AXES - NORMAL AND AXIAL FORCES 

% Calculation of the normal force coefficient
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCz.value = calc_normal_force(obj2, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_Tot_deg.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCl_distr.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCd_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCz.Attributes.unit = "m"; 

% Calculation of the axial force coefficient 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCa.value = calc_axial_force(obj2, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.AoA_Tot_deg.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCl_distr.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCd_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCa.Attributes.unit = "m"; 

% Normal force 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Normal_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCz.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Normal_force.Attributes.unit = "N/m";

% Axial force 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Axial_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cCa.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Axial_force.Attributes.unit = "N/m";

%% SHEAR FORCE CALCULATION 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Shear_distr.value = calc_shear_force(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Normal_force.value)*(1e-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Shear_distr.Attributes.unit = "daN";

%% BENDING MOMENT CALCULATION 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Bend_mom_distr.value = calc_bend_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Shear_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Bend_mom_distr.Attributes.unit = "daN*m";

%% SHEAR, BENDING MOMENT AND TORSION DIAGRAM

% Pitching moment per unit length
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.m_distr.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cm_G.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cm_G.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.m_distr.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cm_G.value(i)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value*((Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value(i))^2);
end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.m_distr.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.cm_G.value.*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.qG.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value.^2);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.m_distr.Attributes.unit = "N";

Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Tors_mom_distr.value = calc_tors_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.m_distr.value)*(1e-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Tors_mom_distr.Attributes.unit = "daN*m";

% Subplots with shear, bending moment and torsion diagram
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Shear_BendMom_diagram.value = Shear_Bending_Torsion_diag(obj2, flip(Aircraft.Geometry.Wing.y.value), ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Shear_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Bend_mom_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Tors_mom_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.point_name.value);
pause(1);
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointG.pdf', 'ContentType', 'vector')
% Saving figures inside correct folder
fprintf('Saving in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile ShearBendingTorsionDiagramPointG.pdf Output

%% POINT E CALCULATIONS                 
% Lift coefficient distribution along the span at the Point E
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cl_E.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CL_E.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cl_E.Attributes.unit = "Non dimensional";

% Interpolation to obtain Cd at the desired global CL
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Cd.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value, ...
                                                                        XI, YI, 'spline');
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Cd.Attributes.unit = 'Non dimensional';

% Pitching moment interpolation along the span
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Cm.value = interp2(XX, YY, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value, ...
                                                                        XI, YI, 'spline'); 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Cm.Attributes.unit = 'Non dimensional';
           
% Integrate along semi-span and assign drag and pitching moment along span
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Global_CD.value = zeros(length(yi), 1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Global_CD.Attributes.unit = "Non dimensional";
for i = 1:length(yi)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Global_CD.value(i) = trapz(Aircraft.Geometry.Wing.half_span_y.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Cd.value(i,:));
    if abs(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Global_CD.value(i) - Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.CD_E.value) < 1e-2
       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cd_E.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Cd.value(i,:)';
       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cm_E.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Interpolated_Cm.value(i,:)';
    end
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cd_E.Attributes.unit = "Non dimensional";
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cm_E.Attributes.unit = "Non dimensional";

%% IMPORT DATA FROM THE CORRECT SOURCE

%% PRODUCT C(y)*Cl(y) AND C(y)*Cd(y)
% In this section of the code two vectors are defined to store the product 
% c(y)*Cl(y) and c(y)*Cd(y) inside the struct variable 'Aircraft'
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCl_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cl_E.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCl_distr.Attributes.unit = "m";

Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCd_distr.value = times(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value', Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cd_E.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCd_distr.Attributes.unit = "m";

%% ALPHA VECTOR ALONG THE MAIN WING SPAN 
cd .. 
cd utilities
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

% Calling the function 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_alongspan_deg.value = alpha_calc(obj1, ...
                                       Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cl_E.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL0.value, ...
                                       Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
                                       Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
                                       p) - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
%                          alpha_calc_lin(obj1, ...
%                                         Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cl_E.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%                                         Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value) - Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_alongspan_deg.Attributes.unit = "degrees"; 

% Convert to radians 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_alongspan_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_alongspan_deg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_alongspan_rad.Attributes.unit = "radians";

% Changing folder 
cd .. 
cd csvla
% The 'dir' variable contains working directory path saved as a
% char value
dir = pwd;
% Store working directory inside the log file
fprintf('-----------------');
fprintf('\n');
fprintf('### Current directory ###');
fprintf('\n');
fprintf('%s\n', dir);

%% TOTAL ANGLE OF ATTACK DATA STORAGE

% AoA_Tot = AoA + Twist_angle of the main wing
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_Tot_deg.value = Aircraft.Geometry.Wing.twist_angle.value + Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.alpha_E.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_Tot_deg.Attributes.unit = "Degrees"; 

% Convert in radians 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_Tot_rad.value = deg2rad(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_Tot_deg.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_Tot_rad.Attributes.unit = "Radians";

%% PROJECTION OF FORCES ALONG WING AXES - NORMAL AND AXIAL FORCES 

% Calculation of the normal force coefficient
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCz.value = calc_normal_force(obj2, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_Tot_deg.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCl_distr.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCd_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCz.Attributes.unit = "m"; 

% Calculation of the axial force coefficient 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCa.value = calc_axial_force(obj2, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.AoA_Tot_deg.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCl_distr.value, ...
                                                                                  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCd_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCa.Attributes.unit = "m"; 

% Normal force 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Normal_force.value = N_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCz.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Normal_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCz.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Normal_force.Attributes.unit = "N/m";

% Axial force 
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.Axial_force.value = A_distr_along_wing(obj2, Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cCa.value, Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.qA.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Axial_force.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cCa.value*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Axial_force.Attributes.unit = "N/m";

%% SHEAR FORCE CALCULATION 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Shear_distr.value = calc_shear_force(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Normal_force.value)*(1e-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Shear_distr.Attributes.unit = "daN";

%% BENDING MOMENT CALCULATION 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Bend_mom_distr.value = calc_bend_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Shear_distr.value);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Bend_mom_distr.Attributes.unit = "daN*m";

%% SHEAR, BENDING MOMENT AND TORSION DIAGRAM

% Pitching moment per unit length
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.m_distr.value = zeros(length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cm_E.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cm_E.value)
    Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.m_distr.value(i) = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cm_E.value(i)*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value*((Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value(i))^2);
end
% Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.m_distr.value = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.cm_E.value.*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.qE.value*(Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value.^2);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.m_distr.Attributes.unit = "N";

Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Tors_mom_distr.value = calc_tors_mom(obj2, Aircraft.Geometry.Wing.y.value, ...
                                                                                Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.m_distr.value)*(1e-1);
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Tors_mom_distr.Attributes.unit = "daN*m";

% Subplots with shear, bending moment and torsion diagram
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Shear_BendMom_diagram.value = Shear_Bending_Torsion_diag(obj2, flip(Aircraft.Geometry.Wing.y.value), ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Shear_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Bend_mom_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Tors_mom_distr.value, ...
                                                                               Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.point_name.value);
pause(1);
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Shear_BendMom_diagram.value, 'ShearBendingTorsionDiagramPointE.pdf', 'ContentType', 'vector')
% Saving figures inside correct folder
fprintf('Saving in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile ShearBendingTorsionDiagramPointE.pdf Output

%% COMPARING SHEAR CURVES 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Shear_Comparison.value = Compare_Shear_curves(flip(Aircraft.Geometry.Wing.y.value), ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Shear_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Shear_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Shear_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Shear_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Shear_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Shear_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Shear_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.point_name.value);
% Saving
pause(1);
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Shear_Comparison.value, 'ShearComparison.pdf', 'ContentType', 'vector')
% Saving figures inside correct folder
fprintf('Saving in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile ShearComparison.pdf Output

%% COMPARING BENDING CURVES 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Bending_Comparison.value = Compare_Bending_curves(flip(Aircraft.Geometry.Wing.y.value), ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Bend_mom_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Bend_mom_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Bend_mom_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Bend_mom_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Bend_mom_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Bend_mom_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Bend_mom_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.point_name.value);
% Saving
pause(1);
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Bending_Comparison.value, 'BendingComparison.pdf', 'ContentType', 'vector')
% Saving figures inside correct folder
fprintf('Saving in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile BendingComparison.pdf Output

%% COMPARING TORSION CURVES 
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Torsion_Comparison.value = Compare_Torsion_curves(flip(Aircraft.Geometry.Wing.y.value), ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.Tors_mom_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.Tors_mom_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.Tors_mom_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.Tors_mom_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.Tors_mom_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.Tors_mom_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.Tors_mom_distr.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.point_name.value, ...
Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.point_name.value);
% Saving
pause(1);
exportgraphics(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Torsion_Comparison.value, 'TorsionComparison.pdf', 'ContentType', 'vector')
% Saving figures inside correct folder
fprintf('Saving in: ');
fprintf('\n'); 
fprintf('%s\n', SaveFolder);
% Moving file inside correct folder
movefile TorsionComparison.pdf Output
