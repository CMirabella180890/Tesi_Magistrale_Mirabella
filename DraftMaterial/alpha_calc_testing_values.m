%% TEST TO TRY THE ALPHA_CALC FUNCTION INSIDE AERO_MODEL CLASS

%% Class instatiantion 
obj = aero_model; 

%% Function call 
% prova_alpha_calc_PuntoF = alpha_calc(obj, ...
%     Aircraft.Certification.Regulation.SubpartC.Final_envelope.PointF.CL_F.value, ...
%     Aircraft.Certification.Aerodynamic_data.CL0.value, ...
%     Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
%     Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
%     Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a, ...
%     Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b, ...
%     Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c);

prova_fromFtoE_alpha = zeros(length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value), 1);
for i = 1:length(Aircraft.Certification.Regulation.SubpartC.Final_envelope.Positive_stall_speed.value)
    
    % A complete documentation of the function alpha_calc(...) used here is
    % inside the class file aero_model.m, which can be found inside the
    % 'utilities' folder of this library.      
    prova_fromFtoE_alpha(i) = alpha_calc(obj, ...
    Aircraft.Certification.Regulation.SubpartC.Balancingloads.CL_fromFtoE.value(i), ...
    Aircraft.Certification.Aerodynamic_data.CL0.value, ...
    Aircraft.Certification.Aerodynamic_data.CL_star.value, ...
    Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value, ...
    Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_a, ...
    Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_b, ...
    Aircraft.Certification.Aerodynamic_data.Alpha_PolCoeff_c);
end