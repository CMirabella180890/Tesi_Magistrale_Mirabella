
%% GLOBAL LIFT CALCULATOR

% STORE RESULTS FROM OPEN VSP DATA INTERPOLATION AND OTHER DATA
cl_interpolated = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_interpolated.value';
cd_interpolated = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Interpolated_Cd.value';
cm_interpolated = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.Interpolated_Cm.value';
half_span       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Yavg.value(1,:)'; 
S               = Aircraft.Geometry.Wing.S.value;
 
%% CHORD DISTRIBUTION 
kink1     = Aircraft.Geometry.Wing.chord_kink_one.value;
kink2     = Aircraft.Geometry.Wing.chord_kink_two.value;
wing_type = Aircraft.Geometry.Wing.type.value;
switch (wing_type)
    case 'Rectangular'
        % CHORD PARAMETERS
        croot       = Aircraft.Geometry.Wing.croot.value;
        ctip        = Aircraft.Geometry.Wing.ctip.value;
        taper_ratio = ctip / croot;

        % STORE INSIDE THE STRUCT VARIABLE
        Aircraft.Geometry.Wing.taper_ratio.value = taper_ratio;
        Aircraft.Geometry.Wing.taper_ratio.Attributes.unit = "Non dimensional";

        % Calculation of a chord distribution with a convenient, simple function.
        % 
        % c(y) = calc_chord(Swing, taper_ratio, span, y)
        % A complete documentation of this function is included inside the class
        % ShearBendingTorsion.m
        S           = Aircraft.Geometry.Wing.S.value;
        b           = Aircraft.Geometry.Wing.b.value;
        chord_distr = calc_chord(obj2, S, taper_ratio, b, half_span);

        % STORE INSIDE THE STRUCT VARIABLE
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value = chord_distr';
        Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.Attributes.unit = "m"; 
        chord_distr = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value;

    case 'With_kinks'
        if (kink1 ~= kink2)
            % WING FIRST SECTION - CHORDS
            ctip_section1        = Aircraft.Geometry.Wing.chord_kink_one.value;
            croot_section1       = Aircraft.Geometry.Wing.croot.value;
            taper_ratio_section1 = ctip_section1 / croot_section1;
            % WING SECOND SECTION - CHORDS
            ctip_section2        = Aircraft.Geometry.Wing.chord_kink_two.value;
            croot_section2       = Aircraft.Geometry.Wing.chord_kink_one.value;
            taper_ratio_section2 = ctip_section2 / croot_section2;
            % WING THIRD SECTION - CHORDS
            ctip_section3        = Aircraft.Geometry.Wing.ctip.value;
            croot_section3       = Aircraft.Geometry.Wing.chord_kink_two.value;
            taper_ratio_section3 = ctip_section3 / croot_section3;
            % STORE INSIDE THE AIRCRAFT STRUCT VARIABLE
            Aircraft.Geometry.Wing.taper_ratio1.value = taper_ratio_section1;
            Aircraft.Geometry.Wing.taper_ratio1.Attributes.unit = "Non dimensional";
            Aircraft.Geometry.Wing.taper_ratio2.value = taper_ratio_section2;
            Aircraft.Geometry.Wing.taper_ratio2.Attributes.unit = "Non dimensional";
            Aircraft.Geometry.Wing.taper_ratio3.value = taper_ratio_section3;
            Aircraft.Geometry.Wing.taper_ratio3.Attributes.unit = "Non dimensional";
            % Calculation of a chord distribution with a convenient, simple function.
            % 
            % c(y) = calc_chord(Swing, taper_ratio, span, y)
            % A complete documentation of this function is included inside the class
            % ShearBendingTorsion.m
            S           = Aircraft.Geometry.Wing.S.value;
            b           = Aircraft.Geometry.Wing.b.value;
            chord_distr1 = calc_chord(obj2, S, taper_ratio_section1, b, half_span(1:ceil(N/3)));
            chord_distr2 = calc_chord(obj2, S, taper_ratio_section2, b, half_span(ceil((N/3)+1):ceil(2*N/3)));
            chord_distr3 = calc_chord(obj2, S, taper_ratio_section3, b, half_span(ceil((2*N/3)+1):ceil(N/3)));
            % STORE INSIDE THE AIRCRAFT STRUCT VARIABLE
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value = [ chord_distr1'; chord_distr2'; chord_distr3' ];
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.Attributes.unit = "m";   
            chord_distr = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value;

        elseif (kink1 == kink2)
            % FROM ROOT TO KINK
            ctip_section1        = Aircraft.Geometry.Wing.chord_kink_one.value;
            croot_section1       = Aircraft.Geometry.Wing.croot.value;
            taper_ratio_section1 = ctip_section1 / croot_section1;
            % FROM KINK TO TIP
            ctip_section2        = Aircraft.Geometry.Wing.ctip.value;
            croot_section2       = Aircraft.Geometry.Wing.chord_kink_one.value;
            taper_ratio_section2 = ctip_section2 / croot_section2;

            % STORE INSIDE THE AIRCRAFT STRUCT VARIABLE                
            Aircraft.Geometry.Wing.taper_ratio1.value = taper_ratio_section1;
            Aircraft.Geometry.Wing.taper_ratio1.Attributes.unit = "Non dimensional";
            Aircraft.Geometry.Wing.taper_ratio2.value = taper_ratio_section2;
            Aircraft.Geometry.Wing.taper_ratio2.Attributes.unit = "Non dimensional";

            % Calculation of a chord distribution with a convenient, simple function.
            % 
            % c(y) = calc_chord(Swing, taper_ratio, span, y)
            % A complete documentation of this function is included inside the class
            % ShearBendingTorsion.m
            S           = Aircraft.Geometry.Wing.S.value;
            b           = Aircraft.Geometry.Wing.b.value;
            chord_distr1 = calc_chord(obj2, S, taper_ratio_section1, b, half_span(1:ceil(N/2)));
            chord_distr2 = calc_chord(obj2, S, taper_ratio_section2, b, half_span(ceil((N/2)+1):N));
            % STORE INSIDE THE AIRCRAFT STRUCT VARIABLE
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value = [ chord_distr1'; chord_distr2' ];
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.Attributes.unit = "m";   
            chord_distr = Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.chord_distr.value;
        end
end
% =================================================================

%% UNIT LIFT DISTRIBUTION 
chord_distr = chord_distr';
cCl_calc    = chord_distr.*cl_interpolated;
global_cl   = ones(length(cl_interpolated(1,:)), 1);
for i = 1:length(cl_interpolated(1,:))
    global_cl(i,1)  = (2/S) * trapz(half_span, cCl_calc(:,i));
end

tol           = 1e-2;
cl_unit_index = find(global_cl > 1.0-tol & global_cl < 1.0+tol);
cl_unit_distr = cl_interpolated(:, cl_unit_index);
if exist('cl_unit_index', 'var') == 0
   tol = 1e-1; 
   cl_unit_index = find(global_cl > 1.0-tol & global_cl < 1.0+tol);
   cl_unit_distr = cl_interpolated(:, cl_unit_index(end));
end
if exist('cl_unit_distr', 'var') == 0  
    error("ERROR: CL distribution along the wing semi-span does not contain the unit CL distribution!")
end
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.value           = cl_unit_distr;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.cl_at_CL1.Attributes.unit = "Non dimensional"; 
