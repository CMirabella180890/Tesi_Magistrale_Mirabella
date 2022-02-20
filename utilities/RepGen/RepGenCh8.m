%% CHAPTER 8 - Manoeuvring and Gust load factors n
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

ch = Chapter();
ch.Title = 'V-n Envelope';
disp(['Chapter 8', (' "'), ch.Title,('" ') ,'writing...' ])
% -----------------------------------------------------------------------------------------------
requirement    = Aircraft.Certification.Regulation.value;
altitude       = Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.value;
altitude_unit  = Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.Attributes.unit;
flight_env_req = char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Attributes.cs);
% -------------------------------------------------------------------------
str = ['The final flight envelope diagram showed in the above figure was drawn' ...
    ' considering limit load factor values related to the maximum wing stress' ...
    ' condition. The determination of wing loads for structural design is important' ...
    ' to establish structural adequacy and to assess structural weight for the aircraft.' ...
    ' Aeroelasticity effects will be discounted; this is acceptable for relatively' ...
    ' slow, small aircraft. From this diagram, the designer can obtain load envelopes' ...
    ' in terms of shear, bending moment and torsion.'...
    ' The diagram is referred to sea level density and the airworthiness' ...
    ' rules applied are relative to ' ...
    (' ') ...
    char(requirement) ...
    '. The selected operational altitude is' ...
    (' ') ...
    num2str(altitude) ...
    (' ') ...
    char(altitude_unit) ...
    '. According to' ... 
    (' ') ...
    char(requirement) ...
    (' ') ...
    flight_env_req(1:end-1) ...
    ', each critical altitude expected during normal operations and' ...
    ' each critical weight and centre of gravity positions has been examined.' ...
    ' Symmetrical flight maneouvring conditions are defined as flight maneouvres' ...
    ' for which the aircraft is to be designed that does not involve' ...
    ' motion about the roll or yaw axis. Except where limited by maximum' ...
    ' static lift coefficients, the airplane is assumed to be subjected to' ...
    ' symmetrical maneouvres resulting in the limit maneouvring load' ...
    ' factors. The aeroplane is in steady-state flight condition,' ...
    ' where the pitching acceleration is assumed negligible.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(ch,para);
% -------------------------------------------------------------------------
        %table gust calculation        
        str = ['(Note: the' ...
            (' ') ...
            char(requirement) ...
            (' ') ...
            ' requirement does not impose to take into account for altitude' ...
            ' changes. Calculating the loads at sea level would be acceptable.'];
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        para.BackgroundColor = "green";
        add(ch,para);  
% -------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------
vs       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.value;
vs_un    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS.Attributes.unit;
ns       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.nS.value;

%gust calculation table - cruise
vc       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
vc_un    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.Attributes.unit;
nc       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.nC.value;

%gust calculation table - dive
vd       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
vd_un    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.Attributes.unit;
nd       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.nD.value;

%gust calculation table - manouevre
if isfield(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope, 'PointA1') == 1
    va       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.VA1.value;
    na       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.nA1.value;
else
    va       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
    na       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
end
% -----------------------------------------------------------------------------------------------
vs_inv   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.VS_inv.value;
ns_inv   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.nS_inv.value;
% -----------------------------------------------------------------------------------------------
vg       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG.value;
ng       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.nG.value;
% -----------------------------------------------------------------------------------------------
vf       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF.value;
nf       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.nF.value;
% -----------------------------------------------------------------------------------------------
ve       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.VE.value;
ne       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.nE.value;
% -----------------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------------
%moving to another path for figure
cd ..
cd ..
 regulation = Aircraft.Certification.Regulation.value;
 results_path = [pwd '\' regulation '\Output\'];

 cd (RepDir);

 fig = FormalImage([results_path,'Finalenvelope.png']);
 fig.Caption = 'Maneuver and Gust load factors diagram';
 fig.Height = '5in';
 fig.LinkTarget='maneuver_ref';
 add(ch,fig);
 % ------------------------------------------------------------------------- 
flight_env_airworthiness_req = char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.Attributes.cs);
str = ['The figure shows the final flight envelope, resulting from' ...
    ' superposition of manoeuvring and gust loads, according to ' ...
    (' ') ...
    char(requirement) ...
    (' ') ...
    flight_env_airworthiness_req ...
    ', obtaining the limit combined loads diagram. For the' ...
    ' calculation of structural design speeds, the stalling speeds' ...
    ' Vs0 and Vs1 should be taken to be the 1-g stalling speeds in' ...
    ' the appropriate flap configuration. '];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(ch,para);
% % -------------------------------------------------------------------------
% -------------------------------------------------------------------------
        %table gust calculation        
        str = ['According to' ...
            (' ') ...
            char(requirement) ...
            ', the load factor n varies linearly from points at cruise' ...
            ' speed VC (or VF for inverted flight) to points at dive speed VD' ...
            ' (D or E for inverted flight) as shown in figure.'];
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
%         para.BackgroundColor = "green";
        add(ch,para);  
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
        %table gust calculation        
        str = ['The gust load factors are determined from the previous' ...
            ' equation, known as the "gust formula"; the gust alleviation' ...
            ' factor Kg was determined as an empirical function of airplane' ...
            ' mass ratio mug. This analysis method is considered acceptable' ...
            ' for relatively slow, small, light airplane; more complicated' ...
            ' gust analysis requires a statistical approach.'];
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
%         para.BackgroundColor = "green";
        add(ch,para);  
% -------------------------------------------------------------------------
%  fig = FormalImage([results_path,'Finalenvelope.png']);
%  fig.Caption = 'Maneuver and Gust load factors and diagram';
%  fig.Height = '5in';
%  fig.LinkTarget='maneuver_ref';
%  add(ch,fig);
% -------------------------------------------------------------------------
        %table gust calculation        
        str = ['TABLE TO BE CHECKED!!!'];
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        para.BackgroundColor = "red";
        add(ch,para);  
        
        header = {'ID', strcat('V(',vc_un,')'), 'n'};
        %each table row needs of a fieldValue
        %gust calculation table
        %1 
        PointS     = {'S', num2str(vs,4), num2str(ns,4)};    
        PointA     = {'A', num2str(va,4), num2str(na,4)};     
        PointC     = {'C', num2str(vc,4), num2str(nc,4)};        
        PointD     = {'D', num2str(vd,4), num2str(nd,4)};
        %2 
        PointS_inv = {'S_inv', num2str(vs_inv,4), num2str(ns_inv,4)};    
        PointG     = {'G', num2str(vg,4), num2str(ng,4)};     
        PointF     = {'F', num2str(vf,4), num2str(nf,4)};        
        PointE     = {'E', num2str(ve,4), num2str(ne,4)};
        fieldValue = [PointS; PointA; PointC; PointD; PointS_inv; PointG; PointF; PointE];
    
          
        tbl = FormalTable(header,fieldValue);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
        tbl = BaseTable(tbl);
        tbl.Title = strcat('Final envelope points');
        tbl.LinkTarget = 'finalEnvelopeRef';
        add(ch,tbl);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% str = ['The following assumption has been made for gust calculations:'];
% para = Paragraph(str);
% para.Style = {HAlign('justify')};
% add(ch,para);
% % -------------------------------------------------------------------------
% % -------------------------------------------------------------------------
%         % AXIAL EQUILIBRIUM
%         %
%         myEq = "$ U = \frac{U_{de}}{2}\cdot \biggl( 1 - \cos{\frac{2\pi d_{gust}}{25\cdot MGC}} \biggr) ";
%         eq = Equation(strcat(myEq));
%         eq.DisplayInline = true;
%         eq.FontSize = 10;
%         eqImg = getImpl(eq,rpt);
%         if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
%             eqImg.Style = {VerticalAlign("-30%")};
%         elseif(rpt.Type == "docx")
%             eqImg.Style = {VerticalAlign("-5pt")};
%         end
%         append(ch,eqImg); 
% % -------------------------------------------------------------------------
% str = ['where'];
% para = Paragraph(str);
% para.Style = {HAlign('justify')};
% add(ch,para);
% % -------------------------------------------------------------------------
%         % d_gust
%         myEq = "$ d_{gust} = \mathrm{distance penetrated into gust; } (m) ";
%         eq = Equation(myEq);
%         eq.DisplayInline = true;
%         eq.FontSize = 12;
%         eqImg1 = getImpl(eq,rpt);
%         if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
%             eqImg.Style = {VerticalAlign("-30%")};
%         elseif(rpt.Type == "docx")
%             eqImg.Style = {VerticalAlign("-5pt")};
%         end
%         ref1 = eqImg1; 
%         
%         % MGC
%         myEq = "$ MGC = \mathrm{mean geometric chord of wing; } (m) ";
%         eq = Equation(myEq);
%         eq.DisplayInline = true;
%         eq.FontSize = 12;
%         eqImg2 = getImpl(eq,rpt);
%         if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
%             eqImg.Style = {VerticalAlign("-30%")};
%         elseif(rpt.Type == "docx")
%             eqImg.Style = {VerticalAlign("-5pt")};
%         end
%         ref2 = eqImg2; 
%         
%         % MGC
%         myEq = "$ U_{de} = \mathrm{derived gust velocity; } (m/s) ";
%         eq = Equation(myEq);
%         eq.DisplayInline = true;
%         eq.FontSize = 12;
%         eqImg3 = getImpl(eq,rpt);
%         if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
%             eqImg.Style = {VerticalAlign("-30%")};
%         elseif(rpt.Type == "docx")
%             eqImg.Style = {VerticalAlign("-5pt")};
%         end
%         ref3 = eqImg3;
%         
%         ol = UnorderedList({ref1, ref2, ref3});
% %         ol = UnorderedList({ref1,ref2,ref3,... 
% %             ref4,ref5,ref6, ref7,ref8});
% %         ol = UnorderedList({ref1, ref2, ref3,...
% %             ref4,ref5,ref6, ref7, ref8, ref9});
%         append(ch,ol);
% % ------------------------------------------------------------------------- 
% flight_env_airworthiness_req = char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Flight_envelope.Attributes.cs);
% str = ['The figure shows the final flight envelope, resulting from' ...
%     ' superposition of manoeuvring and gust loads, according to ' ...
%     (' ') ...
%     char(requirement) ...
%     (' ') ...
%     flight_env_airworthiness_req ...
%     ', obtaining the limit combined loads diagram. For the' ...
%     ' calculation of structural design speeds, the stalling speeds' ...
%     ' Vs0 and Vs1 should be taken to be the 1-g stalling speeds in' ...
%     ' the appropriate flap configuration. '];
% para = Paragraph(str);
% para.Style = {HAlign('justify')};
% add(ch,para);
% % % -------------------------------------------------------------------------
% Max_Lift_Coefficient          = Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value;
% Max_Lift_Coefficient_unit     = Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.Attributes.unit;
% Max_Lift_Coefficient_inverted = Aircraft.Certification.Aerodynamic_data.Min_Lift_Coefficient.value;
% CLMAX_takeoff                 = Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_takeoff.value;
% CLMAX_takeoff_unit            = Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_takeoff.Attributes.unit;
% CLMAX_landing                 = Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_landing.value;
% CLMAX_landing_unit            = Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_landing.Attributes.unit;
% CL0                           = Aircraft.Certification.Aerodynamic_data.CL0.value;
% CL0_unit                      = Aircraft.Certification.Aerodynamic_data.CL0.Attributes.unit;
% CLALFA_rad                    = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
% CLALFA_rad_unit               = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.Attributes.unit;
% CLALFA_deg                    = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;
% CLALFA_deg_unit               = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.Attributes.unit;
% CD0                           = Aircraft.Certification.Aerodynamic_data.CD0.value;
% CD0_unit                      = Aircraft.Certification.Aerodynamic_data.CD0.Attributes.unit;
% CDGEAR                        = Aircraft.Certification.Aerodynamic_data.CD_landing_gear.value;
% CDGEAR_UNIT                   = Aircraft.Certification.Aerodynamic_data.CD_landing_gear.Attributes.unit;
% CM0                           = Aircraft.Certification.Aerodynamic_data.CM0.value;
% CM0_unit                      = Aircraft.Certification.Aerodynamic_data.CM0.Attributes.unit;
% CMCL                          = Aircraft.Certification.Aerodynamic_data.CMCL.value;
% CMCL_unit                     = Aircraft.Certification.Aerodynamic_data.CMCL.Attributes.unit;
% CM_landing_gear               = Aircraft.Certification.Aerodynamic_data.CM_landing_gear.value;
% CM_landing_gear_unit          = Aircraft.Certification.Aerodynamic_data.CM_landing_gear.Attributes.unit;
% Alpha_zero_lift               = Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.value;
% Alpha_zero_lift_unit          = Aircraft.Certification.Aerodynamic_data.Alpha_zero_lift.Attributes.unit;
% Oswald_efficiency_factor      = Aircraft.Certification.Aerodynamic_data.e.value;
% Oswald_efficiency_factor_un   = Aircraft.Certification.Aerodynamic_data.e.Attributes.unit;
% % -------------------------------------------------------------------------
%         %table gust calculation        
%         str = ['TABLE TO BE CHECKED!!!'];
%         para = Paragraph(str);
%         para.Style = {HAlign('justify')};
%         para.BackgroundColor = "red";
%         add(ch,para);  
%         
%         header = {'Quantity', 'Value', 'Units'};
%         %each table row needs of a fieldValue
%         %gust calculation table
%         %1 
%         cl_wb_max_clean = {'CLMAX wb',         num2str(Max_Lift_Coefficient,4),          char(Max_Lift_Coefficient_unit)};    
%         cl_wb_inv_clean = {'CLMAXwb inverted', num2str(Max_Lift_Coefficient_inverted,4), char(Max_Lift_Coefficient_unit)};
%         cl_wb_takeoff   = {'CLMAX takeoff',    num2str(CLMAX_takeoff,4),                 char(CLMAX_takeoff_unit)};
%         cl_wb_landing   = {'CLMAX landing',    num2str(CLMAX_landing,4),                 char(CLMAX_landing_unit)};
%         cl_0            = {'CL0',              num2str(CL0,4),                           char(CL0_unit)};
%         cl_alfa_rad     = {'CLALFA (rad)',     num2str(CLALFA_rad,4),                    char(CLALFA_rad_unit)};
%         cl_alfa_deg     = {'CLALFA (deg)',     num2str(CLALFA_deg,4),                    char(CLALFA_deg_unit)};
%         cd0             = {'CD0',              num2str(CD0,4),                           char(CD0_unit)};
%         cd_gear         = {'CD landing gear',  num2str(CDGEAR,4),                        char(CDGEAR_UNIT)};
%         cm0             = {'CM0',              num2str(CM0,4),                           char(CM0_unit)};
%         cm_gear         = {'CM landing gear',  num2str(CM_landing_gear,4),               char(CM_landing_gear_unit)};
%         cmcl            = {'CMCL',             num2str(CMCL,4),                          char(CMCL_unit)};
%         alfa0l          = {'Alfa zero lift',   num2str(Alpha_zero_lift,4),               char(Alpha_zero_lift_unit)};
%         osw_eff_factor  = {'e',                num2str(Oswald_efficiency_factor,4),      char(Oswald_efficiency_factor_un)};
% 
%         fieldValue = [cl_wb_max_clean; cl_wb_inv_clean; cl_wb_takeoff; ...
%             cl_wb_landing; cl_0; cl_alfa_rad; cl_alfa_deg; cd0; cd_gear; ...
%             cm0; cm_gear; cmcl; alfa0l; osw_eff_factor];
%     
%           
%         tbl = FormalTable(header,fieldValue);
%         % In order to put a table with a caption, the API Report denomination should
%         % be used, the other options are from API DOM. In order to solve the problem,
%         % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
%         tbl = BaseTable(tbl);
%         tbl.Title = strcat('Aerodynamic reference values.');
%         tbl.LinkTarget = 'aerodynamictableRef';
%         add(ch,tbl);
% % -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% add(ch,para)
%% END chapter
%Adding chapters
add(rpt,ch);