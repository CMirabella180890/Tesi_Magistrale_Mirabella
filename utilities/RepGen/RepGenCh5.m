%% CHAPTER 5 - Design Airspeeds
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

%% DECLARATION, DATA AND ASSUMPTIONTS CH5:

%reference AMC
requirement = strcat('-',Aircraft.Certification.Regulation.value);
% requirement = Aircraft.Certification.Regulation.value;

%%VS-CLMAX
CLMAX   = Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value; % CLMAX clean, CLMAX TO, CLMAX LAN
CLMAXLAN   = Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_landing.value; % CLMAX clean, CLMAX TO, CLMAX LAN
CLMAXTO   = Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_takeoff.value; % CLMAX clean, CLMAX TO, CLMAX LAN

%VS
WMTOM   = 9.81*Aircraft.Weight.I_Level.W_maxTakeOff.value;
rho     = Aircraft.Certification.ISA_Condition.Sea_Level.rho0.value; % sea level
S       = Aircraft.Geometry.Wing.S.value;
vs      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Stall_speed_VS.value;
vs_un      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Stall_speed_VS.Attributes.unit;

%VS0 landing stall speed
%DATA:
vslan = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointS.VS.value; %
vslan_un = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointS.VS.Attributes.unit; %

%VS1
%DATA:
vsto = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointS.VS.value; %
vsto_un = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointS.VS.Attributes.unit; %

%VD
if isnumeric(Aircraft.Certification.Regulation.SubpartC.Flightloads.Min_Design_Cruise_Speed.value)==1
    vcmin=Aircraft.Certification.Regulation.SubpartC.Flightloads.Min_Design_Cruise_Speed.value;
    disp("VCmin has been assigned.")
    disp(strcat("VCmin = ",num2str(vcmin),Aircraft.Certification.Regulation.SubpartC.Flightloads.Min_Design_Cruise_Speed.Attributes.unit))
else
    vcmin = 0.0; % m/s <-- ASSEGNED IF NOT AVAILABLE
    disp("VCmin has not been assigned.")
    disp("VCmin will be set to zero ")
end

%design airspeeds data structure
desspeeds.VS     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS;
desspeeds.VS0    = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointS.VS;
desspeeds.VS1    = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.PointS.VS;
desspeeds.VA     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA;
desspeeds.VC     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC;
desspeeds.VD     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD;
desspeeds.VE     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.VE;
%desspeeds.VF    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF;
desspeeds.VG     = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG;
desspeeds.VS_inv = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.VS_inv;
desspeeds.VF     = Aircraft.Certification.Regulation.SubpartC.Flapsloads.Landing.PointF.VF;

%% END DECLARATION

ch = Chapter();
ch.Title = 'Design Airspeeds';
disp(['Chapter 5', (' "'), ch.Title,('" ') ,'writing...' ])

%%
switch requirement
    
    % CASE 1: Very Light Aircraft
    case '-CSVLA'
        
        str = ['This chapter defines the operating and design airspeed as required for certification' ...
            (' ') ...
            requirement(2:end) ...
            '.'];
%         str = strcat('This chapter defines the operating and design airspeeds as required for certification ',...
%             requirement);
        para = Paragraph(str);
        % append(para,InternalLink('tlarTableRef','refTabella'));
        para.Style = {HAlign('justify')};
        add(ch,para)
        
        %         %% CHAPTER 5 - SECTION 1
        %         % VH
        sec = Section();
        sec.Title = 'Maximum speed in level flight VH';
        
        %VH
        % if not available assigned vh=0.0;    % to be checked
        if isnan(Aircraft.Certification.Regulation.SubpartC.Flightloads.Max_Continuous_Power_Speed_VH.value) == 1
            disp("VH is Not a Number")
            vh = 0.0; % m/s <-- ASSEGNED IF NOT AVAILABLE
            disp("VH has not been assigned.")
            disp("VH will be set to zero ")
            para = Paragraph(strcat('Data not yet available...to be added Available and Required Power.'));
            
            % or:
            % para = Paragraph(strcat('According to flight performance analysis, aerodynamic drag, ',...
            %     ' installed power at sea level conditions, the maximum speed in level flight has been determined:',...
            %     'V_H = xxxx m/s'));
            para.Style = {HAlign('justify')};
            add(sec,para)
            
        elseif isnumeric(Aircraft.Certification.Regulation.SubpartC.Flightloads.Max_Continuous_Power_Speed_VH.value)==1
            vh=Aircraft.Certification.Regulation.SubpartC.Flightloads.Max_Continuous_Power_Speed_VH.value;
            % to be checked
            disp("VH has been assigned.")
            disp(strcat("VH = ",num2str(vh),Aircraft.Certification.Regulation.SubpartC.Flightloads.Max_Continuous_Power_Speed_VH.Attributes.unit))
        else
            vh = 0.0; % m/s <-- ASSEGNED IF NOT AVAILABLE
            disp("VH has not been assigned.")
            disp("VH will be set to zero ")
            para = Paragraph(strcat('Data not yet available...to be added Available and Required Power.'));
            
            % or:
            % para = Paragraph(strcat('According to flight performance analysis, aerodynamic drag, ',...
            %     ' installed power at sea level conditions, the maximum speed in level flight has been determined:',...
            %     'V_H = xxxx m/s'));
            para.Style = {HAlign('justify')};
            add(sec,para)
        end
        %vh=0.0; %m/s
                
        add(ch,sec);
        
        %         %% CHAPTER 5 - SECTION 2
        %         % VS, VS0, VS1
        sec = Section();
        sec.Title = 'Stall speeds VS, VS0, VS1';
        para = Paragraph(strcat('These speeds will be verified by flight test according to certification requirements.',...
            ' In order to calculate the stall speed, the maximum lift coefficient of the aeroplane as a whole is determined first. ',...
            ' The maximum lift coefficient of the aeroplane has been calculated from high fidelity CFD. ',...
            ' In landing configuration computed with full flap, CLMAX landing = ',...
            num2str(CLMAXLAN),...
            ' in take-off configuration leading to CLMAX takeoff = ',...
            num2str(CLMAXTO),...
            ', and in clean configuration, leading to CLMAX clean = ',...
            num2str(CLMAX),...
            ' , also considering the horizontal tail balancing force.'));
        % append(para,InternalLink('tlarTableRef','refTabella'));
        para.Style = {HAlign('justify')};
        add(sec,para)
        
        myNumEq = strcat (' = \sqrt{\frac{2*',...
            num2str(WMTOM),...
            '}{',...
            num2str(rho),...
            '*',...
            num2str(CLMAX),...
            '*',...
            num2str(S),...
            '}} =',...
            num2str(vs),...
            num2str(vs_un));
        % latex interprete with $ simbol
        para = Paragraph('Flaps retracted(clean configuration):');
        para.Style = {HAlign('left')};
        add(sec,para)
        myEq = "$ V_{S} = \sqrt{\frac{2 W_{MTOM}}{\rho_{0}C_{L_{MAX_{Clean}}}S}}";
        eq = Equation(strcat(myEq, myNumEq));
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg);
        
        % latex interprete with $ simbol
        para = Paragraph('Flaps extended(Landing configuration):');
        % append(para,InternalLink('tlarTableRef','refTabella'));
        para.Style = {HAlign('left')};
        add(sec,para)
        
        myNumEq = strcat (' = \sqrt{\frac{2*',...
            num2str(WMTOM),...
            '}{',...
            num2str(rho),...
            '*',...
            num2str(CLMAXLAN),...
            '*',...
            num2str(S),...
            '}} =',...
            num2str(vslan),...
            num2str(vslan_un));
        % latex interprete with $ simbol
        
        myEq = "$ V_{S_0} = \sqrt{\frac{2 W_{MTOM}}{\rho_{0}C_{L_{MAX_{Landing}}}S}}";
        eq = Equation(strcat(myEq, myNumEq));
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg);
        
        % latex interprete with $ simbol
        para = Paragraph('Flaps extended(Take-off configuration):');
        % append(para,InternalLink('tlarTableRef','refTabella'));
        para.Style = {HAlign('left')};
        add(sec,para)
        
        myNumEq = strcat (' = \sqrt{\frac{2*',...
            num2str(WMTOM),...
            '}{',...
            num2str(rho),...
            '*',...
            num2str(CLMAXTO),...
            '*',...
            num2str(S),...
            '}} =',...
            num2str(vsto),...
            num2str(vsto_un));
        % latex interprete with $ simbol
        
        myEq = "$ V_{S_1} = \sqrt{\frac{2 W_{MTOM}}{\rho_{0}C_{L_{MAX_{Takeoff}}}S}}";
        eq = Equation(strcat(myEq, myNumEq));
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg);
        
        %
        para = Paragraph(strcat('Add here comments if necessary'));
        % append(para,InternalLink('tlarTableRef','refTabella'));
        para.Style = {HAlign('justify')};
        add(sec,para)
        
        %
        para = Paragraph(strcat('Note: These speeds are estimates. The methods for the estimation can be various.',...
            'It is important that these estimations are as precise as possible. Flight tests will be used to validate ',...
            'the stall speeds. In case the flight tests show different values, this might have an impact on the speeds ',...
            'used for design and ultimately might impair the compliance to the',...
            requirement));
        % append(para,InternalLink('tlarTableRef','refTabella'));
        para.Style = {HAlign('justify')};
        para.BackgroundColor = "green";% = {Underline('yellow')};
        add(sec,para)
        
        add(ch,sec);
        
        %         %% CHAPTER 5 - SECTION 3
        %         % VA
        sec = Section();
        sec.Title = 'Design manoeuvring speed VA';
airspeed_req = char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Cruise_Speed_VC.Attributes.cs);
str = ['According to ' ...
       (' ') ...
       requirement(2:end) ...
       (' ') ...
       airspeed_req(2:4) ...
       (' ') ...
       ', the maneuvring speed VA cannot be less then:'];        
%         para = Paragraph(strcat('According to requirement ',...
%             ' ',...
%             requirement,...
%             ' 335', ...
%             ' the maneuvering speed VA ', ...
%             ' cannot be less then:'));
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(sec,para)
        
%         para = Paragraph(strcat (' the maneuvering speed VA ',...
%             ' cannot be less then:'));
%         para.Style = {HAlign('justify')};
%         add(sec,para)
        
        %VA
        %DATA:
        nmax = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value;
        va   = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value; % VA
        va_un = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.Attributes.unit; %
        
        myNumEq = strcat (' = ',...
            num2str(vs),...
            ' * \sqrt{ ',...
            num2str(nmax),...
            ' } = ',...
            num2str(va),...
            num2str(va_un));
        % latex interprete with $ simbol
        myEq = "$ V_{A} \geq  V_{S}\sqrt{n_{max}} ";
        eq = Equation(strcat(myEq, myNumEq));
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg);
        para = Paragraph(strcat('Add here comments if necessary'));
        para.Style = {HAlign('justify')};
        % append(para,InternalLink('tlarTableRef','refTabella'));
        add(sec,para)
        
        
        add(ch,sec);
        
        
        %         %% CHAPTER 5 - SECTION 4
        %         %
        sec = Section();
        sec.Title = 'Flaps maximum operating speed VF';
        flaps_req = char(Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.VF.Attributes.cs);      
        str = ['According to ' ...
            (' ') ...
            requirement(2:end) ...
            (' ') ...
            flaps_req(1:4) ...
            ', such speed shall be not less than the greater of' ...
            ' 1.4*VS and 1.8*VS0.'];
%         para = Paragraph(strcat ('According to requirement ',...
%             requirement,...
%             ' -345, '));
%         add(sec,para)
%         %csvla and cs23
%         para = Paragraph(strcat ( ' such speed shall be not less than the greater of ',...
%             ' 1.4VS and',...
%             ' 1.8VS0'));
        para       = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(sec,para)
        %
        para = Paragraph(strcat('The speed has been selected as the greater between ',...
            '  1.4VS = ',...
            num2str(1.4*vs),...
            num2str(vs_un),...
            ' and 1.8 VSF = ',...
            num2str(1.4*vslan),...
            num2str(vslan_un),...
            ',',...
            ' where VSF is the computed stalling speed with flaps fully extended at the design weight.'));
        para.Style = {HAlign('justify')};
        add(sec,para)
        %
        para = Paragraph(strcat('The flaps operating speeds is: '));
        para.Style = {HAlign('justify')};
        add(sec,para)
        myNumEq = strcat ('V_{F} = ',...
            num2str(max(1.4*vs,1.4*vslan)),...
            num2str(vs_un));
        eq = Equation(myNumEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg);
        
        
        add(ch,sec);
        
        
        %         %% CHAPTER 5 - SECTION 5
        %         %
        sec = Section();
        sec.Title = 'Flaps maximum extension speed VFE';
        
        para = Paragraph(strcat('On this aeroplane the maximum flap extension speed is identical to the flap operating speed VF.',...
            ' This speed is the maximum speed for flaps in take-off and landing configuration.'));
        para.Style = {HAlign('justify')};
        add(sec,para)
        %
        myNumEq = strcat ('V_{FE} = ',...
            num2str(max(1.4*vs,1.4*vslan)),...
            num2str(vs_un));
        eq = Equation(myNumEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg);
        
        add(ch,sec);
        
        %         %% CHAPTER 5 - SECTION 6
        %         %
        sec = Section();
        sec.Title = 'Design cruising speed VC';
        
%         para = Paragraph(strcat('According to requirement ',...
%             requirement,...
%             '-335.'));
%         para.Style = {HAlign('justify')};
%         add(sec,para)
        %check on regultation for VC
        vc1 = 2.4*sqrt(Aircraft.Weight.I_Level.W_maxTakeOff.value*Aircraft.Constants.g.value/Aircraft.Geometry.Wing.S.value);
        vc2 = 0.9*vh;
%         para = Paragraph(strcat('- VC (in m/s) may not be less than -'));
        
        str        = ['According to ' ...
            (' ') ...
            requirement(2:end) ...
            (' ') ...
            airspeed_req(2:4) ...
            ', VC (in m/s) may not be less than'];
        para       = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(sec,para)
        % latex interprete with $ simbol
        myEq = "$2.4\sqrt{\frac{Mg}{S}}\left(V_{C}(kt)=4.7\sqrt{\frac{Mg}{S}}\right)";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        %    append(sec,eqImg);
        %
        %    myEq = "$2.4\sqrt{\frac{Mg}{S}}=";
        myNumEq = strcat('\rightarrow2.4\ast\sqrt{\frac{',...
            num2str(Aircraft.Weight.I_Level.W_maxTakeOff.value),...
            '\ast',...
            num2str(Aircraft.Constants.g.value),...
            '}{',...
            num2str(Aircraft.Geometry.Wing.S.value),...
            '}} = ',...
            num2str(vc1),...
            'm/s');
        eq = Equation(strcat(myEq,myNumEq));
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg);
        
        para = Paragraph(strcat('where M/S is the wing loading in kg/m^2 and g is the acceleration due to gravity in m/s^2.'));
        add(sec,para)
        %VCmin
        para = Paragraph(strcat('VC need not be more than 0.9*VH at sea level.'));
        add(sec,para);
        
        
        para = Paragraph(strcat('VH must be available. Otherwise previous value is considered!!!'));
        add(sec,para)
        
        %vc
        %
        myNumEq = strcat ('V_{C} = ',...
            num2str(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value),...
            num2str(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.Attributes.unit));
        eq = Equation(myNumEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        append(sec,eqImg);
        
        add(ch,sec);
        
        
        %         %% CHAPTER 5 - SECTION 7
        %         %
        sec = Section();
        sec.Title = 'Design dive speed VD';
        % -----------------------------------------------------------------
        str        = ['According to ' ...
            (' ') ...
            requirement(2:end) ...
            (' ') ...
            airspeed_req(2:4) ...
            ];
        para       = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(sec,para)
        % -----------------------------------------------------------------
%         para = Paragraph(strcat('According to requirement ',...
%             requirement,...
%             '-335.'));
%         add(sec,para)
        para = Paragraph(strcat('(1) VD may not be less than 1.25*VC; and',...
            ' (2) with VCmin, the required minimum design cruising speed, VD may not be less than 1.40*VCmin.'));
        add(sec,para)
        %1.25VC
        para = Paragraph(strcat('',...
            '(1)  1.25*VC = ',...
            num2str(1.25*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value),...
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.Attributes.unit));
        add(sec,para)
        %1.4VCmin
        para = Paragraph(strcat('',...
            '(2) 1.4VCmin = ',...
            num2str(1.25*vcmin),...
            Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.Attributes.unit));
        add(sec,para)
        
        if vcmin <= 1.25*Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value
            %vd
            myNumEq = strcat ('V_{C}=1.25\ast',...
                num2str(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value),...
                '=',...
                num2str(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value),...
                num2str(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.Attributes.unit));
            eq = Equation(myNumEq);
            eq.DisplayInline = true;
            eq.FontSize = 12;
            eqImg = getImpl(eq,rpt);
            if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
                eqImg.Style = {VerticalAlign("-30%")};
            elseif(rpt.Type == "docx")
                eqImg.Style = {VerticalAlign("-5pt")};
            end
            append(sec,eqImg);
        else
            %vd
            myNumEq = strcat ('V_{D}=1.25\ast',...
                num2str(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value),...
                '=',...
                num2str(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value),...
                num2str(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.Attributes.unit));
            eq = Equation(myNumEq);
            eq.DisplayInline = true;
            eq.FontSize = 12;
            eqImg = getImpl(eq,rpt);
            if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
                eqImg.Style = {VerticalAlign("-30%")};
            elseif(rpt.Type == "docx")
                eqImg.Style = {VerticalAlign("-5pt")};
            end
            append(sec,eqImg);
            
        end
        
        
        add(ch,sec);
        
        
        %         %% CHAPTER 5 - SECTION 8
        %         %
        sec = Section();
        sec.Title = 'Demonstrated dive speed VDF';
        
        para = Paragraph('VDF is not a design airspeeds for this category.');
        para.Style = {HAlign('justify')};
        add(sec,para)
        
        
        add(ch,sec);
        
        %         %% CHAPTER 5 - SECTION 9
        %         %
        sec = Section();
        sec.Title = 'Never exceed speed VNE';
        
        para = Paragraph('VNE is not a design airspeeds. It must be checked into sec. CSVLA 1505 Airspeed limitations.');
        para.Style = {HAlign('justify')};
        add(sec,para)
        
        
        add(ch,sec);
        
        %         %% CHAPTER 5 - SECTION 10
        %         %
        sec = Section();
        sec.Title = 'Design Airspeeds summary';
        
        para = Paragraph('Design airspeeds summary is resumed in Table:');
        append(para,InternalLink('speedsTableRef','Design airspeeds'));
        add(sec,para)
       %table
        speeds = fieldnames(desspeeds);
        fieldValue = cell(length(speeds),1);
        fieldUnit = cell(length(speeds),1);
        for i = 1:length(speeds)
            fieldValue{i} = desspeeds.(speeds{i}).value;
            % significant digits
            if isnumeric(fieldValue{i})
                fieldValue{i} = num2str(fieldValue{i},4);
            end
            % Not every field has Attributes. If not, they have only one field
            fieldUnit{i} = desspeeds.(speeds{i}).Attributes.unit;
        end
        
        speeds = [speeds, fieldValue, fieldUnit];
        header = {'Design airspeeds', 'Value', 'Measure unit'};
        tbl = FormalTable(header,speeds);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
        tbl = BaseTable(tbl);
        tbl.Title = 'Design airspeeds';
        tbl.LinkTarget = 'speedsTableRef';
        add(sec,tbl);
        
        
        add(ch,sec);
        
        % CASE: CS23
    case '-CS23'
        str = strcat('Check regulation to update ');
        str2 = strcat('Check regulation to update ');
        %1
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(ch,para)
        %2
        para2 = Paragraph(str2);
        para2.Style = {HAlign('justify')};
        add(ch,para2)
        
    % CASE 3: CS25
    case '-CS25'
        str = strcat('Check regulation to update ');
        str2 = strcat('Check regulation to update ');
        %1
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(ch,para)
        %2
        para2 = Paragraph(str2);
        para2.Style = {HAlign('justify')};
        add(ch,para2)
        
        % CASE CS22
    case '-CS22'
        str = strcat('Check regulation to update ');
        str2 = strcat('Check regulation to update ');
        %1
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(ch,para)
        %2
        para2 = Paragraph(str2);
        para2.Style = {HAlign('justify')};
        add(ch,para2)
        
        % CASE CSLSA
    case '-CSLSA'
        str = strcat('Check regulation to update ');
        str2 = strcat('Check regulation to update ');
        %1
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(ch,para)
        %2
        para2 = Paragraph(str2);
        para2.Style = {HAlign('justify')};
        add(ch,para2)
        
        
end



%% END chapter
%Adding chapters
add(rpt,ch);