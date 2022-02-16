%% CHAPTER 7 - Manoeuvring and Gust load factors n
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

%% DECLARATION, DATA AND ASSUMPTIONTS CH5:

%reference AMC
requirement = strcat(' ', Aircraft.Certification.Regulation.value);

nmax = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.value;
nmin = Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.value;

altitude = Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.value;
altitude_un  = Aircraft.Certification.ISA_Condition.Operative_ceiling.Altitude.Attributes.unit;

% a CNA slope
a = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.value;
a_un = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope.Attributes.unit;
adeg = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.value;
adeg_un = Aircraft.Certification.Aerodynamic_data.Normal_Force_Curve_Slope_deg.Attributes.unit;

%gust
gustc = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.value;
gustc_un = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_cruise.Attributes.unit;
gustd = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.value;
gustd_un = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_speed_dive.Attributes.unit;

%gust calculation table - cruise
vc       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.value;
vc_un    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC.Attributes.unit;
mtow     = Aircraft.Weight.I_Level.W_maxTakeOff.value;
mtow_un  = Aircraft.Weight.I_Level.W_maxTakeOff.Attributes.unit;
s        = Aircraft.Geometry.Wing.S.value;
s_un     = Aircraft.Geometry.Wing.S.Attributes.unit;
m_over_s = mtow/s;
rho      = Aircraft.Certification.ISA_Condition.Operative_ceiling.rho0.value;
rho_un   = Aircraft.Certification.ISA_Condition.Operative_ceiling.rho0.Attributes.unit;
mg       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Mass_ratio.value;
kg       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_alleviation_factor.value;
nc       = max (Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.value);

%gust calculation table - dive
vd       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.value;
vd_un    = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD.Attributes.unit;
nd       = max(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_dive.value);

%gust calculation table - manouevre
if isfield(Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope, 'PointA1') == 1
    va       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.VA1.value;
    na       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA1.nA1.value;
else
    va       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA.value;
    na       = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.nA.value;
end
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
    
%% END DECLARATION

ch = Chapter();
ch.Title = 'Manoeuvring and Gust load factors n';
disp(['Chapter 7', (' "'), ch.Title,('" ') ,'writing...' ])

switch requirement
    % CASE 1: Very Light Aircraft
    case 'CSVLA'
        %1
        airworth_rules1 = char(Aircraft.Certification.Regulation.SubpartC.Flightloads.nmax.Attributes.cs);
        airworth_rules2 = char(Aircraft.Certification.Regulation.SubpartC.Flightloads.nmin.Attributes.cs);
        str = ['According to ' ...
               (' ') ...
               requirement ...
               (' ') ...
               airworth_rules1 ...
               ', the positive limit moeuvring load factor n may not be less than 3.8, while ' ...
               'according to' ...
               (' ') ...
               requirement ...
               airworth_rules2 ...
               ', the negative limit manoeuvring load factor may not be less than -1.5.'];
        %2
        str2 = strcat('The following value will be considered: ');
        %ordered list
        ref1 = strcat('nmax = ', num2str(nmax));
        ref2 = strcat('nmin = ', num2str(nmin));
        ol = OrderedList({ref1, ref2});
        
        %1
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(ch,para)
        %2
        para2 = Paragraph(str2);
        para2.Style = {HAlign('justify')};
        add(ch,para2)
        %list
        append(ch,ol);
        
        
        %sec
        sec = Section();
        sec.Title = 'Gust envelope';
        %1
        str = ['Gust load factors need to be considered because they can exceed'...
            (' ')...
            'the prescribed maximum load factors at different weights and altitudes. '...
            (' ')...
            'Since gust loads depend on air density and aircraft mass they will be '...
            (' ')...
            'calculated for Compliance with the flight load requirements of this subpart to show:'];
        %
%         str2 = strcat('(1) At each critical altitude within the range in which the aeroplane may be expected to operate,',...
%             ' from sea level up to maximum operative altitude equal to:',...
%             strcat(num2str(altitude), altitude_un),...
%             '.');
        str2 = [' (1) At each critical altitude within the range in which the aeroplane may be expected to operate ' ...
                ' from sea level up to maximum operative altitude equal to: ' ...
                (' ') ...
                strcat(num2str(altitude)) ...
                (' ') ...
                char(altitude_un) ...
                '.'];
        str3 =['(2) At each practicable combination of weight and disposable load within the operating limitations specified in the Flight Manual'...
            (' ')...
            'according to requirement '...
            (' ') ...
            requirement...
            (' ') ...
            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Attributes.cs)...
            ' and fully extended (requirement'...
            (' ') ...
            requirement...
            (' ') ...
            char(Aircraft.Certification.Regulation.SubpartC.Flapsloads.Takeoff.VF.Attributes.cs) ...
            ' at VF).'];
        %1
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        %            ...       
        add(sec,para);
                %2
        para2 = Paragraph(str2);
        para2.Style = {HAlign('justify')};
        %            ...       
        add(sec,para2);
        %2
        para3 = Paragraph(str3);
        para3.Style = {HAlign('justify')};
        %            ...       
        add(sec,para3);
        
        
        str = ['The calculation is based on' ...
            (' ') ...
            requirement ...
            (' ') ...
            char(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.Attributes.cs) ...
            '. To calculate the gust '...
            'loads at altitudes other than at sea level the following equation is '...
            (' ')...
            ' altered to include the density at any altitude.'];
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(sec,para);
        
        %Gust equation
        %$ for tex interpreter
        myEq = "$ n = 1 + \frac{1/2 \rho_{0} V a K_{g} U_{de}}{Mg/S}";
        eq = Equation(myEq);
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
        str = ['where: '];
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(sec,para);
        %unordered lists
        %kg
        myEq = "$ K_{g} = \frac{0.88\mu_{g}}{5.3+\mu_{g}}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg1 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref1 = eqImg1;
        
        %mg
        myEq = "$ \mu_{g} = \frac{2(M/S)}{\rho\bar{C}a}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg2 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref2 = eqImg2;
        
        %ude
        myEq = "$U_{de} = \mathrm{derived gust velocities referred to in CSVLA 333(c) (m/s)}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg3 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref3 = eqImg3;
        
        %rho0
        myEq = "$\rho_{0} = \mathrm{density of air at sea level (kg/m3)}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg4 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref4 = eqImg4;
        
        %rho
        myEq = "$\rho = \mathrm{density of air (kg/m3)}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg5 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref5 = eqImg5;

        %M/S
        myEq = "$M/S = \mathrm{wing loading (kg/m2)}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg6 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref6 = eqImg6;   
       
        %cbar
        myEq = "$\bar{C} = \mathrm{mean geometric chord (m)}; g = \mathrm{acceleration due to gravity (m/s2);}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg7 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref7 = eqImg7; 
        
%cbar
        myEq = "$V = \mathrm{aeroplane equivalent speed (m/s); and}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg8 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref8 = eqImg8;         

        %a
        myEq = "$a = \mathrm{slope of the aeroplane normal force coefficient curve CNA per radian}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg8 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref8 = eqImg8;        
        
        ol = UnorderedList({ref1,ref2,ref3,...
            ref4,ref5,ref6, ref7,ref8});
%         ol = UnorderedList({ref1, ref2, ref3,...
%             ref4,ref5,ref6, ref7, ref8, ref9});
        append(sec,ol);
       
        str = strcat('Since the gust loads on the wing and tail have been chosen to be', ...
            ' treated together, a is the slope of the lift-curve of the aeroplane is equal to a = ',...
            num2str(a,4),...
            a_un,...
            ' and ',...
            num2str(adeg,4),...
            adeg_un,...
            '.');        
%         str = [ 'Since the gust loads on the wing and tail have been chosen to be trated together, ' ...
%                 'a is the slope of the lift-curve of the aeroplane which is equal to ' ...
%                 (' ') ...
%                 strcat(num2str(a,4)) ...
%                 (' ')
%                 strcat(char(a_un)) ...
%                 ' and ' ...
%                 (' ') ...
%                 strcat(num2str(adeg,4)) ...
%                 (' ') ...
%                 strcat(char(adeg_un)) ...
%                 '.' ];
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(sec,para);
        
        %gust@vc
        str = ['The gust speed at VC is equal to:'...
            (' ')...
            num2str(gustc)...
            gustc_un];
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(sec,para);
        %gust@vd
        str = ['The gust speed at VD is equal to:'...
            (' ')...
            num2str(gustd)...
            gustd_un];
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(sec,para);                
        
        %table gust calculation        
        str = ['TABLE TO BE CHECKED!!!'];
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        para.BackgroundColor = "red";
        add(sec,para);  
        
        header = {'ID',strcat('V(',vc_un,')'), strcat('M(',mtow_un,')'),...
                strcat('M/S(',mtow_un,'/',s_un,')'), strcat('Altitude(',altitude_un,')'),...
                strcat('rho(',rho_un,')'),'mug', 'Kg', strcat('Ude(',gustc_un,')') , 'n'};
        %each table row needs of a fieldValue
        %gust calculation table
        %1 
        PointA     = {'A', num2str(va,4), num2str(mtow,4), ...
                    num2str(m_over_s,4), num2str(altitude,4),...
                    num2str(rho,4), num2str(mg,4), num2str(kg,4), ...
                    num2str(gustc,4), num2str(na,4)};     
        %2 
        PointC     = {'C', num2str(vc,4), num2str(mtow,4), ...
                    num2str(m_over_s,4), num2str(altitude,4),...
                    num2str(rho,4), num2str(mg,4), num2str(kg,4), ...
                    num2str(gustc,4), num2str(nc,4)};        
        %3 
        PointD     = {'D', num2str(vd,4), num2str(mtow,4), ...
                    num2str(m_over_s,4), num2str(altitude,4), ...
                    num2str(rho,4), num2str(mg,4), num2str(kg,4), ...
                    num2str(gustd,4), num2str(nd,4)};
        %1 
        PointG     = {'G', num2str(vg,4), num2str(mtow,4), ...
                    num2str(m_over_s,4), num2str(altitude,4),...
                    num2str(rho,4), num2str(mg,4), num2str(kg,4), ...
                    num2str(gustc,4), num2str(ng,4)};     
        %2 
        PointF     = {'F', num2str(vf,4), num2str(mtow,4), ...
                    num2str(m_over_s,4), num2str(altitude,4),...
                    num2str(rho,4), num2str(mg,4), num2str(kg,4), ...
                    num2str(gustc,4), num2str(nf,4)};        
        %3 
        PointE     = {'E', num2str(ve,4), num2str(mtow,4), ...
                    num2str(m_over_s,4), num2str(altitude,4), ...
                    num2str(rho,4), num2str(mg,4), num2str(kg,4), ...
                    num2str(gustd,4), num2str(ne,4)};
        fieldValue = [PointA; PointC; PointD; PointG; PointF; PointE];
    
          
        tbl = FormalTable(header,fieldValue);
        % In order to put a table with a caption, the API Report denomination should
        % be used, the other options are from API DOM. In order to solve the problem,
        % the table is created as FormalTable (DOM) but it is inserted in a BaseTable (Report).
        tbl = BaseTable(tbl);
        tbl.Title = strcat('Gust load factor, different Speeds and Altitude');
        tbl.LinkTarget = 'gustcTableRef';
        add(sec,tbl);
        
        %note
        str = [  '(Note: the applicant should provide the method for the calculation of the slope of the lift-curve of the aeroplane) '];
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        para.BackgroundColor = "green";% = {Underline('yellow')};
        add(sec,para)
        
        
        %moving to another path for figure
        cd ..
        cd ..
        regulation = Aircraft.Certification.Regulation.value;
        results_path = [pwd '\' regulation '\Output\'];
        
        cd (RepDir);
        
        fig = FormalImage([results_path,'Vndiagram.png']);
        fig.Caption = 'V-n diagram';
        fig.Height = '4in';
        fig.Width = '4in';
        fig.LinkTarget='V-n diagram';
        add(sec,fig);
        
        fig = FormalImage([results_path,'Gustenvelope.png']);
        fig.Caption = 'Gust diagram';
        fig.Height = '4in';
        fig.Width = '4in';
        fig.LinkTarget='Gust diagram';
        add(sec,fig);
        % -----------------------------------------------------------------
        Gust_load_pos_cruise = max(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_pos_cruise.value);
        Gust_load_neg_cruise = min(Aircraft.Certification.Regulation.SubpartC.Flightloads.Gustloads.Gust_load_neg_cruise.value);
        str = ['In the last figure, it is possible to observe the gust' ...
            ' load factor varying linearly along the flight envelope' ...
            ' airspeed; the maximum and minimum gust load factor are,' ...
            ' respectively,' ...
            (' ') ... 
            num2str(Gust_load_pos_cruise, 4) ...
            ' and ' ...
            num2str(Gust_load_neg_cruise, 4) ...
            '.'];
        para = Paragraph(str);
        para.Style = {HAlign('justify')};
        add(sec,para);  
        % -----------------------------------------------------------------
        add(ch,sec);
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        % CASE 2: CS23
        % -----------------------------------------------------------------
    case 'CS23'
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
    case 'CS25'
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
    case 'CS22'
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
    case 'CSLSA'
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