%% CHAPTER 5 - Design Airspeeds
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

ch5 = Chapter();
ch5.Title = 'Design Airspeeds';

str = ['This chapter defines the operating and design airspeeds as required for certification REFFFF'];
para = Paragraph(str);
% append(para,InternalLink('tlarTableRef','refTabella'));
add(ch5,para)

%         %% CHAPTER 5 - SECTION 1
%         % VH
sec1 = Section();
sec1.Title = 'Maximum speed in level flight VH';

para = Paragraph(strcat('According to flight tests [5] at maximum weight and maximum continuous ',...
    'power at sea level conditions, the maximum speed in level flight has been determined:',...
    'V_H=130\ kts'));
% append(para,InternalLink('tlarTableRef','refTabella'));
para.Style = {HAlign('justify')};
add(sec1,para)

add(ch5,sec1);

%         %% CHAPTER 5 - SECTION 2
%         % VS, VS0, VS1
sec2 = Section();
sec2.Title = 'Stall speeds VS, VS0, VS1';
para = Paragraph(strcat('These speeds will be verified by flight test according to requirement 4.4.1 [1]',... 
'In order to calculate the stall speed, the maximum lift coefficient of the aeroplane as a whole is determined first.',...
'The maximum lift coefficient of the aeroplane has been calculated starting from the polar curve of the wing profile',...
'taken form ref. [6] (p. 236, Re=2.9E6 flaps retracted c_(L⁡_profile_max)=1.35 and p.237, δ_f=40 deg for the flaps',...
'in landing configurationc_(L⁡_profile_flapped_max)=2.15, and δ_f=10 deg in take-off configurationc_(L⁡_profile_flapped_to)=1.70).',...
'Considering the horizontal tail balancing force and the lower total wing lift due to wing lift distribution,',...
'the total aeroplane lift coefficient has been lowered by 15% with respect to the one of the profile.'));
% append(para,InternalLink('tlarTableRef','refTabella'));
para.Style = {HAlign('justify')};
add(sec2,para)

%VS
WMTOM   = 9.81*Aircraft.Weight.I_Level.W_maxTakeOff.value;
rho     = 1.225;
CLMAX   = max(Aircraft.Certification.Aerodynamic_data.CLMAX.value);
S       = Aircraft.Geometry.Wing.S.value; 
vs      = Aircraft.Certification.Regulation.SubpartC.Flightloads.print_positive_vs.value;
vs_un      = Aircraft.Certification.Regulation.SubpartC.Flightloads.print_positive_vs.Attributes.unit;

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
para = Paragraph('Flaps retracted(cleam configuration):');
para.Style = {HAlign('left')};
add(sec2,para)
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
append(sec2,eqImg);

%VS0
% latex interprete with $ simbol
para = Paragraph('Flaps extended(Landing configuration):');
% append(para,InternalLink('tlarTableRef','refTabella'));
para.Style = {HAlign('left')};
add(sec2,para)
myEq = "$ V_{S_0} = \sqrt{\frac{2 W_{MTOM}}{\rho_{0}C_{L_{MAX_{Landing}}}S}}";
eq = Equation(myEq);
eq.DisplayInline = true;
eq.FontSize = 12;
eqImg = getImpl(eq,rpt);
if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
    eqImg.Style = {VerticalAlign("-30%")};
elseif(rpt.Type == "docx") 
    eqImg.Style = {VerticalAlign("-5pt")};
end
append(sec2,eqImg);

%VS1
% latex interprete with $ simbol
para = Paragraph('Flaps extended(Take-off configuration):');
% append(para,InternalLink('tlarTableRef','refTabella'));
para.Style = {HAlign('left')};
add(sec2,para)
myEq = "$ V_{S_1} = \sqrt{\frac{2 W_{MTOM}}{\rho_{0}C_{L_{MAX_{Takeoff}}}S}}";
eq = Equation(myEq);
eq.DisplayInline = true;
eq.FontSize = 12;
eqImg = getImpl(eq,rpt);
if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
    eqImg.Style = {VerticalAlign("-30%")};
elseif(rpt.Type == "docx") 
    eqImg.Style = {VerticalAlign("-5pt")};
end
append(sec2,eqImg);

para = Paragraph(strcat('Therefore aeroplane lift coefficient is estimated to',...
        'c_(L_clean⁡_max)=0.85*1.35=1.15',...
        'and for the landing configuration (since the span extension of the flaps is half of the span of the wing):',...
        'c_(L_flaps_max)=(c_(L⁡_profile_flapped_max)+ c_(L_profile_max))/2*0.85=(2.15+1.35)/2*0.85=1.49'));
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec2,para)

%
para = Paragraph(strcat('The stall speed in landing configuration (flaps fully extended to ', 'xxx degrees', ' degrees)',...
            'is', 'XXX kts. ', 'Therefore it is In accordance with CS-LSA.5 [4]. ',...
            'In Take-Off configuration (flaps extended to , ', 'xxx degrees', 'degrees) the stall speed is , ', 'xxx', 'kts.'));
% append(para,InternalLink('tlarTableRef','refTabella'));
para.Style = {HAlign('left')};
add(sec2,para)

%
para = Paragraph(strcat('(Note: These speeds are estimates. The methods for the estimation can be various.',...
    'It is important that these estimations are as precise as possible. Flight tests will be used to validate ',...
    'the stall speeds. In case the flight tests show different values, this might have an impact on the speeds ',...
    'used for design and ultimately might impair the compliance to the CS-LSA.5.)'));
% append(para,InternalLink('tlarTableRef','refTabella'));
para.Style = {HAlign('left')};
para.BackgroundColor = "green";% = {Underline('yellow')};
add(sec2,para)



add(ch5,sec2);

%         %% CHAPTER 5 - SECTION 3
%         % VA
sec3 = Section();
sec3.Title = 'Design manoeuvring speed VA         ';

para = Paragraph('ADD TEXTS:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec3,para)



add(ch5,sec3);


%         %% CHAPTER 5 - SECTION 4
%         %
sec4 = Section();
sec4.Title = 'Flaps maximum operating speed VF';

para = Paragraph('ADD TEXTS:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec4,para)


add(ch5,sec4);

%         %% CHAPTER 5 - SECTION 5
%         %
sec5 = Section();
sec5.Title = 'Flaps maximum operating speed VFE';

para = Paragraph('ADD TEXTS:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec5,para)



add(ch5,sec5);

%         %% CHAPTER 5 - SECTION 6
%         %
sec6 = Section();
sec6.Title = 'Design cruising speed VC';

para = Paragraph('ADD TEXTS:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec6,para)


add(ch5,sec6);


%         %% CHAPTER 5 - SECTION 7
%         %
sec7 = Section();
sec7.Title = 'Design dive speed VD';

para = Paragraph('ADD TEXTS:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec7,para)


add(ch5,sec7);


%         %% CHAPTER 5 - SECTION 8
%         %
sec8 = Section();
sec8.Title = 'Demonstrated dive speed VDF';

para = Paragraph('ADD TEXTS:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec8,para)


add(ch5,sec8);

%         %% CHAPTER 5 - SECTION 9
%         %
sec9 = Section();
sec9.Title = 'Never exceed speed VNE';

para = Paragraph('ADD TEXTS:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec9,para)


add(ch5,sec9);



%% END chapter
%Adding chapters
add(rpt,ch5);