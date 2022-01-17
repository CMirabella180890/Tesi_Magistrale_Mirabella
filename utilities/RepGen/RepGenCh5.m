%% CHAPTER 5 - Design Airspeeds
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

ch = Chapter();
ch.Title = 'Design Airspeeds';
%design airspeeds data structure

desspeeds.VS = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS.VS;
desspeeds.VA = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointA.VA;
desspeeds.VC = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointC.VC;
desspeeds.VD = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointD.VD;
desspeeds.VE = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointE.VE;
%desspeeds.VF = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointF.VF;
desspeeds.VG = Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointG.VG;
desspeeds.VS_inv =  Aircraft.Certification.Regulation.SubpartC.Flightloads.Final_envelope.PointS_inv.VS_inv;


%reference AMC
requirement = strcat('-',Aircraft.Certification.Regulation.value);

str = strcat('This chapter defines the operating and design airspeeds as required for certification ',...
            requirement);
para = Paragraph(str);
% append(para,InternalLink('tlarTableRef','refTabella'));
add(ch,para)

%         %% CHAPTER 5 - SECTION 1
%         % VH
sec = Section();
sec.Title = 'Maximum speed in level flight VH';
para = Paragraph(strcat('Data not yet available...to be added Available and Required Power'));
% or:
% para = Paragraph(strcat('According to flight performance analysis, aerodynamic drag, ',...
%     ' installed power at sea level conditions, the maximum speed in level flight has been determined:',...
%     'V_H = xxxx m/s'));
para.Style = {HAlign('justify')};
add(sec,para)

add(ch,sec);

%         %% CHAPTER 5 - SECTION 2
%         % VS, VS0, VS1
sec = Section();
sec.Title = 'Stall speeds VS, VS0, VS1';
para = Paragraph(strcat('These speeds will be verified by flight test according to certification requirements',... 
'In order to calculate the stall speed, the maximum lift coefficient of the aeroplane as a whole is determined first. ',...
' The maximum lift coefficient of the aeroplane has been calculated from high fidelity CFD. ',...
' ',...
' In landing configuration computed with full flap, CLMAX landing = 2.1, in take-off configuration leading to CLMAX takeoff = 1.9, and in clean configuration, leading to CLMAX clean = 1.61',...
' , also considering the horizontal tail balancing force.'));
% append(para,InternalLink('tlarTableRef','refTabella'));
para.Style = {HAlign('justify')};
add(sec,para)

%VS
WMTOM   = 9.81*Aircraft.Weight.I_Level.W_maxTakeOff.value;
rho     = 1.225; % sea level
CLMAX   = Aircraft.Certification.Aerodynamic_data.Max_Lift_Coefficient.value; % CLMAX clean, CLMAX TO, CLMAX LAN
S       = Aircraft.Geometry.Wing.S.value; 
vs      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Stall_speed_VS.value;
vs_un      = Aircraft.Certification.Regulation.SubpartC.Flightloads.Positive_Stall_speed_VS.Attributes.unit;

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

%VS0 landing stall speed
%DATA:
CLMAXLAN   = Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_landing.value; % CLMAX clean, CLMAX TO, CLMAX LAN
vslan = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value; %to be updated
vslan_un = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.Attributes.unit; %to be updated

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

%VS1
%DATA:
CLMAXTO   = Aircraft.Certification.Aerodynamic_data.Flaps.CLMAX_takeoff.value; % CLMAX clean, CLMAX TO, CLMAX LAN
vsto = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.value; %to be updated
vsto_un = Aircraft.Certification.Regulation.SubpartC.Flapsloads.PointS.VS.Attributes.unit; %to be updated

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

para = Paragraph(strcat('Add here comments if necessary'));
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec,para)

%
para = Paragraph(strcat('Add here comments if necessary'));
% append(para,InternalLink('tlarTableRef','refTabella'));
para.Style = {HAlign('left')};
add(sec,para)

%
para = Paragraph(strcat('(Note: These speeds are estimates. The methods for the estimation can be various.',...
    'It is important that these estimations are as precise as possible. Flight tests will be used to validate ',...
    'the stall speeds. In case the flight tests show different values, this might have an impact on the speeds ',...
    'used for design and ultimately might impair the compliance to the CS-LSA.5.)'));
% append(para,InternalLink('tlarTableRef','refTabella'));
para.Style = {HAlign('left')};
para.BackgroundColor = "green";% = {Underline('yellow')};
add(sec,para)

add(ch,sec);

%         %% CHAPTER 5 - SECTION 3
%         % VA
sec = Section();
sec.Title = 'Design manoeuvring speed VA         ';

para = Paragraph(strcat ('According to requirement ',...
                    ' ',...
                    requirement,...
                    ' -333, '));
add(sec,para)
%
para = Paragraph(strcat (' the maneuvering speed VA ',...
                 ' is computed as follow'));
add(sec,para)
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
myEq = "$ V_{A} = V_{S} \sqrt{n_{max}}";
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
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec,para)




add(ch,sec);


%         %% CHAPTER 5 - SECTION 4
%         %
sec = Section();
sec.Title = 'Flaps maximum operating speed VF';

para = Paragraph(strcat ('According to requirement ',...
                    requirement,...
                    ' -345, '));
add(sec,para)
%csvla and cs23
para = Paragraph(strcat ( ' such speed shall be not less than the greater of ',...
                        ' 1.4VS and',...
                        ' 1.8VS0'));                  
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
add(sec,para)
%
para = Paragraph(strcat('The flaps operating speeds is: '));
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

if requirement == "-CSVLA"
    para = Paragraph(strcat('On this aeroplane the maximum flap extension speed is identical to the flap operating speed VF.',...
        ' This speed is the maximum speed for flaps in take-off and landing configuration.'));  
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
    
else
    para = Paragraph(strcat('Check regulations.'));
end
add(sec,para)

add(ch,sec);

%         %% CHAPTER 5 - SECTION 6
%         %
sec = Section();
sec.Title = 'Design cruising speed VC';

para = Paragraph(strcat('According to requirement ',...
    'ADD req. par',...
    ' V_C may not be less than: ',...	
    'V_{Cmin}=4.77\sqrt{\frac{\ W_{MTOM}}{S}}=\ 4.77\ \sqrt{\frac{600\ast9.81}{15.1\ }}=94\ kts',...
    'and need not be greater than: ',...
    'V_{Cmax}={0.9\ V}_H=0.9\ast140\ =126\ kts',...
    'The speed has been selected as: ',...
    ' V_C=120\ kts '));
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec,para)


add(ch,sec);


%         %% CHAPTER 5 - SECTION 7
%         %
sec = Section();
sec.Title = 'Design dive speed VD';

para = Paragraph(strcat('According to requirement ',...
    'ADD req par',...
    'V_D=1.4\ V_{Cmin}=1.4\ \ast94\ =132\ kts',...
    'For V_D a higher value than the one above has been chosen: ',...
    'V_D=160\ kts'));
add(sec,para)


add(ch,sec);


%         %% CHAPTER 5 - SECTION 8
%         %
sec = Section();
sec.Title = 'Demonstrated dive speed VDF';

para = Paragraph('ADD TEXTS:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec,para)


add(ch,sec);

%         %% CHAPTER 5 - SECTION 9
%         %
sec = Section();
sec.Title = 'Never exceed speed VNE';

para = Paragraph('ADD TEXTS:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec,para)


add(ch,sec);

%         %% CHAPTER 5 - SECTION 10
%         %
sec = Section();
sec.Title = 'Design Airspeeds summary';

para = Paragraph('Design airspeeds summary is resumed in Table:');
% append(para,InternalLink('tlarTableRef','refTabella'));
add(sec,para)

speeds = fieldnames(desspeeds);
fieldValue = cell(length(speeds),1);
fieldUnit = cell(length(speeds),1);
for i = 1:length(speeds)
    fieldValue{i} = desspeeds.(speeds{i}).value;
    % significant digits
    if isnumeric(fieldValue{i})
        fieldValue{i} = num2str(fieldValue{i},5);
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



%% END chapter
%Adding chapters
add(rpt,ch);