%% CHAPTER 7 - Manoeuvring and Gust load factors n
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

ch = Chapter();
ch.Title = 'Manoeuvring and Gust load factors n';

str = ['Summary of limit load factors according to certification specifications and gust requirements.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};

add(ch,para)

%sec
sec = Section();
sec.Title = 'Gust envelope';
str = ['Gust load factors need to be considered because they can exceed'...
    'the prescribed maximum load factors at different weights and altitudes. '...
    'Since gust loads depend on air density and aircraft mass they will be '...
    'calculated for all twelve cases (sea level and 10000ft=FL100, maximum, '...
    'minimum flying weight and minimum flying weight with full wing fuel tanks) '...
    'according to requirement 5.2.3.3 [1] with flaps retracted (requirement 5.2.6.1 [1])'...
    'and fully extended (requirement 5.2.6.2 [1]) at V_F.'];
para = Paragraph(str);
para.Style = {HAlign('justify')};

add(sec,para);

str = ['The calculation is based on appendix X3 [1]. To calculate the gust '...
      'loads at altitudes other than at sea level the formula X3.1 [1] is '...
      ' altered to include the density at sea level \rho_0 as well:'];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);

str = [  'n_{3/4}=1\pm\frac{\frac{1}{2}\ \rho_0\ V\ K_g\ a\ U_{de}}{\left(\frac{w}{s}\right)}'];
para = Paragraph(str);
add(sec,para);

str = ['The corresponding weights are defined within [7]. '...
       'Since the gust loads on the wing and tail have been chosen to be' ...
       ' treated together, a is the slope of the lift-curve of the aeroplane '...
       '(\left(\frac{{dc}_L}{\alpha}\right)_{aeroplane}=4.77\frac{1}{rad}). '...
'(Note: the applicant should provide the method for the calculation of the slope of the lift-curve of the aeroplane) '];
para = Paragraph(str);
para.Style = {HAlign('justify')};
add(sec,para);
	
%moving to another path for figure
cd ..
cd ..
 regulation = Aircraft.Certification.Regulation.value;
 results_path = [pwd '\' regulation '\Output\'];

 cd (RepDir);
 
 fig = FormalImage([results_path,'Gustenvelope.png']);
 fig.Caption = 'Maneuver and Gust load factors and diagram';
 fig.Height = '5in';
 fig.LinkTarget='maneuver_ref';
 add(sec,fig);
 
 
 add(ch,sec);

%add(ch7,sec1);




%% END chapter
%Adding chapters
add(rpt,ch);