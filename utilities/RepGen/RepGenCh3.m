%% CHAPTER 3 - List of Abbreviations
import mlreportgen.report.*  % import report API(report related methods
% @see https://it.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html?searchHighlight=mlreportgen.report&s_tid=doc_srchtitle#mw_63820826-78dc-459b-a646-67d4d77f91e5 )
import mlreportgen.dom.*     % import document object model DOM API (DOM related method
% @see https://it.mathworks.com/help/search.html?qdoc=mlreportgen.dom&submitsearch=)
import mlreportgen.utils.*

ch = Chapter();
ch.Title = 'List of Abbreviations';
disp(['Chapter 3', (' "'), ch.Title,('" ') ,'writing...' ])

str = ['ADD HERE list of abbreviations as a formatted table....to be created'];

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -------------------------------------------------------------------------
        % C.G.
        myEq = "$ \mathbf{c.g.} = \mathrm{centre of gravity;}";
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

        % Leading edge
        myEq = "$ \mathbf{l.e.} = \mathrm{leading edge;}";
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
        
        % ISA Atmosphere
        myEq = "$ \mathbf{i.s.a.} = \mathrm{international standard atmosphere;}";
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
        
        % keas 
        myEq = "$ \mathbf{keas} = \mathrm{knots equivalent speed;}";
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
        
        % mtow
        myEq = "$ \mathbf{mtow} = \mathrm{maximum takeoff weight;}";
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
        
        % mac
        myEq = "$ \mathbf{m.a.c.} = \mathrm{mean aerodynamic chord;}";
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
        
        % alpha
        myEq = "$ \alpha = \mathrm{angle of attack of the wing;}";
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
        
        % alpha_i
        myEq = "$ \alpha_{i} = \mathrm{local angle of attack of the station } i\mathrm{;}";
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
        
        % alpha_0
        myEq = "$ \alpha_{0} = \mathrm{wing zero lift angle of attack;}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg9 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref9 = eqImg9;
        
        % CLalfa
        myEq = "$ C_{L_{\alpha}} = \mathrm{lift curve slope } [1/\mathrm{deg}]";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg10 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref10 = eqImg10;
        
        % aspect ratio
        myEq = "$ \mathit{AR} = \mathrm{wing aspect ratio}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg11 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref11 = eqImg11;
        
        % wing span
        myEq = "$ b = \mathrm{wing span}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg12 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref12 = eqImg12;
        
        % average chord at station i
        myEq = "$ b = \mathrm{wing span}";
        eq = Equation(myEq);
        eq.DisplayInline = true;
        eq.FontSize = 12;
        eqImg13 = getImpl(eq,rpt);
        if (rpt.Type == "html" || rpt.Type == "html-file" || rpt.Type == "pdf")
            eqImg.Style = {VerticalAlign("-30%")};
        elseif(rpt.Type == "docx")
            eqImg.Style = {VerticalAlign("-5pt")};
        end
        ref13 = eqImg13;
        
        ol = UnorderedList({ref1,ref2,ref3,ref4,ref5,ref6,ref7, ...
            ref8,ref9,ref10,ref11,ref12,ref13});
%         ol = UnorderedList({ref1,ref2,ref3,...
%             ref4,ref5,ref6, ref7,ref8});
%         ol = UnorderedList({ref1, ref2, ref3,...
%             ref4,ref5,ref6, ref7, ref8, ref9});
%        append(subsec,ol);
% ------------------------------------------------------------------------- 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ------------------------------------------------------------------------- 
% % %List un-ordered
%  item1 = 'CL = lift coefficient';
%  item2 = 'CD....';
%  item3 = '...';
%  item4 = '...';
%  item5 = '...';
%  item6 = '...';
%  item7 = '...';
% 
% ol = UnorderedList({item1, item2, item3,...
%     item4, item5,item6,item7});

append(ch,ol);

para = Paragraph(str);
% append(para,InternalLink('tlarTableRef','refTabella'));
add(ch,para)
%% END chapter
%Adding chapters
add(rpt,ch);