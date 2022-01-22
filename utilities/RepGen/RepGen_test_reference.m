<<<<<<< HEAD
%% bio reference example
import mlreportgen.dom.*
import mlreportgen.utils.*
d = Document('mydoc-reference','pdf');

% Append a link target to a heading
h = Heading(1,'Author''s Biography');
h.Style = {PageBreakBefore(true)};
linkID = normalizeLinkID('bio');
append(h,LinkTarget(linkID));

% Link to the target
append(d,InternalLink(linkID,'About the Author'));

% Append the heading 
append(d,h);



close(d);
=======
%% bio reference example
import mlreportgen.dom.*
import mlreportgen.utils.*
d = Document('mydoc-reference','pdf');

% Append a link target to a heading
h = Heading(1,'Author''s Biography');
h.Style = {PageBreakBefore(true)};
linkID = normalizeLinkID('bio');
append(h,LinkTarget(linkID));

% Link to the target
append(d,InternalLink(linkID,'About the Author'));

% Append the heading 
append(d,h);



close(d);
>>>>>>> eb25f8fea66ec6552ebf7b0bde47154c0a1a06a0
rptview(d);