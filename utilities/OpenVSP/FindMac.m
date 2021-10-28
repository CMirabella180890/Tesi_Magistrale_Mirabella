function [mac,xloc] = FindMac(wing)
% Finds mean aerodynamic chord and its location along X global axis.

npanels = length(wing.span);  % number of profile sections

for i = 1:npanels
    %#ok<*AGROW>
    lambda(i) = wing.ctip(i)/wing.croot(i);         % taper ratio(s)
    c(i) = 2*wing.croot(i)*(1+lambda(i)...          % panel mac
        +lambda(i)^2)/(3*(1+lambda(i)));
%     x(i)      = wing.span(i)*(1+...             % panel mac l.e. position
%         2*lambda(i))*tand(wing.sweep(i))/(6*(1+lambda(i)));
end

mac = wing.S * c' / sum(wing.S);    % mean aerodynamic chord

if npanels == 1
    if wing.ctip == wing.croot  % add exception for rectangular wing
        xloc = wing.xloc;
    else                        % tapered wing
        a = (wing.ctip-wing.croot)/(wing.span/2);
        b = wing.croot;
        y = (mac-b)/a;
        xloc = y*tand(wing.sweep(1)) + wing.xloc;
    end
elseif npanels == 2  % cranked wing
    if wing.S(1) > wing.S(2)
        a = (wing.ctip(1)-wing.croot(1))/(wing.span(1)/2);
        b = wing.croot(1);
        y = (mac-b)/a;
        xloc = y*tand(wing.sweep(1)) + wing.xloc;
    else
        a = (wing.ctip(2)-wing.croot(2))/wing.span(2)/2;
        b = wing.croot(2);
        % y = (mac-b)/a+wing.span(1)/2;
        xloc = wing.span(1)*tand(wing.sweep(1))/2 + ...
            (mac-b)/a*tand(wing.sweep(2)) + wing.xloc;
    end
else % if there are more than two panels, solve the integral
    if ~isfield(wing,'yloc')
        wing.yloc = 0;
    end
    y = [wing.yloc, wing.span];  % wing stations
    x(1) = wing.xloc;
    for i = 1:npanels
        x(i+1) = (y(i+1)-y(i)) * tand(wing.sweep(i)) + x(i);
    end
    c = [wing.croot, wing.ctip(end)];   % chord at y stations
    xloc = 2/sum(wing.S) * trapz(y,x.*c);
end

end