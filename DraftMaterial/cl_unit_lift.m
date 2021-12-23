function cl = cl_unit_lift(CL1, CL2, cl1, cl2)

% LIFT COEFFICIENT DISTRIBUTION ALONG THE SPAN, WITH GLOBAL CL = 1.0 
% A simple function to evaluate the lift coefficient distribution along the
% span cl = cl(y) when the associated global lift coefficient of the whole
% wing is equal to 1.0; the function use a method similar to that suggested
% by Abbott in Theory of Wing Section. 
%   
%   To evaluate the requested lift coefficient distribution along the span
%   the function use the following formula: 
%   
%   IF CL1 IS CLOSER TO CL = 1.0 THAN CL2: 
%   
%   cl_unit_lift = 0.5*(k2*cl1 + k1*cl2)
%
%   IFT CL2 IS CLOSER TO CL = 1.0 THAN CL1: 
% 
%   cl_unit_lift = 0.5*(k2*cl1 + k1*cl2)
%   
%   INPUT 
%   CL1 --> Wing global lift coefficient, with CL1 < 1.0
%   CL2 --> Wing global lift coefficient, with CL2 > 1.0
%   cl1 --> Lift coefficient distribution along the span associated with
%           the global lift coefficient CL1
%   cl2 --> Lift coefficient distribution along the span associated with
%           the global lift coefficient CL2
%   OUTPUT 
%   cl  --> Lift coefficient distribution along the span associated with 
%           the global lift coefficient CL = 1.0

% Distance between CL2 and CL1
check0 = abs(CL2 - CL1); 

% Variables to check different cases for the interpolation
check1 = abs(1.0 - CL1); 
check2 = abs(1.0 - CL2);

% Coefficient that must be applied to the known lift distr. along the span
k1 = 1.0 + check1;
k2 = 1.0 + check2; 
    
% Condition to evaluate lift distribution along the span when CL = 1.0 
cl = 0.5*(k2*cl1 + k1*cl2);

end

