% PROVA INTERPOLANTE
clear all; close all; clc; 

a           = -0.021837866;
b           = 0.436064773;
c           = -0.56312855; 
alpha_plus  = @(CL) (-b + sqrt(b^2 - 4*a*(c - CL)))/(2*a);

numb        = 1e3;
tol         = 1e-2;

alfa_data   = [-4; 0; 4; 8; 10; 12; 13];
CL_data     = [0.43483; 0.84054; 1.21541; 1.52777; 1.61373; 1.52500; 1.47190];

alfa_dot1   = [-4;       8];
CL_dot1     = [0.43483;  1.52777];

p1          = polyfit(alfa_dot1, CL_dot1, 1);

alfa_i      = - 4.0; 
alfa_f      = 13.0; 
alfa_interp = linspace(alfa_i, alfa_f, numb);

CL_star     = 1.47;
% CL_star     = 1.48;
alfa_star   = alpha_plus(CL_star);

alfa_dot2   = [alfa_star; 8;       10;      12;      13];
CL_dot2     = [CL_star;   1.52777; 1.61373; 1.52500; 1.47190];
p2          = polyfit(alfa_dot2, CL_dot2, 2);

CL_interp   = ones(length(alfa_interp), 1);

for i = 1:length(CL_interp)
    if alfa_interp(i) < alfa_star 
        CL_interp(i) = polyval(p1, alfa_interp(i));
    elseif alfa_interp(i) > alfa_star
        CL_interp(i) = polyval(p2, alfa_interp(i));
    end
end
        

figure; 
hold on; grid on; grid minor; 

plot(alfa_data, CL_data, 'k.', 'MarkerSize', 10)
plot(alfa_interp, CL_interp, 'r', 'LineWidth', 1.5)