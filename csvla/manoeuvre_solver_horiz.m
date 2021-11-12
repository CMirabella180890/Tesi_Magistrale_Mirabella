function [out1, out2, out3, out4, out5] = manoeuvre_solver_horiz(time_vector, d2thetadt2, ...
    A0, alpha_prime_horiz, delta_theta, dthetadt, time_step, delta_v, l_horiz_tail, ...
    damping_factor, V, alpha_new_horiz)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

     for i = 2:length(time_vector)
        d2thetadt2(i) = A0.value*(alpha_prime_horiz(i) - delta_theta(i-1));
        dthetadt(i) = dthetadt(i-1) + 0.5*(d2thetadt2(i-1)+d2thetadt2(i))*(time_step);
        delta_v(i) = dthetadt(i)*l_horiz_tail;
        delta_theta(i) = delta_v(i)*(damping_factor/V);
        alpha_new_horiz(i) = alpha_prime_horiz(i)  - delta_theta(i);
     end
     out1 = d2thetadt2;
     out2 = dthetadt;
     out3 = delta_v;
     out4 = delta_theta;
     out5 = alpha_new_horiz; 

end

