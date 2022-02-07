function [a,b] = chord_linear(c0 ,c1, span1, span0)
%CHORD_LINEAR Summary of this function goes here
%   linear chord law variation
% a = known terms (initial chord)
% b = slope --> c2-c1/(span2 - span1)
a = c0;
b = (c1 - c0)/(span1 - span0);
end

