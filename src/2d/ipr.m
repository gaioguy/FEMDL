function [r,l,a] = ipr(c)

%% IPR iosperimetric ratio
%
% [r,l,a] = IPR(c) computes the isoperimetric ratio for the curve c
%   c: (n x 2), discrete curve, one node per row
%   r: scalar, isoperimetric ratio l/a
%   l: scalar, length of the curve
%   a: scalar, area enclosed by the curve
%
% (C) 2019 by O. Junge, see COPYRIGHT 

c = poly2cw(c);                 % order points clockwise
l = clength(c);                 % length of the curve
a = polyarea(c(:,1),c(:,2));    % enclosed area
r = l/a;






