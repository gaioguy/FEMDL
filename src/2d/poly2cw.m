function p = poly2cw(p)

% orders the points in p clockwise

x = p(:,1); y = p(:,2);
cx = mean(x); cy = mean(y);
a = atan2(y-cy, x-cx);
[~, ord] = sort(a);
p(:,1) = x(ord); 
p(:,2) = y(ord);
