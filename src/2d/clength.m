function l = clength(c)

% returns the length of a curve c in the plane

dx = diff(c(:,1));
dy = diff(c(:,2));
l = sum(sqrt(dx.^2 + dy.^2));