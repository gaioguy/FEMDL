function closed = is_closed(c)

% a contour c is closed if the distance between the first and the last
% point is less than twice the average distance between the points in the
% contour

md = mean(sum(diff(c).^2,2));  % mean distance between points
ed = sum((c(1,:)-c(end,:)).^2,2);
closed = (ed < 2*md) & size(c,1) > 2;
