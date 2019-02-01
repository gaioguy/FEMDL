function t3 = delaunay_C2_new(p, box)

n = size(p,1);
dx = box(2)-box(1); 
dp = [dx 0; -dx 0];             
P = [p; p+dp(1,:); p+dp(2,:)];  % copy vertices to the left and the right
p_map = [1:n, 1:n, 1:n];

t = delaunay(P(:,1),P(:,2));    % construct triangulation
I = find(sum(t > 2*n,2) == 0);  % find triangles with at least one vertex in left copy of domain
t1 = t(I,:);                    % remove them from the triangulation
I = find(sum(t1 < n,2) > 0);    % find triangles with no vertex in original domain
t2 = t1(I,:);                   % remove them from the triangulation
t3 = p_map(t2);                 % map node numbers
