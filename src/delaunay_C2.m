function t = delaunay_C2(p, dx)

%% DELAUNAY_C2 construct Delaunay triangulation on 2d cylinder
%
% t = delaunay_C2(p, dx) constructs a Delaunay triangulation on the 2d flat
%   cylinder (i.e. periodic in x-direction)
%   p: (n x 2), nodes of the triangulation, one node per row
%   dx: size of the domain in x-direction
%   t: (m3 x 3), triangulation, each row defines one triangle by indexing into p
%
% (C) 2019 by O. Junge, see COPYRIGHT 

n = size(p,1);            
P = [p; p+[dx 0]; p-[dx 0]];    % copy vertices to the left and the right
p_map = [1:n, 1:n, 1:n];        % map from extended nodes to original ones

t0 = delaunay(P(:,1),P(:,2));   % construct triangulation
% find triangles with no vertex in left copy and at least one in original domain
I = find(sum(t0 > 2*n,2) == 0 & sum(t0 < n,2) > 0);  
t = p_map(t0(I,:));             % remove rest from triangulation and map node numbers

