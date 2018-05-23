function [p1,t3,pb] = delaunay_T3(p, box)

%% DELAUNAY_T3 construct Delaunay triangulation on 3-torus
%
% [P1,t3,pb] = delaunay_T3(p, box) constructs a Delaunay
% triangulation on the (flat) 3-torus
%   p: (n x 3), nodes of the triangulation, one node per row
%   box: (1 x 6), specifies the domain (like in axis)
%   p1: (n1 x 3), extended set of nodes
%   t3: (m3 x 4), triangulation, each row defines one triangle by indexing
%   into p1
%   pb: (n1 x 2), boundary map: node pb(i,1) is identified with pb(:,2)
%
% (C) 2017 by O. Junge, see COPYRIGHT 

n = size(p,1);
xmin = box(1); xmax = box(2); ymin = box(3); ymax = box(4); zmin = box(5); zmax = box(6);

%% copy vertices around
dx = xmax-xmin; dy = ymax-ymin; dz = zmax-zmin;
dp = [ 0 0 0; dx 0 0; 0 dy 0; 0 0 dz; dx dy 0; dx 0 dz; 0 dy dz; dx dy dz; ...
       0 0 -dz; dx 0 -dz; dx dy -dz; 0 dy -dz; -dx dy -dz; -dx 0 -dz; -dx -dy -dz; 0 -dy -dz; dx -dy -dz; ...
       -dx dy 0; -dx 0 0; -dx -dy 0; 0 -dy 0; dx -dy 0; ...
       -dx dy dz; -dx 0 dz; -dx -dy dz; 0 -dy dz; dx -dy dz];  

%dp = [dx 0 0; dx dy 0; 0 dy 0; -dx dy 0; -dx 0 0; -dx -dy 0; 0 -dy 0; dx -dy 0; ...
%      dx 0 dz; dx dy dz; 0 dy dz; -dx dy dz; -dx 0 dz; -dx -dy dz; 0 -dy dz; dx -dy dz; ...
%      dx 0 -dz; dx dy -dz; 0 dy -dz; -dx dy -dz; -dx 0 -dz; -dx -dy -dz; 0 -dy -dz; dx -dy -dz; ...
%      0 0 dz; 0 0 -dz];    
P = []; for k = 1:27, P = [P; p+dp(k,:)]; end        
%clf; scatter3(P(:,1),P(:,2),P(:,3),10,'b','filled'); hold on
PB = repmat(1:n,1,27);

%% construct triangulation and sort points
% aS = alphaShape(P(:,1),P(:,2),P(:,3),1);
% t = alphaTriangulation(aS); 
tr = delaunayTriangulation(P);
t = tr.ConnectivityList;
% clf; tetramesh(tr,'FaceAlpha',0); hold on; axis tight

%% extract triangles with at least one vertex in original domain
in = 1:n;           % vertices in original domain
% scatter3(P(in,1),P(in,2),P(in,3),12,'r'); 

%% tr = triangulation(t,P);
ti = vertexAttachments(tr,in');   % indices of triangles adjacent to each vertex in original domain
adj = [];
for i = 1:length(in), adj = [adj ti{i}]; end
adj = unique(adj);
t1 = t(adj,:);
%tetramesh(t1,P);

%% remove duplicate triangles
flags = sum(t1 > 8*n,2) == 0;  % test if some vertex of a triangle is outside
t2 = t1(find(flags),:);
% tetramesh(t2,P,'FaceColor','red');

%% extract required points and corresponding triangles
all = unique(t2(:)); 
p1 = P(all,:);
Kt(all) = 1:length(all); Kt = Kt(:);        % mapping from old vertex numbers to 1:no_of_vertices_in_p1
t3 = [Kt(t2(:,1)) Kt(t2(:,2)) Kt(t2(:,3)) Kt(t2(:,4))];
% tetramesh(t3,p1,'FaceColor','cyan');

%% compute boundary mapping
out = setdiff(all,in)'; 
in = setdiff(all,out)';
pb = ([1:length(all); in PB(out)])';
% scatter(P(out,1),P(out,2),10,'r','filled');  axis tight



