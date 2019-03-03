function [p1,t3,pb] = delaunay_T2(p, box)

%% DELAUNAY_T2 construct Delaunay triangulation on 2-torus
%
% [P1,t3,pb] = delaunay_T2(p, box) constructs a Delaunay
% triangulation on the (flat) 2-torus
%   p: (n x 2), nodes of the triangulation, one node per row
%   box: (1 x 4), specifies the domain (like in axis)
%   p1: (n1 x 2), extended set of nodes
%   t3: (m3 x 3), triangulation, each row def. one triangle by indexing
%   into p1
%   pb: (n1 x 2), boundary map: node pb(i,1) is identified with pb(:,2)
%
% (C) 2017 by O. Junge, see COPYRIGHT 

n = size(p,1);
xmin = box(1); xmax = box(2); ymin = box(3); ymax = box(4);

%% copy vertices around
dx = xmax-xmin; dy = ymax-ymin;
dp = [dx 0; dx dy; 0 dy; -dx dy; -dx 0; -dx -dy; 0 -dy; dx -dy];
P = [p; p+dp(1,:); p+dp(2,:); p+dp(3,:); p+dp(4,:); p+dp(5,:); p+dp(6,:); p+dp(7,:); p+dp(8,:);];
% clf; scatter(P(:,1),P(:,2),10,'b','filled'); hold on
PB = repmat(1:n,1,9);

%% construct triangulation and sort points
aS = alphaShape(P(:,1),P(:,2),1);
t = alphaTriangulation(aS);
% triplot(t,P(:,1),P(:,2),'b','linewidth',1); hold on; axis tight

%% extract triangles with at least one vertex in original domain
in = 1:n;           % vertices in original domain
% scatter(P(in,1),P(in,2),3,'r','filled'); 
tr = triangulation(t,P);
ti = vertexAttachments(tr,in');   % indices of triangles adjacent to each vertex in original domain
adj = [];
for i = 1:length(in), adj = [adj ti{i}]; end
adj = unique(adj);
t1 = t(adj,:);
% triplot(t1,P(:,1),P(:,2),'g');

%% remove duplicate triangles
flags = sum(t1 > 4*n,2) == 0;  % test if some vertex of a triangle is outside
t2 = t1(find(flags),:);
% triplot(t2,P(:,1),P(:,2),'m');

%% extract required points and corresponding triangles
all = unique(t2(:)); 
p1 = P(all,:);
Kt(all) = 1:length(all); Kt = Kt(:);        % mapping from old vertex numbers to 1:no_of_vertices_in_p1
t3 = [Kt(t2(:,1)) Kt(t2(:,2)) Kt(t2(:,3))];
% triplot(t3,P1(:,1),P1(:,2),'k'); hold on

%% compute boundary mapping
out = setdiff(all,in)'; 
in = setdiff(all,out)';
pb = ([1:length(all); in PB(out)])';
% scatter(P(out,1),P(out,2),10,'r','filled');  axis tight



