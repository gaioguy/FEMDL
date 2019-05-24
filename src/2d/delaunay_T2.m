function mesh = delaunay_T2(p0, dx, dy)

%% DELAUNAY_T2 construct Delaunay triangulation on 2-torus
%
% mesh = delaunay_T2(p0, dx) constructs a Delaunay mesh on the 2d flat
%   torus (i.e. periodic in x- and y-direction)
%   p0: (n x 2), nodes of the triangulation, one node per row
%   dx,dy: size of the domain in x- and y-direction
%   mesh: triangle mesh as produced by trimesh
%
% (C) 2019 by O. Junge, see COPYRIGHT 

n = size(p0,1);

%% copy vertices around
dp = [dx 0; dx dy; 0 dy; -dx dy; -dx 0; -dx -dy; 0 -dy; dx -dy];
P = [p0; p0+dp(1,:); p0+dp(2,:); p0+dp(3,:); p0+dp(4,:); p0+dp(5,:); p0+dp(6,:); p0+dp(7,:); p0+dp(8,:);];
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
mesh.p = P(all,:);
Kt(all) = 1:length(all); Kt = Kt(:);        % mapping from old vertex numbers to 1:no_of_vertices_in_p1
mesh.t = [Kt(t2(:,1)) Kt(t2(:,2)) Kt(t2(:,3))];
% triplot(t3,P1(:,1),P1(:,2),'k'); hold on

%% compute boundary mapping
out = setdiff(all,in)'; 
in = setdiff(all,out)';
mesh.pb = ([1:length(all); in PB(out)])';
% scatter(P(out,1),P(out,2),10,'r','filled');  axis tight

mesh.b = [];


