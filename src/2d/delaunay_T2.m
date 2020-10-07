function mesh = delaunay_T2(p0, dx, dy)

%% DELAUNAY_T2 construct Delaunay triangulation on 2d torus
%
% mesh = delaunay_T2(p0, dx, dy) constructs a Delaunay mesh on the 2d flat
%   torus [0,dx] x [0,dy] (i.e. periodic in x- and y-direction)
%   p0: (n x 2), nodes of the triangulation, one node per row
%   dx, dy: size of the domain in x- and y-direction
%   mesh: triangle mesh as produced by trimesh
%
% (C) 2020 by O. Junge, see COPYRIGHT 

n = size(p0,1);           
p = [p0;           p0+[dx 0];   p0+[0 dy]; ...
     p0+[dx dy];   p0+[-dx dy]; p0+[-dx 0]; ...
     p0+[-dx -dy]; p0+[0 -dy];  p0+[dx -dy]];    % copy vertices
pb = repmat(1:n,1,9);               % map from extended nodes to original ones

Dt = delaunayTriangulation(p);      % construct triangulation on extended nodes
t = Dt.ConnectivityList;        
t = t(find(sum(t <= n,2) > 0),:);   % triangles with at least one vertex in original domain
[~,I] = unique(sort(pb(t),2),'rows'); % remove duplicate triangles
t = t(I,:);                         % resulting reduced triangulation
J = unique(t);                      % used nodes
Jinv(J) = find(J);                  % inverse of J as a map on node numbers

mesh.p = p(J,:);                    % corresponding nodes
mesh.t = Jinv(t);
mesh.pb = [1:length(J); pb(J)]';    % and corresponding map for the boundary nodes
mesh.b = [];

% scatter(p1(:,1),p1(:,2),10,'filled','b'); hold on
% plot([0 dx dx 0 0], [0 0 dy dy 0],'k','linewidth',2);
% triplot(t1,p1(:,1),p1(:,2),'b');
% triplot(mesh.t,mesh.p(:,1),mesh.p(:,2),'r','linewidth',2); 
% scatter(p0(:,1),p0(:,2),40,'g','linewidth',2); 
% hold off; axis tight

