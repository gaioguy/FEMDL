function mesh = delaunay_T2(p0, dx, dy)

%% DELAUNAY_T2 construct Delaunay triangulation on 2d torus
%
% mesh = delaunay_T2(p0, dx, dy) constructs a Delaunay mesh on the 2d flat
%   torus [0,dx] x [0,dy] (i.e. periodic in x- and y-direction)
%   p0: (n x 2), nodes of the triangulation, one node per row
%   dx, dy: size of the domain in x- and y-direction
%   mesh: triangle mesh as produced by trimesh
%
% (C) 2019 by O. Junge, see COPYRIGHT 

n = size(p0,1);            
p1 = [p0;           p0+[dx 0];   p0+[0 dy]; ...
      p0+[dx dy];   p0+[-dx dy]; p0+[-dx 0]; ...
      p0+[-dx -dy]; p0+[0 -dy];  p0+[dx -dy]];    % copy vertices
pb1 = repmat(1:n,1,9);              % map from extended nodes to original ones

Dt = delaunayTriangulation(p1);     % construct triangulation on extended nodes
t1 = Dt.ConnectivityList; 
f1 = sum(t1 > 4*n,2) == 0;          % triangles with no vertex in left or lower copies
f2 = sum(t1 <= n,2) > 0;            % triangles with at least one vertex in original domain
I = find(f1 & f2);                  % triangles with both properties
d = 1e-15;
crn = [dx-d dy-d];
crnI = pointLocation(Dt,crn);
I = unique([I;crnI]);
t2 = t1(I,:);                       % remove all triangles not in I from triangulation
J = unique(t2);                     % node numbers appearing in t2
Jinv(J) = find(J);                  % inverse of J as a map on node numbers
mesh.p = p1(J,:);                   % corresponding nodes
mesh.t = Jinv(t2);
mesh.pb = [1:length(J); pb1(J)]';   % and corresponding map for the boundary nodes
mesh.b = [];

% scatter(p1(:,1),p1(:,2),10,'filled','b'); hold on
% plot([0 dx dx 0 0], [0 0 dy dy 0],'k','linewidth',2);
% triplot(t1,p1(:,1),p1(:,2),'b');
% scatter(mesh.p(:,1),mesh.p(:,2),10,'filled','r'); 
% scatter(p0(:,1),p0(:,2),30,'k','linewidth',2); 
% triplot(mesh.t,mesh.p(:,1),mesh.p(:,2),'r','linewidth',2); hold off

