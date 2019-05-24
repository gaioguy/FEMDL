function mesh = delaunay_C2(p0, dx)

%% DELAUNAY_C2 construct Delaunay triangulation on 2d cylinder
%
% mesh = delaunay_C2(p0, dx) constructs a Delaunay mesh on the 2d flat
%   cylinder (i.e. periodic in x-direction)
%   p0: (n x 2), nodes of the triangulation, one node per row
%   dx: size of the domain in x-direction
%   mesh: triangle mesh as produced by trimesh
%
% (C) 2019 by O. Junge, see COPYRIGHT 

n = size(p0,1);            
p1 = [p0; p0+[dx 0]; p0-[dx 0]];    % copy vertices to the left and the right
pb1 = [1:n, 1:n, 1:n];              % map from extended nodes to original ones

t1 = delaunay(p1(:,1),p1(:,2));     % construct triangulation on extended nodes
f1 = sum(t1 > 2*n,2) == 0;          % triangles with no vertex in left copy
f2 = sum(t1 <= n,2) > 0;            % triangles with at least one vertex in original domain
I = find(f1 & f2);                  % triangles with both properties
t2 = t1(I,:);                       % remove all triangles not in I from triangulation
J = unique(t2);                     % node numbers appearing in t1
Jinv(J) = find(J);                  % inverse of J as a map on node numbers
mesh.p = p1(J,:);                        % corresponding nodes
mesh.t = Jinv(t2);
mesh.pb = [1:length(J); pb1(J)]';        % and corresponding map for the boundary nodes
mesh.b = [];

% scatter(p1(:,1),p1(:,2),100,'filled','b'); hold on
% triplot(t1,p1(:,1),p1(:,2),'b');
% scatter(p(:,1),p(:,2),100,'filled','r'); 
% triplot(t,p(:,1),p(:,2),'r'); hold off

