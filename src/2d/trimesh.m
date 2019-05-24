function [mesh,tri] = trimesh(p)

n = size(p,1);
tri = alphaShape(p(:,1),p(:,2),1);
mesh.t = alphaTriangulation(tri);
mesh.b = unique(boundaryFacets(tri));
mesh.p = p;
mesh.pb = [1:n; 1:n]';
