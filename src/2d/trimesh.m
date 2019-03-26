function [p,t,pb,b,tri] = trimesh(p)

n = size(p,1);
pb = [1:n; 1:n]';
tri = alphaShape(p(:,1),p(:,2),1);
t = alphaTriangulation(tri);
b = unique(boundaryFacets(tri));
