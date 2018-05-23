function [dphi,vol] = gradbasis_3d(p,t)

% GRADBASIS_3D compute gradients of shape functions
%
% [dphi,area] = gradbasis_3d(p,t) computes the gradients of the shape
% functions on the standard simplex and the areas of the triangles
%   p: (n x 3), one node per row
%   t: (m x 4), integers, each row defines a triangle by indexing into p
%
% based on code from ifem by Long Chen
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

m = size(t,1);
face = uint32([t(:,[2 4 3]);t(:,[1 3 4]);t(:,[1 4 2]);t(:,[1 2 3])]);
v12 = p(face(:,2),:)-p(face(:,1),:);
v13 = p(face(:,3),:)-p(face(:,1),:);
nrml = cross(v12,v13,2);
v12 = v12(3*m+1:4*m,:); 
v13 = v13(3*m+1:4*m,:);
v14 = p(t(:,4),:)-p(t(:,1),:);
vol = dot(cross(v12,v13,2),v14,2);
vvol = vol*[1,1,1];
dphi = zeros(m,3,4);
dphi(:,:,1) = nrml(1:m,:)./vvol;
dphi(:,:,2) = nrml(m+1:2*m,:)./vvol;
dphi(:,:,3) = nrml(2*m+1:3*m,:)./vvol;
dphi(:,:,4) = nrml(3*m+1:4*m,:)./vvol;
vol = abs(vol/6);
