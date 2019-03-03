function [dphi,area] = gradbasis(p,t)

% [dphi,area] = GRADBASIS(p,t) computes the gradients of the shape
% functions on the standard simplex and the areas of the triangles
%   p: (n x 2), one node per row
%   t: (m x 3), integers, each row defines a triangle by indexing into p
%
% based on code from ifem by Long Chen
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

v(:,:,1) = p(t(:,3),:)-p(t(:,2),:);
v(:,:,2) = p(t(:,1),:)-p(t(:,3),:);
v(:,:,3) = p(t(:,2),:)-p(t(:,1),:);
area = 0.5*(-v(:,1,3).*v(:,2,2) + v(:,2,3).*v(:,1,2));
dphi(:,:,:) = [-v(:,2,:)./(2*area), v(:,1,:)./(2*area)];
area = abs(area);
