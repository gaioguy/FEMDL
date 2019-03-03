function Alpha = get_Alpha(tri,T1,Phi)

%% GET_ALPHA compute transfer operator matrix
%
% Alpha = get_Alpha(tri,T1,Phi) computes the matrix Alpha of coefficients
% of the transfer operator for the map T with inverse T1
%   tri: triangulation
%   T1: inverse of the underlying map
%   Phi: given basis, cell array of (3 x 3) matrices
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

p = tri.Points;
t = tri.ConnectivityList;
n = size(p,1);
Alpha = sparse(n,n);
T1p = T1(p);
Te = pointLocation(tri,T1p);    % numbers of triangles the points T1p are located in

for k = 1:size(p,1)             % loop over all points
    Tns = t(Te(k),:);           % node number of triangle Te(k)
    Alpha(k,Tns) = (Phi{Te(k)}*[1 T1p(k,:)]')';
end
