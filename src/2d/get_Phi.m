function Phi = get_Phi(m)

%% GET_PHI compute linear shape functions on elements
%
% Phi = get_Phi(m) computes a cell array Phi of basis functions on each
% element
%   m: triangle mesh as produced by trimesh
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

p = m.p; t = m.t; pb = m.pb;
nt = size(t,1);                         % number of triangles
for e = 1:nt                            % integration over each element 
    ns = t(e,1:3);                      % nodes of triangle e
    Pe = [ones(3,1),p(ns,:)];           % 3 x 3 matrix with rows = [1 xnode ynode]
    B = inv(Pe);                        % columns of C are coeffs in a+bx+cy 
    Phi{e} = B';
end
