function Phi = get_Phi(p,t,pb)

%% GET_PHI compute linear shape functions on elements
%
% Phi = get_Phi(p,t,pb) computes a cell array Phi of basis functions on each
% element
%   p: (n x 2), one node per row
%   t: (m x 3), integers, each row defines a triangle by indexing into p
%   pb: (n x 2), node pb(i,2) maps to pb(i,1) (for perodic boundaries)
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

m = size(t,1);                      % number of triangles
for e = 1:m                         % integration over each element 
    ns = t(e,1:3);                  % nodes of triangle e
    Pe = [ones(3,1),p(ns,:)];       % 3 x 3 matrix with rows = [1 xnode ynode]
    B = inv(Pe);                    % columns of C are coeffs in a+bx+cy 
    Phi{e} = B';
end
