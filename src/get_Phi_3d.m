function Phi = get_Phi3(p,t,pb)

%% GET_PHI3 compute linear shape functions on elements
%
% Phi = get_Phi(p,t,pb) computes a cell array Phi of basis functions on each
% element
%   p: (n x 3), one node per row
%   t: (m x 4), integers, each row defines a triangle by indexing into p
%   pb: (n x 2), node pb(i,2) maps to pb(i,1) (for perodic boundaries)
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

m = size(t,1);                      % number of triangles
for e = 1:m                         % integration over each element 
    ns = t(e,1:4);                  % nodes of triangle e
    Pe = [ones(4,1),p(ns,:)];       % 4 x 4 matrix with rows = [1 xnode ynode znode]
    B = inv(Pe);                    % columns of C are coeffs in a+bx+cy+dz
    Phi{e} = B';
end
