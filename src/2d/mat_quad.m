function F = mat_quad(mesh,f,deg)

%% MAT_QUAD integrate matrix valued function on triangles
%
% G = mat_quad(mesh,f,deg)
%   mesh; mesh data structure
%   f: a 2 x 2 matrix valued function, takes (k x 2) matrix of 2-vectors x as input
%   deg: degree of quadrature rule
% returns a (m x 2 x 2) tensor F, row i is f integrated on triangle t(i,:)
%
% based on code from ifem by Long Chen
%
% (C) 2018 by O. Junge and G. Froyland, see COPYRIGHT 

p = mesh.p; t = mesh.t; pb = mesh.pb; b = mesh.b;

m = size(t,1);
[lam,w] = quadpts(deg); 
nq = size(lam,1);
F = zeros(2,2,m);
for k = 1:nq
    x = lam(k,1)*p(t(:,1),:) + lam(k,2)*p(t(:,2),:) + lam(k,3)*p(t(:,3),:);
    F = F + w(k)*f(x);
end

