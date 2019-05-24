function F = triquad(m,f,deg)

%% TRIQUAD integrate function on 2d triangle mesh
%
% F = triquad(m,f,deg)
%   m: triangle mesh as produced by trimesh
%   f: function (with arbitrary values), must take (k x 2) matrix of points as input
%   deg: degree of quadrature rule
% return
%   F(i,:,...,:) is f integrated on triangle t(i,:)
%
% based on code from ifem by Long Chen
%
% (C) 2019 by O. Junge, see COPYRIGHT 

p = m.p; t = m.t;

x = p(t(1,1),:); 
tmp = f(x);
F = zeros([size(tmp) size(t,1)]);

[lam,w] = quadpts(deg); 
nq = size(lam,1);
for k = 1:nq
    x = lam(k,1)*p(t(:,1),:) + lam(k,2)*p(t(:,2),:) + lam(k,3)*p(t(:,3),:);
    F = F + w(k)*f(x);
end

