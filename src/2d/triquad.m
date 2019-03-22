function F = triquad(p,t,f,deg)

%% TRIQUAD integrate function on triangles
%
% F = triquad(p,t,f,deg)
%   p: (n x 2), nodes of the triangle mesh, one node per row
%   t: (m x 3), triangles, each row defines a triangle by indexing into p
%   f: function (possibly tensor valued), must take (k x 2) matrix of points as input
%   deg: degree of quadrature rule
%   F: F(i,:,...,:) is f integrated on triangle t(i,:)
%
% based on code from ifem by Long Chen
%
% (C) 2019 by O. Junge, see COPYRIGHT 

m = size(t,1);
[lam,w] = quadpts(deg); 
x = p(t(1,1),:); tmp = f(x);
F = zeros([size(tmp) m]);
nq = size(lam,1);
for k = 1:nq
    x = lam(k,1)*p(t(:,1),:) + lam(k,2)*p(t(:,2),:) + lam(k,3)*p(t(:,3),:);
    F = F + w(k)*f(x);
end

