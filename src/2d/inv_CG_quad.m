function G = inv_CG_quad2(p,t,CG,deg)

%% INV_CG_QUAD integrate Cauchy-Green tensor on triangles
%
% G = compute_G(p,t,CG,deg)
%   p: (n x 2), one node per row
%   t: (m x 3), integers, each row defines a triangle by indexing into p
%   CG: a 2 x 2 matrix valued function
%   deg: degree of quadrature rule
% returns a (m x 3) matrix G with G(i,:) = [G(1,1) G(1,2) G(2,2)]
% and G the averaged inverse Cauch-Green tensor on triangle t(i,:)
%
% based on code from ifem by Long Chen
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

m = size(t,1);
[lam,w] = quadpts(deg); 
nq = size(lam,1);
G = zeros(m,3);
for k = 1:nq
    x = lam(k,1)*p(t(:,1),:) + lam(k,2)*p(t(:,2),:) + lam(k,3)*p(t(:,3),:);
    G = G + w(k)*CG(x);
end

