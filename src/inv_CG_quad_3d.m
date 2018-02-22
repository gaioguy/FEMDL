function K = invCG_quad_3d(p,t,CG,deg)

%% G = invCG_quad_3d(p,t,CG,deg) quadrature of time averaged inverse Cauchy-Green 
%% tensor on triangles
%   p: (n x 3), one node per row
%   t: (m x 4), integers, each row defines a triangle by indexing into p
%   CG: function handle, computes the time averaged inverse Cauchy-Green tensor
%   deg: degree of quadrature rule
% returns a (n x 6) matrix G with G(i,:) = [CG(1,1) CG(1,2) CG(1,3) CG(2,2) CG(2,3) CG(3,3)]
% and CG the integral of the time averaged inverse Cauch-Green tensor on triangle t(i,:)
%
% based on code from ifem by Long Chen
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

m = size(t,1);
[lam,w] = quadpts_3d(deg);
nq = size(lam,1);
K = zeros(m,6);
for k = 1:nq
	x = lam(k,1)*p(t(:,1),:) ...
      + lam(k,2)*p(t(:,2),:) ...
      + lam(k,3)*p(t(:,3),:) ...
      + lam(k,4)*p(t(:,4),:);
    K = K + w(k)*CG(x);           
end

