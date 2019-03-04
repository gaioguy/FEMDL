function fx = eval_P2(p,t,f_P2,x)

%% eval_P2 evaluate P2 Lagrange functions
%
% EVAL_P2(t,p,pb,v,mesh) plots the scalar function given by 
%   p: (n x 2), one node per row
%   t: (m x 3), integers, each row defines a triangle by indexing into p
%   f_P2: (n x k), values of k functions on the P2 dofs
%   x: (l x 2) l points in which the functions will be evaluated
%   fx: (l x k) values of f_P2 at x
%
% (C) 2018 by O. Junge, see COPYRIGHT 

elem2dof = dofP2(t);
tr = triangulation(t,p);
id = pointLocation(tr,x); 
%I = find(isnan(id)), x(I,:), id(I)
%id = pointLocation(tr,x(I,:)) 

lam = cartesianToBarycentric(tr,id,x);

phi(:,1) = lam(:,1).*(2*lam(:,1)-1);
phi(:,2) = lam(:,2).*(2*lam(:,2)-1);
phi(:,3) = lam(:,3).*(2*lam(:,3)-1);
phi(:,4) = 4*lam(:,2).*lam(:,3);
phi(:,5) = 4*lam(:,3).*lam(:,1);
phi(:,6) = 4*lam(:,1).*lam(:,2);

dof = elem2dof(id,:);

fx = zeros(size(x,1),size(f_P2,2));
for j = 1:size(f_P2,2)
    f = f_P2(:,j);
    f_dof = f(dof); 
    fx(:,j) = sum(phi.*f_dof,2);
end
