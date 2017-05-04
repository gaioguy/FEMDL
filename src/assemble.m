function [D,M] = assemble(p,t,pb,G)

%% ASSEMBLE stiffness and mass matrices
%
% [D,M] = ASSEMBLE(p,t,pb,G) computes the stiffness matrix D and the mass matrix M
%   p: (n x 2), one node per row
%   t: (m x 3), integers, each row defines a triangle by indexing into p
%   pb: (n x 2), node pb(i,2) maps to pb(i,1) (for perodic boundaries)
%   G: (m x 3), each row defines a tensor on the correspondig triangle
%
% based on code from ifem by Long Chen
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

n = size(p,1); m = size(t,1);

% gradients of shape functions
v1 = p(t(:,3),:)-p(t(:,2),:);
v2 = p(t(:,1),:)-p(t(:,3),:);
v3 = p(t(:,2),:)-p(t(:,1),:);
area = 0.5*(-v3(:,1).*v2(:,2) + v3(:,2).*v2(:,1));       % areas of triangles
dphi(:,:,1) = [-v1(:,2)./(2*area), v1(:,1)./(2*area)];
dphi(:,:,2) = [-v2(:,2)./(2*area), v2(:,1)./(2*area)];
dphi(:,:,3) = [-v3(:,2)./(2*area), v3(:,1)./(2*area)];
area = abs(area);

% assembly
D = sparse(n,n); M = sparse(n,n);
for i = 1:3
    for j = i:3
        Aij = -area.*(dphi(:,1,i).*G(:,1).*dphi(:,1,j) ...
                    + dphi(:,1,i).*G(:,2).*dphi(:,2,j) ...
                    + dphi(:,2,i).*G(:,2).*dphi(:,1,j) ...
                    + dphi(:,2,i).*G(:,3).*dphi(:,2,j));
        Mij = area/12.*ones(size(dphi,1),1);
        I = pb(t(:,i),2); J = pb(t(:,j),2);
        if (j==i)
            D = D + sparse(I,J,Aij,n,n);
            M = M + sparse(I,J,Mij+area/12,n,n);
        else
            D = D + sparse([I;J],[J;I],[Aij; Aij],n,n);   
            M = M + sparse([I;J],[J;I],[Mij; Mij],n,n);
        end        
    end
end


