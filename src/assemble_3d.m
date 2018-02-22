function [D,M] = assemble_3d(p,t,pb,A)

%% ASSEMBLE_3D assemble stiffness and mass matrix 
%
% [D,M] = assemble_3d(p,t,pb,G) computes the stiffness matrix D and the mass matrix M
%   p: (n x 3), one node per row
%   t: (m x 4), integers, each row defines a triangle by indexing into p
%   pb: (n x 2), node pb(i,2) maps to pb(i,1) (for perodic boundaries)
%   A: (m x 6), row k is the space integral of some diffusion tensor
%               on the triangle t(k,:)
%
% based on code from ifem by Long Chen
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

n = max(pb(:,2)); 
m = size(t,1);

[dphi,volume] = gradbasis_3d(p,t);

D = sparse(n,n); 
M = sparse(n,n);
for i = 1:4
    for j = i:4
        Aij = -volume.*dotA3_sym(dphi(:,:,i),A,dphi(:,:,j));
        Mij = volume/12.*ones(m,1);
        I = pb(t(:,i),2); J = pb(t(:,j),2);
        if (j==i)
            D = D + sparse(I,J,Aij,n,n);
            M = M + sparse(I,J,Mij+volume/12,n,n);
        else
            D = D + sparse([I;J],[J;I],[Aij; Aij],n,n); 
            M = M + sparse([I;J],[J;I],[Mij; Mij],n,n);
        end        
    end
end
