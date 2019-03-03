function [D,M] = assemble2(mesh,G)

%% ASSEMBLE stiffness and mass matrices
%
% [D,M] = ASSEMBLE(p,t,pb,G) computes the stiffness matrix D and the mass matrix M
%   p: (n x 2), one node per row
%   t: (m x 3), integers, each row defines a triangle by indexing into p
%   pb: (n x 2), node pb(i,2) maps to pb(i,1) (for perodic boundaries)
%   G: (2 x 2 x m), where G(:,:,i) defines a symmetric 2x2 matrix on the 
%      correspondig triangle
%
% based on code from ifem by Long Chen
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

p = mesh.p; t = mesh.t; pb = mesh.pb; b = mesh.b;
n = max(pb(:,2)); m = size(t,1);
[dphi,area] = gradbasis(p,t);

% assembly
PG = permute(G,[3 1 2]);
D = sparse(n,n); M = sparse(n,n);
for i = 1:3
    for j = i:3
        Dij = -area.*(dphi(:,1,i).*PG(:,1,1).*dphi(:,1,j) ...
                    + dphi(:,1,i).*PG(:,1,2).*dphi(:,2,j) ...
                    + dphi(:,2,i).*PG(:,2,1).*dphi(:,1,j) ...
                    + dphi(:,2,i).*PG(:,2,2).*dphi(:,2,j));
        Mij = area/12.*ones(size(dphi,1),1);
        I = pb(t(:,i),2); J = pb(t(:,j),2);
        if (j==i)
            D = D + sparse(I,J,Dij,n,n);
            M = M + sparse(I,J,Mij+area/12,n,n);
        else
            D = D + sparse([I;J],[J;I],[Dij; Dij],n,n);   
            M = M + sparse([I;J],[J;I],[Mij; Mij],n,n);
        end        
    end
end


