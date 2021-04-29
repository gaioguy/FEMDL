function [K,M] = assemble2(mesh,G)

%% ASSEMBLE stiffness and mass matrices
%
% [D,M] = ASSEMBLE(mesh,G) computes the stiffness matrix D and the mass matrix M
%   mesh: triangle mesh as produced by trimesh
%   G: (2 x 2 x m), where G(:,:,i) defines a symmetric 2x2 matrix on the 
%      correspondig triangle
%
% based on code from ifem by Long Chen
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

p = mesh.p; t = mesh.t; pb = mesh.pb;
np = max(pb(:,2)); nt = size(t,1);
[dphi,area] = gradbasis(p,t);
if length(size(G))==2 & size(G)==[2,2]
    G = repmat(G, [1,1,size(mesh.t,1)]);
end
PG = permute(G,[3 1 2]);
K = sparse(np,np); M = sparse(np,np);
for i = 1:3
    for j = i:3
        Kij = -area.*(dphi(:,1,i).*PG(:,1,1).*dphi(:,1,j) ...
                    + dphi(:,1,i).*PG(:,1,2).*dphi(:,2,j) ...
                    + dphi(:,2,i).*PG(:,2,1).*dphi(:,1,j) ...
                    + dphi(:,2,i).*PG(:,2,2).*dphi(:,2,j));
        Mij = area/12.*ones(size(dphi,1),1);
        I = pb(t(:,i),2); 
        J = pb(t(:,j),2);
        if (j==i)
            K = K + sparse(I,J,Kij,np,np);
            M = M + sparse(I,J,Mij+area/12,np,np);
        else
            K = K + sparse([I;J],[J;I],[Kij; Kij],np,np);   
            M = M + sparse([I;J],[J;I],[Mij; Mij],np,np);
        end        
    end
end


