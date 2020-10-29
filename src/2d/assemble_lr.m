function L = assemble_lr(mesh,w)

%% ASSEMBLE_LR stiffness matrix for linear response
%
% L = ASSEMBLE_LR(mesh,w) computes the stiffness matrix L for the linear
% response
%   mesh: triangle mesh as produced by trimesh
%   w: (n x 2), each row is the value of w at the corresponding node
%
% based on code from ifem by Long Chen
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

p = mesh.p; t = mesh.t; pb = mesh.pb;
n = max(pb(:,2)); m = size(t,1);

[dphi,area] = gradbasis(p,t);

% assembly
L = sparse(n,n); 
for i = 1:3
    for j = 1:3
        for k = 1:3
            K = pb(t(:,k),2);   % global number of the k-th shape function
            Dij = area.*(dphi(:,1,i).*dphi(:,1,k) + dphi(:,2,i).*dphi(:,2,k)) ...
                      .*(w(K,1).*dphi(:,1,j) + w(K,2).*dphi(:,2,j));
            I = pb(t(:,i),2);
            J = pb(t(:,j),2);
            L = L + sparse(I,J,Dij,n,n);
        end
    end
end


