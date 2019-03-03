function [D,M,ell,w,area] = assemble(p,t,pb,G,a)

%% ASSEMBLE stiffness and mass matrices
%
% [D,M] = ASSEMBLE(p,t,pb,G) computes the stiffness matrix D and the mass matrix M
%   p: (n x 2), one node per row
%   t: (m x 3), integers, each row defines a triangle by indexing into p
%   pb: (n x 2), node pb(i,2) maps to pb(i,1) (for perodic boundaries)
%   G: (m x 3), each row defines a tensor on the correspondig triangle
%   a: (n x 1), weighting of nodes
%
% based on code from ifem by Long Chen
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

n = max(pb(:,2));                   % number of nodes
m = size(t,1);                      % number of elements (triangles)

[dphi,area] = gradbasis(p,t);

lt = repmat(area,1,3);
ell = accumarray(pb(t(:),2),lt(:),[n 1])/3;      % l(i) = \int \phi_i^t dl
al = a./ell;
w = zeros(m,1); 
for i = 1:3                                      % weights 
    I = pb(t(:,i),2); 
    w = w + 1/3*al(I);          
end            

D = sparse(n,n); M = sparse(n,n);
for k = 1:3
    for l = 1:3
        D_kl = -area.*dotA2(dphi(:,:,k),G,dphi(:,:,l)).*w; 
        K = pb(t(:,k),2); 
        L = pb(t(:,l),2);  
        if (k==l)
            j = setdiff([1 2 3],k); 
            I = pb(t(:,j(1)),2); 
            J = pb(t(:,j(2)),2);
            M_kl = area.*2/60.*(3*al(K) + al(I) + al(J));
        else
            j = setdiff([1 2 3],[k l]); 
            J = pb(t(:,j),2);
            M_kl = area.*2/120.*(2*al(K) + 2*al(L) + al(J)); 
        end
        D = D + sparse(K,L,D_kl,n,n);
        M = M + sparse(K,L,M_kl,n,n);
    end
end


