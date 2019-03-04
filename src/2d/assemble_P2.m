function [D,M] = assemble_P2(node,t,pb,G)

%% ASSEMBLE_P2 stiffness and mass matrices for P2 Lagrange elements
%
% [D,M] = ASSEMBLE_P2(p,t,pb,G) computes the stiffness matrix D and the mass matrix M
%   p: (n x 2), one node per row
%   t: (m x 3), integers, each row defines a triangle by indexing into p
%   pb: (n x 2), node pb(i,2) maps to pb(i,1) (for perodic boundaries)
%   G: (m x 3), each row defines a symmetric 2x2 matrix on the correspondig triangle
%
% based on code from ifem by Long Chen
%
% (C) 2018 by O. Junge, see COPYRIGHT 

[elem2dof,edge,bdDof] = dofP2(t);

n = size(node,1);  
m = size(t,1); 
ne = size(edge,1); 
ndof = n + ne;

[Dlam,area] = gradbasis(node,t);

% generate sparse pattern
ii = zeros(21*m,1); 
jj = zeros(21*m,1); 
I = 0;
for i = 1:6
    for j = i:6
        ii(I+1:I+m) = double(elem2dof(:,i)); 
        jj(I+1:I+m) = double(elem2dof(:,j));  
        I = I + m;
    end
end

%% quadrature points
quadorder = 2;          % default order
[lam, w] = quadpts(quadorder);
nq = size(lam,1);       % number of quad points

sD = zeros(21*m,nq);
sM = zeros(21*m,nq);

% phi at quadrature points
phi(:,1) = lam(:,1).*(2*lam(:,1)-1);
phi(:,2) = lam(:,2).*(2*lam(:,2)-1);
phi(:,3) = lam(:,3).*(2*lam(:,3)-1);
phi(:,4) = 4*lam(:,2).*lam(:,3);
phi(:,5) = 4*lam(:,3).*lam(:,1);
phi(:,6) = 4*lam(:,1).*lam(:,2);

PG = permute(G,[3 1 2]);
for p = 1:nq
    % Dphip at quadrature points
    Dphip(:,:,1) = (4*lam(p,1)-1).*Dlam(:,:,1);            
    Dphip(:,:,2) = (4*lam(p,2)-1).*Dlam(:,:,2);            
    Dphip(:,:,3) = (4*lam(p,3)-1).*Dlam(:,:,3);            
    Dphip(:,:,4) = 4*(lam(p,2)*Dlam(:,:,3)+lam(p,3)*Dlam(:,:,2));
    Dphip(:,:,5) = 4*(lam(p,3)*Dlam(:,:,1)+lam(p,1)*Dlam(:,:,3));
    Dphip(:,:,6) = 4*(lam(p,1)*Dlam(:,:,2)+lam(p,2)*Dlam(:,:,1));
    I = 0;
    for i = 1:6
        for j = i:6
            Dij = 0; Mij = 0;
%            Dij = Dij + w(p)*dot(Dphipp(:,:,i),Dphipp(:,:,j),2).*pde.d(pxy);
            % Dij = Dij + w(p)*dot(Dphipp(:,:,i),Dphipp(:,:,j),2);
            Dij = Dij + w(p)*(Dphip(:,1,i).*PG(:,1,1).*Dphip(:,1,j) ...
                            + Dphip(:,1,i).*PG(:,1,2).*Dphip(:,2,j) ...
                            + Dphip(:,2,i).*PG(:,2,1).*Dphip(:,1,j) ...
                            + Dphip(:,2,i).*PG(:,2,2).*Dphip(:,2,j));
            Mij = Mij + w(p)*phi(p,i).*phi(p,j);
            Dij = Dij.*area; Mij = Mij.*area;
            sD(I+1:I+m,p) = -Dij;
            sM(I+1:I+m,p) = Mij;
            I = I + m;
        end
    end
end
sD = sum(sD,2); sM = sum(sM,2);

dI = (ii == jj);   
uI = ~dI;
D = sparse(ii(dI),jj(dI),sD(dI),ndof,ndof);
DU = sparse(ii(uI),jj(uI),sD(uI),ndof,ndof);
D = D + DU + DU';
M = sparse(ii(dI),jj(dI),sM(dI),ndof,ndof);
MU = sparse(ii(uI),jj(uI),sM(uI),ndof,ndof);
M = M + MU + MU';


