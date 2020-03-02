function [V,lam,K,M] = solve_CG(mesh, G, deg)

A = triquad(mesh,G,deg);            % integrate inverse Cauchy Green tensor
[K,M] = assemble2(mesh,A);          % assemble matrices

b = mesh.b;
if ~isempty(b)                      % apply Dirichlet boundary condition              
    K(b,:) = 0; K(:,b) = 0; M(b,:) = 0; M(:,b) = 0;
    K(b,b) = speye(length(b),length(b));
end

nev = 20;                           % compute eigenvectors
[V,lam] = eigs(K+1e-4*speye(size(K)),M,nev,'SM'); 
[lam,ord] = sort(diag(lam),'descend'); 
V = V(:,ord);

for i = 1:nev                       % fix sign 
    V(:,i) = V(:,i)*sign(sum(V(:,i)));
end
