function [V,lambda,K,M] = solve_TO(mesh0,mesh1)

p = mesh0.p; tri = mesh0.t; pb = mesh0.pb; b = mesh0.b;
p1 = mesh1.p; tri1 = mesh1.t; pb1 = mesh1.pb; b1 = mesh1.b;

I = repmat([1 0; 0 1],[1, 1, size(tri,1)]);
[K,M] = assemble2(mesh0,I);     

I1 = repmat([1 0; 0 1],[1, 1, size(tri1,1)]);
[K1,~] = assemble2(mesh1,I1);

K = 0.5*(K+K1);

% if applicable, apply Dirichlet boundary condition   
if ~isempty(b)                        
    K(b,:) = 0; K(:,b) = 0; M(b,:) = 0; M(:,b) = 0;
    K(b,b) = speye(length(b),length(b));
end

% compute eigenvectors
[V,lambda] = eigs(K+1e-8*speye(size(K)),M,20,'SM'); 
[lambda,ord] = sort(diag(lambda),'descend'); 
V = V(:,ord);

% fix signs of eigenvectors
for i = 1:20
    V(:,i) = V(:,i)*sign(sum(V(:,i)));
end
