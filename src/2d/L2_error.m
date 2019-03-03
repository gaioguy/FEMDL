function error = L2_error(p,t,u,uh)
  
m = size(t,1);   
quadOrder = 3; 

%% compute errors element-wise 
error = zeros(m,1); 
[lambda,w] = quadpts(quadOrder);
phi = lambda;                             % basis function at quadrature points
nQuad = size(lambda,1);
for j = 1:nQuad  
    uhp = uh(t(:,1))*phi(j,1) + ...       % evaluate uh at quadrature point
          uh(t(:,2))*phi(j,2) + ...
          uh(t(:,3))*phi(j,3); 
    pxy = lambda(j,1)*p(t(:,1),:) ...     % quadrature points in the x-y coordinate
        + lambda(j,2)*p(t(:,2),:) ...
        + lambda(j,3)*p(t(:,3),:);
    error = error + w(j)*(u(pxy) - uhp).^2;
end
[~,area] = gradbasis(p,t);                % compute area of triangles
error = sqrt(sum(area.*error));
