addpath('../../src'); clear all

%% identity map 
t0 = 0; tf = 1; nt = 2;
CG = @(x) cg_tensor(@identity_vf,x,linspace(t0,tf,nt));

%% triangulation
n = 50; x = linspace(0,1,n);
[X,Y] = meshgrid(x,x); p = [X(:) Y(:)];     % nodes
pb = [1:n^2; 1:n^2]';               % nonperiodic boundary
t = delaunay(p);                    % t is m by 3 matrix of triangles

%% assembly 
deg = 1;                            % degree of quadrature
G = compute_G(p,t,CG,deg);          % compute Cauchy-Green tensors
[D,M] = assemble(p,t,pb,G);     % assemble stiffness and mass matrices

%% solve eigenproblem
[V,L] = eigs(D,M,10,'SM');
[lam,order] = sort(diag(L),'descend'); V = V(:,order);

%% plot spectrum
figure(1); clf; plot(lam,'s','markerfacecolor','b'); axis tight
xlabel('$k$'); ylabel('$\lambda_k$');

%% plot eigenvector
figure(2); clf; plotev(t,p,pb,V(:,6),1); axis([0 1 0 1]); colorbar
xlabel('$x$'); ylabel('$y$');

%% compute partition
n1 = 200; x1 = linspace(0,1,n1); [X1,Y1] = meshgrid(x1,x1); 
V1 = eval_p1(p,V(:,1:3),[X1(:) Y1(:)]);     % evaluate eigenvectors on grid
idx = kmeans(V1, size(V1,2),'Replicates',10);           % kmeans clustering

%% plot partition
figure(3); clf; surf(X1,Y1,reshape(idx,n1,n1)); view(2); shading flat
axis equal; axis tight; xlabel('$x$'); ylabel('$y$'); colorbar
colormap(70*[ 1 1 1; 2 2 2; 3 3 3]/255);

