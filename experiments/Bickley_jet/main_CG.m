addpath('../../src'); clear all

%% flow map 
t0 = 0; tf = 60*60*24*40; nt = 10;
CG = @(x) cg_tensor(@bickleyjet_vf,x,linspace(t0,tf,nt));

%% triangulation
nx = floor(200);  ny = nx/20*6;  n = nx*ny;
[X,Y] = meshgrid(linspace(0,20,nx),linspace(-3,3,ny)); p = [X(:) Y(:)];
pb = [1:n; [1:((nx-1)*ny), 1:ny]]';             % boundary periodic in x
t = delaunay(p); 

%% assembly
deg = 1;                                        % degree of quadrature
G = compute_G(p,t,CG,deg);
[D,M] = assemble(p,t,pb,G);

%% solve eigenproblem
[V,L] = eigs(D+1e-4*speye(size(D)),M,20,'SM'); 
[lam,order] = sort(diag(L),'descend'); V = V(:,order);

%% plot spectrum
figure(1); clf; plot(lam,'s','markerfacecolor','b'); axis tight
xlabel('$k$'); ylabel('$\lambda_k$')

%% plot eigenvector
figure(2), clf; plotev(t,p,pb,V(:,8),0); colorbar
xlabel('$x$'); ylabel('$y$'); 

%% compute partition
nx1 = 200; ny1 = nx1/20*6; x1 = linspace(0,20,nx1); y1 = linspace(-3,3,ny1);
[X1,Y1] = meshgrid(x1,y1); 
V1 = eval_p1(p,V(pb(:,2),1:8),[X1(:) Y1(:)]);       % evaluate eigenvectors on grid
idx = kmeans(V1, size(V1,2),'Replicates',10);       % kmeans clustering

%% plot partition
figure(3); clf; surf(X1,Y1,reshape(idx,ny1,nx1)); view(2); shading flat
axis equal; axis tight; xlabel('$x$'); ylabel('$y$'); colorbar


