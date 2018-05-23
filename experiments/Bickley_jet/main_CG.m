addpath('../../src'); clear all; 

%% flow map 
t0 = 0; days = 60*60*24; tf = 40*days; nt = 2; 
tspan = linspace(t0,tf,nt);
CG = @(x) inv_CG(@bickleyjet,x,tspan);

%% nodes
nx = 25;  ny = nx/20*6;  n = nx*ny;
[X,Y] = meshgrid(linspace(0,20,nx),linspace(-3,3,ny)); 
p = [X(:) Y(:)];
pb = [1:n; [1:((nx-1)*ny), 1:ny]]';             % boundary periodic in x

%% triangulation
t = delaunay(p); 

%% assembly
deg = 1;                                        % degree of quadrature
tic; G = inv_CG_quad(p,t,CG,deg); toc
tic; [D,M] = assemble(p,t,pb,G); toc

%% solve eigenproblem
tic; [V,L] = eigs(D,M,15,'SM'); toc
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord);

%% plot spectrum
figure(1); clf; plot(lam,'*'); axis tight, axis square
xlabel('$k$'); ylabel('$\lambda_k$')

%% plot eigenvector
figure(2), plotf(p,t,pb,normed(-V(:,3)),0); caxis([-1,1]); colorbar
ylabel('$y$');  % xlabel('$x$'); 

%% compute partition
nx1 = 400; ny1 = nx1/20*6; x1 = linspace(0,20,nx1); y1 = linspace(-3,3,ny1);
[X1,Y1] = meshgrid(x1,y1); 
nc = 8;
V1 = eval_p1(p,V(pb(:,2),1:nc),[X1(:) Y1(:)]);       % evaluate eigenvectors on grid
idx = kmeans(V1, size(V1,2),'Replicates',10);       % kmeans clustering

%% plot partition
%figure(3); clf; 
surf(X1,Y1,reshape(idx,ny1,nx1)); view(2); shading flat
axis equal; axis tight; xlabel('$x$'); ylabel('$y$'); colorbar


