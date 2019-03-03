addpath('../../src/2d'); clear all; colormap jet; 

%% rotating double gyre map 
t0 = 0; tf = 1; nt = 2; tspan = linspace(t0,tf,nt);
CG = @(x) inv_CG(@double_gyre,x,tspan);

%% triangulation
n = 25; x = linspace(0,1,n);
[X,Y] = meshgrid(x,x); p = [X(:) Y(:)];     % nodes
% p = rand(n^2,2);                          % for Fig. 3
pb = [1:n^2; 1:n^2]';                       % nonperiodic boundary
t = delaunay(p);                            % t is m by 3 matrix of triangles

%% assembly 
deg = 2;                                    % degree of quadrature
tic; G = inv_CG_quad(p,t,CG,deg);  toc      % compute inverse Cauchy-Green tensors
tic; [D,M] = assemble(p,t,pb,G);  toc       % assemble stiffness and mass matrices

%% solve eigenproblem
tic; [V,L] = eigs(D,M,6,'SM'); toc
[lam,ord] = sort(diag(L),'descend'); lam, V = V(:,ord);

%% plot spectrum
figure(1); plot(lam,'*'); axis tight, axis square
xlabel('$k$'); ylabel('$\lambda_k$');

%% plot eigenvector
figure(2); clf; w = V(:,2); 
plotf(p,t,pb,w/norm(w),1); 
axis([0 1 0 1]); colorbar; 
xlabel('$x$'); ylabel('$y$');

%% compute partition
n1 = 200; x1 = linspace(0,1,n1); 
[X1,Y1] = meshgrid(x1,x1); 
tic; V1 = eval_p1(p,V(:,1:4),[X1(:) Y1(:)]);  toc   % evaluate eigenvectors on grid
tic; idx = kmeans(V1, size(V1,2),'Replicates',10);  toc         % kmeans clustering

%% plot partition
figure(3); clf; surf(X1,Y1,reshape(idx,n1,n1)); view(2); shading flat
axis equal; axis tight; xlabel('$x$'); ylabel('$y$'); colorbar
colormap(50*[ 1 1 1; 2 2 2; 3 3 3; 4 4 4]/255);

%% plot CG tensor
figure(4); clf
trisurf(t,p(:,1),p(:,2),zeros(size(p,1),1),log10(abs(G(:,1)+G(:,3)))); shading flat
view(2); axis([0 1 0 1]); axis equal; axis tight; caxis([0,4]); colormap jet; 
c = logspace(0,4,5); cbar = colorbar('YTick',log10(c),'YTickLabel',c);


