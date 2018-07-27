addpath('../../src/');  clear all

%% identity map 
t0 = 0; tf = 1; nt = 2;
CG = @(x) inv_CG(@double_gyre,x,linspace(t0,tf,nt));

%% triangulation
n = 15; x = linspace(0,1,n);
[X,Y] = meshgrid(x,x); p = [X(:) Y(:)];     % nodes
pb = [1:n^2; 1:n^2]';                       % nonperiodic boundary
t = delaunay(p);                            % t is m by 3 matrix of triangles

%% assembly 
deg = 1;                                    % degree of quadrature
G = inv_CG_quad(p,t,CG,deg);
[D,M] = assemble_P2(p,t,pb,G);                % assemble stiffness and mass matrices

%% solve eigenproblem
[V,L] = eigs(D,M,6,'SM');
[lam,order] = sort(real(diag(L)),'descend'); lam, V = V(:,order);

%% plot eigenvector
n1 = 200; x1 = linspace(0,1,n1); [X1,Y1] = meshgrid(x1,x1); 
V1 = eval_P2(p,t,V(:,1:4),[X1(:) Y1(:)]);
figure(2); clf; surf(X1,Y1, reshape(V1(:,2),n1,n1)); shading flat
axis([0 1 0 1]); colorbar; xlabel('$x$'); ylabel('$y$');

%% compute partition
idx = kmeans(V1, size(V1,2),'Replicates',10);           % kmeans clustering

%% plot partition
figure(3); clf; surf(X1,Y1,reshape(idx,n1,n1)); view(2); shading flat
axis equal; axis tight; xlabel('$x$'); ylabel('$y$'); colorbar
colormap(50*[ 1 1 1; 2 2 2; 3 3 3; 4 4 4]/255);

