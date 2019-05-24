addpath('../../src/2d'); init

%% flow map etc.
t0 = 0; tf = 1; 
v = @double_gyre;
T = @(x) at_tf(flowmap(v, x, [t0 tf]));             % flow map
DL = @(DT) 0.5*(eye(2) + inv(DT)*inv(DT)');         % dynamic Laplacian
DLx = @(x) fapply1(DL, D(T,x));                     % evaluate DL at each row of x

%% triangulation
n = 50; 
p = grid2(n,n);                                     % nodes
m = trimesh(p);                                     % triangulation

%% assembly 
deg = 2;                                            % degree of quadrature
A = triquad(m,DLx,deg);                             % integrate DL on triangles
[K,M] = assemble2(m,A);                             % assemble stiffness and mass matrices

%% eigenproblem
[V,L] = eigs(K,M,6,'SM'); 
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord);

figure(1); clf; plot(lam,'*'); 
figure(2); clf; plotf(m,normed(V(:,2)),0); colorbar; 

%% coherent partition
nc = 3;                                             % number of clusters
W = kmeans(V(:,1:nc),nc,'Replicates',20);           % kmeans clustering
figure(3); clf; scatter(p(:,1),p(:,2),30,W,'filled');
axis equal; axis tight; colormap(jet(nc));

%% advected partition
Tp = T(p);                                          % advect grid
figure(4); clf; scatter(Tp(:,1),Tp(:,2),30,W,'filled'); 
axis equal; axis tight; colormap(jet(nc));

