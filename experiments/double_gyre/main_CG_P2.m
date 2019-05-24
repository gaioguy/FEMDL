addpath('../../src/2d'); init

%% rotating double gyre map 
t0 = 0; tf = 1; 
vf = @double_gyre;
T  = @(x) at_tf(flowmap(vf, x, [t0 tf]));   % flow map
DL = @(DT) 0.5*(eye(2) + inv(DT)*inv(DT)'); % dynamic Laplacian
DLx = @(x) fapply1(DL, D(T,x));             % evaluate DL at each row of x

%% triangulation
n = 25; 
p = grid2(n,n);                                     % nodes
m = trimesh(p);                                     % triangulation

%% assembly 
deg = 2;                                            % degree of quadrature
A = triquad(m,DLx,deg);                             % integrate DL on triangles
[K,M] = assemble_P2(m,A);                           % assemble stiffness and mass matrices

%% solve eigenproblem
[V,L] = eigs(K,M,6,'SM');
[lam,order] = sort(real(diag(L)),'descend'); V = V(:,order);

%% plot eigenvector
p1 = grid2(100,100); 
m1 = trimesh(p1);                              
V1 = eval_P2(m,V(:,1:4),m1.p);
figure(1); clf; plotf(m1,normed(-V1(:,3)),0); colorbar

%% coherent partition
nc = 3;                                             % number of clusters
W = kmeans(V1(:,1:nc),nc,'Replicates',20);           % kmeans clustering
figure(2); clf; scatter(p1(:,1),p1(:,2),30,W,'filled');
axis equal; axis tight; colormap(jet(nc));

%% advected partition
Tp1 = T(p1);                                          % advect grid
figure(3); clf; scatter(Tp1(:,1),Tp1(:,2),30,W,'filled'); 
axis equal; axis tight; colormap(jet(nc));

