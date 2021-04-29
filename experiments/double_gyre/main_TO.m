addpath('../../src/2d'); clear all; init

%% flow map 
t0 = 0; tf = 1; 
v = @double_gyre;
T  = @(x) at_tf(flowmap(v,x,[t0 tf]));
T1 = @(x) at_tf(flowmap(v,x,[tf t0]));

%% triangulation
nx = 25; ny = nx; n = nx*ny; 
m.p = grid2(nx,ny); m.pb = [1:n; 1:n]'; m.b = []; 
tri = delaunayTriangulation(m.p);
m.t = tri.ConnectivityList; 

%% transfer operator approximation
Phi = get_Phi(m);                               % Phis map elements to standard simplex
Alpha = get_alpha(tri,T1,Phi);                  % alpha matrix approximates transfer op.

%% assembly
I = repmat(eye(2),[1 1 size(m.t,1)]);           % id tensor
[K,M] = assemble2(m,I);                         % assembly of stiffness and mass

%% eigenproblem
DL = 0.5*(K + Alpha'*K*Alpha); 
[V,L] = eigs(DL,M,10,'SM'); 
[lam,order] = sort(diag(L),'descend'); lam; V = V(:,order);

figure(1); plot(lam,'*'); axis tight, axis square
figure(2); plotf(m,normed(V(:,2)),1); colorbar; 

%% coherent partition
nc = 3;                                         % number of clusters
W = kmeans(V(:,1:nc),nc,'Replicates',20);       % kmeans clustering
figure(3); clf; scatter(m.p(:,1),m.p(:,2),30,W,'filled');
axis equal; axis tight; colormap(jet(nc));

%% advected partition
Tp = T(m.p);                                      % advect grid
figure(4); clf; scatter(Tp(:,1),Tp(:,2),30,W,'filled'); 
axis equal; colormap(jet(nc));


