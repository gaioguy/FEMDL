addpath('../../src/2d'); init

%% flow map 
t0 = 0; tf = 1; 
v = @double_gyre;
T  = @(x) at_tf(flowmap(v,x,[t0 tf]));
T1 = @(x) at_tf(flowmap(v,x,[tf t0]));

%% triangulation
nx = 25; ny = nx; n = nx*ny; 
p = grid2(nx,ny); pb = [1:n; 1:n]'; b = []; 
tri = delaunayTriangulation(p);
t = tri.ConnectivityList; 

%% transfer operator approximation
Phi = get_Phi(p,t,pb);                      % Phis map elements to standard simplex
Alpha = get_alpha(tri,T1,Phi);              % alpha matrix approximates transfer op.

%% assembly
I = repmat(eye(2),[1 1 size(t,1)]);         % id tensor
[K,M] = assemble2(p,t,pb,I);       % assembly of stiffness and mass

%% eigenproblem
DL = 0.5*(K + Alpha'*K*Alpha); 
[V,L] = eigs(DL,M,10,'SM'); 
[lam,order] = sort(diag(L),'descend'); lam; V = V(:,order);

figure(1); plot(lam,'*'); axis tight, axis square
figure(2); plotf(p,t,pb,normed(V(:,2)),1); colorbar; 

%% coherent partition
nc = 3;                                         % number of clusters
W = kmeans(V(:,1:nc),nc,'Replicates',20);       % kmeans clustering
figure(3); scatter(p(:,1),p(:,2),30,W,'filled');
axis equal; axis tight; colormap(jet(nc));

%% advected partition
Tp = T(p);                                      % advect grid
figure(4); scatter(Tp(:,1),Tp(:,2),30,W,'filled'); 
axis equal; colormap(jet(nc));


