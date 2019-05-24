# FEMDL

A MATLAB package for a finite element based discretization of dynamic Laplacians

## Usage
A typical script which computes coherent sets for a rotating double gyre is
```Matlab
addpath('FEMDL/src/2d'); init

t0 = 0; tf = 1;                                     % initial and final time
v = @double_gyre;                                   % vector field
T = @(x) at_tf(flowmap(v, x, [t0 tf]));             % flow map
DL = @(DT) 0.5*(eye(2) + inv(DT)*inv(DT)');         % dynamic Laplacian
DLx = @(x) fapply1(DL, D(T,x));                     % evaluate DL at each row of x
p = grid2(50,50);                                   % nodes on a grid
m = trimesh(p);                                     % triangulation
A = triquad(m,DLx,2);                               % integrate DL on triangles
[K,M] = assemble2(m,A);                             % assemble stiffness and mass matrix
[V,L] = eigs(K,M,6,'SM');                           % solve eigenproblem
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord);  % sort eigenvalues
figure; plotf2(m,normed(V(:,2)),0); colorbar;       % plot eigenvector
W = kmeans(V(:,1:3),3,'Replicates',20);             % kmeans clustering
figure; scatter(p(:,1),p(:,2),30,W,'filled');       % plot clustering
Tp = T(p);                                          % advect grid
figure; scatter(Tp(:,1),Tp(:,2),30,W,'filled');     % plot advected clustering
colormap(jet(nc));
```
There are more demo computations provided in experiments/ which should be
executed after changing into the corresponding folder.

## Reference
Froyland, G.; Junge, O.: Robust FEM-based extraction of finite-time coherent sets using scattered, sparse, and incomplete trajectories, SIAM J. Appl. Dyn. Syst. 17 (2), p. 1891-1924, 2018. https://arxiv.org/abs/1705.03640


