# FEMDL

A MATLAB package for a finite element based discretization of dynamic Laplacians

## Usage
A typical script which computes coherent sets for a rotating double gyre is
```Matlab
addpath('../../src/2d'); init

t0 = 0; tf = 1;                                     % initial and final time
v = @double_gyre;                                   % vector field
T = @(x) at_tf(flowmap(v, x, [t0 tf]));             % flow map
DL = @(DT) 0.5*(eye(2) + inv(DT)*inv(DT)');         % dynamic Laplacian
DLx = @(x) fapply1(DL, D(T,x));                     % evaluate DL at each row of x

nx = 50; ny = nx;                                   % number of grid points
[p,t,pb,b] = trimesh(grid2(nx,ny));                 % triangulation

deg = 2;                                            % degree of quadrature
A = triquad(p,t,DLx,deg);                           % integrate DL on triangles
[K,M] = assemble2(p,t,pb,A);                        % assemble stiffness and mass matrix

[V,L] = eigs(K,M,6,'SM');                           % solve eigenproblem
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord);
figure(1); plot(lam,'*');                           % plot spectrum
figure(2); plotf(p,t,pb,normed(V(:,2)),0);colorbar; % plot eigenvector 

nc = 3;                                             % number of clusters
W = kmeans(V(:,1:nc),nc,'Replicates',20);           % kmeans clustering
figure(3); scatter(p(:,1),p(:,2),30,W,'filled');    % plot clustering
axis equal; axis tight; colormap(jet(nc));

Tp = T(p);                                          % advect grid
figure(4); scatter(Tp(:,1),Tp(:,2),30,W,'filled');  % plot advected clustering
axis equal; axis tight; colormap(jet(nc));
```
There are more demo computations provided in experiments/ which should be
executed after changing into the corresponding folder.

## Reference
Froyland, G.; Junge, O.: Robust FEM-based extraction of finite-time coherent sets using scattered, sparse, and incomplete trajectories, SIAM J. Appl. Dyn. Syst. 17 (2), p. 1891-1924, 2018. https://arxiv.org/abs/1705.03640


