addpath('../../src/2d'); init 

%% flow map etc.
t0 = 0; tf = 1; 
v = @(t,x) zeros(size(x));
T  = @(x) at_tf(flowmap(v,x,[t0 tf]));      % flow map
DL = @(DT) 0.5*(eye(2) + inv(DT)*inv(DT)'); % dynamic Laplace
DLx = @(x) fapply1(DL, D(T,x));             % evaluate DL at each row of x

%% triangulation: square
nx = 20; ny = nx;                           % number of grid points
mesh = trimesh(grid2(nx,ny));         % triangular mesh

%% assembly 
deg = 1;                                    % degree of quadrature
A = triquad(mesh,DLx,deg);                   % integrate DL on triangles
[K,M] = assemble2(mesh,A);                % assemble stiffness and mass matrices

%% solve eigenproblem
[V,L] = eigs(K,M,10,'SM');
[lam,order] = sort(diag(L),'descend'); V = V(:,order);

figure(1); plot(lam,'*'); 
figure(2); plotf(mesh,V(:,2),0);  colorbar


