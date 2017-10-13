addpath('../../src'); clear all

%% identity map 
t0 = 0; tf = 1; nt = 2;
CG = @(x) cg_tensor(@identity_vf,x,linspace(t0,tf,nt));

%% triangulation
n = 100; x = linspace(0,1,n);
[X,Y] = meshgrid(x,x); p = [X(:) Y(:)];     % nodes
pb = [1:n^2; 1:n^2]';               % nonperiodic boundary
t = delaunay(p);                    % t is m by 3 matrix of triangles

%% assembly 
deg = 1;                            % degree of quadrature
G = compute_G(p,t,CG,deg);          % compute Cauchy-Green tensors
[D,M] = assemble(p,t,pb,G);         % assemble stiffness and mass matrices

%% solve eigenproblem
[V,L] = eigs(D,M,10,'SM');
[lam,order] = sort(diag(L),'descend'); lam/pi^2, V = V(:,order);

%% plot eigenvector
figure(2); clf; plotev(t,p,pb,V(:,6),1); axis([0 1 0 1]); colorbar
xlabel('$x$'); ylabel('$y$');

