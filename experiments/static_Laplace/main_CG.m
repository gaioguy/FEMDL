addpath('../../src'); clear all

%% identity map 
t0 = 0; tf = 1; nt = 2;
CG = @(x) inv_CG(@identity_vf,x,linspace(t0,tf,nt));

%% triangulation: square
n = 100; x = linspace(0,1,n);
[X,Y] = meshgrid(x,x); p = [X(:) Y(:)];     % nodes
pb = [1:n^2; 1:n^2]';               % nonperiodic boundary
t = delaunay(p);                    % t is m by 3 matrix of triangles

%% triangulation: circle
n = 50; phi = linspace(0,2*pi,n); 
m = 20; r = linspace(0.01,1,m); 
[R,Phi] = meshgrid(r,phi); 
X = R.*cos(Phi); 
Y = R.*sin(Phi);
p = [X(:) Y(:)];     % nodes
pb = [1:n*m; 1:n*m]';               % nonperiodic boundary
t = delaunay(p);                    % t is m by 3 matrix of triangles

%% assembly 
deg = 1;                            % degree of quadrature
G = inv_CG_quad(p,t,CG,deg);          % compute Cauchy-Green tensors
[D,M] = assemble(p,t,pb,G);         % assemble stiffness and mass matrices

%% solve eigenproblem
[V,L] = eigs(D+1e-8*speye(size(D)),M,10,'SM');
[lam,order] = sort(diag(L),'descend'); lam/pi^2, V = V(:,order);

%% plot eigenvector
figure(2); clf; plotf(p,t,pb,V(:,1),1);  colorbar
xlabel('$x$'); ylabel('$y$');
view(30,40);

