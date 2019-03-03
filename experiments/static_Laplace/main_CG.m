addpath('../../src','../../src/2d'); clear all

%% identity map 
t0 = 0; tf = 1; nt = 2;
CG = @(x) inv_CG(@identity_vf,x,linspace(t0,tf,nt));

%% triangulation: square
n = 100; x = linspace(0,1,n);
[X,Y] = meshgrid(x,x); p = [X(:) Y(:)];     % nodes
pb = [1:n^2; 1:n^2]';               % nonperiodic boundary
t = delaunay(p);                    % t is m by 3 matrix of triangles

%% triangulation: circle
d = 2;
n = 50; phi = linspace(0,2*pi,n); 
m = 20; r = linspace(0.01,1,m); 
[R,Phi] = meshgrid(r,phi); 
X = R.*cos(Phi); 
Y = R.*sin(Phi);
p = [X(:) Y(:)];     % nodes
pb = [1:n*m; 1:n*m]';               % nonperiodic boundary
t = delaunay(p);                    % t is m by 3 matrix of triangles

%% triangulation: strip
d = 2;
n = 100; x = linspace(0,1,n); 
m = 20; y = linspace(0,0.1,m); 
[X,Y] = meshgrid(x,y); 
p = [X(:) Y(:)];     % nodes
p = rand(n*m,2)*diag([1,0.1]);
pb = [1:n*m; 1:n*m]';               % nonperiodic boundary
t = delaunay(p);                    % t is m by 3 matrix of triangles

%% assembly 
deg = 2;                            % degree of quadrature
G = inv_CG_quad(p,t,CG,deg);          % compute Cauchy-Green tensors
[D,M] = assemble(p,t,pb,G);         % assemble stiffness and mass matrices

%% solve eigenproblem
[V,L] = eigs(D,M,20,'SM');
[lam,order] = sort(diag(L),'descend'); V = V(:,order);

%% plot spectrum
figure(1); clf; plot(2:r,lam(2:end),'*-'); 

%% plot eigenvector
figure(2); clf; plotf(p,t,pb,V(:,3),0);  colorbar
xlabel('$x$'); ylabel('$y$');


