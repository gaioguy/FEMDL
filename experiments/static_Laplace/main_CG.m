addpath('../../src','../../src/2d'); clear all
at_tf = @(A) squeeze(A(:,:,end,:));

%% identity map 
t0 = 0; tf = 1; 
vf = @identity_vf;
DT  = @(x) at_tf(Dflow_map(vf,x,[t0 tf]));  % space derivative of flow map
DL = @(DT) 0.5*(eye(2) + inv(DT)*inv(DT)'); % dynamic Laplace
DLx = @(x) fapply1(DL, DT(x));              % evaluate DL at each row of x

%% triangulation: square
n = 100; x = linspace(0,1,n);
[X,Y] = meshgrid(x,x); p = [X(:) Y(:)];     % nodes
pb = [1:n^2; 1:n^2]';                       % nonperiodic boundary
t = delaunay(p);                            % t is m by 3 matrix of triangles

%% triangulation: circle
d = 2;
n = 50; phi = linspace(0,2*pi,n); 
m = 20; r = linspace(0.01,1,m); 
[R,Phi] = meshgrid(r,phi); 
X = R.*cos(Phi); 
Y = R.*sin(Phi);
p = [X(:) Y(:)];                            % nodes
pb = [1:n*m; 1:n*m]';                       % nonperiodic boundary
t = delaunay(p);                            % t is m by 3 matrix of triangles

%% triangulation: strip
d = 2;
n = 100; x = linspace(0,1,n); 
m = 20; y = linspace(0,0.1,m); 
[X,Y] = meshgrid(x,y); 
p = [X(:) Y(:)];                            % nodes
p = rand(n*m,2)*diag([1,0.1]);
pb = [1:n*m; 1:n*m]';                       % nonperiodic boundary
t = delaunay(p);                            % t is m by 3 matrix of triangles

%% assembly 
deg = 2;                                    % degree of quadrature
tic; A = triquad(p,t,DLx,deg); toc          % integrate DL on triangles
tic; [D,M] = assemble2(p,t,pb,A); toc       % assemble stiffness and mass matrices

%% solve eigenproblem
[V,L] = eigs(D,M,20,'SM');
[lam,order] = sort(diag(L),'descend'); V = V(:,order);

%% plot spectrum
figure(1); clf; plot(lam,'*'); 

%% plot eigenvector
figure(2); clf; plotf(p,t,pb,V(:,3),0);  colorbar
xlabel('$x$'); ylabel('$y$');


