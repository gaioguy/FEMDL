addpath('../../src'); clear all; clc; colormap jet
extract = @(A,I) A(:,I);

%% flow map 
t0 = 0; tf = 1; nt = 2; tspan = linspace(tf,t0,nt);
T1 = @(x) extract(flow_map(@double_gyre,x,tspan),[2,4]);

%% nodes 
n = 25; x = linspace(0,1,n);
[X,Y] = meshgrid(x,x); p = [X(:) Y(:)];
%p = rand(n,2)*diag([1,1]); 
pb = [1:n^2; 1:n^2]';

%% triangulation
tri = delaunayTriangulation(p);
t = tri.ConnectivityList; 
m = size(t,1);

%% transfer operator approximation
Phi = get_Phi(p,t,pb);                      % Phis map elements to standard simplex
tic; Alpha = get_alpha(tri,T1,Phi); toc     % alpha matrix approximates transfer op.

%% assembly
A = kron([1 0 1],ones(size(t,1),1));       % 2 x 2 identity matrix
tic; [D,M] = assemble(p,t,pb,A); toc       % stiffness and mass matrix

%% eigenproblem
DL = 0.5*(D + Alpha'*D*Alpha); 
tic, [V,L] = eigs(DL,M,10,'SM'); toc
[lam,order] = sort(diag(L),'descend'); lam; V = V(:,order);

%% plot spectrum
figure(1); plot(lam,'*'); axis tight, axis square
xlabel('$k$'); ylabel('$\lambda_k$');

%% plot eigenvector
figure(2); clf; w = V(:,2); 
plotf(p,t,pb,w/norm(w),1); 
axis([0 1 0 1]); colorbar; 
xlabel('$x$'); ylabel('$y$');

%% compute partition
n1 = 200; x1 = linspace(0,1,n1); 
[X1,Y1] = meshgrid(x1,x1); 
V1 = eval_p1(p,V(:,1:3),[X1(:) Y1(:)]);     % evaluate eigenvectors on grid
idx = kmeans(V1, size(V1,2));               % kmeans clustering

%% plot partition
figure(3); clf; surf(X1,Y1,reshape(idx,n1,n1)); view(2); shading flat
axis equal; axis tight; xlabel('$x$'); ylabel('$y$'); colorbar
colormap(70*[ 1 1 1; 2 2 2; 3 3 3]/255);


