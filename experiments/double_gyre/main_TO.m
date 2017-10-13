addpath('../../src'); clear all
extract = @(A,I) A(:,I);

%% flow map parameters
v = @(t,x) rotating_double_gyre_vf(t,x);
t0 = 0; tf = 1; nt = 2;  % number of time steps for Laplacian
T1 = @(x) extract(flow_map(v,x,linspace(tf,t0,nt)),[2,4]);

%% triangulation
n = 25; x = linspace(0,1,n);
[X,Y] = meshgrid(x,x); p = [X(:) Y(:)];
%p = rand(n,2)*diag([1,1]); 
pb = [1:n^2; 1:n^2]';
tri = delaunayTriangulation(p);
t = tri.ConnectivityList; m = size(t,1);

%% assembly and transfer operator approximation
Phi = get_Phi(p,t,pb);              % Phis map elements to standard simplex
Alpha = get_alpha(tri,T1,Phi);      % alpha matrix approximates transfer op.
CG = kron([1 0 0 1],ones(size(t,1),1));      % 2 x 2 identity matrix
[D,M] = assemble(p,t,pb,CG);        % stiffness and mass matrix
B = 0.5*(D + Alpha'*D*Alpha); 

%% solve eigenproblem
tic, [V,L] = eigs(B,M,10,'SM'); toc
[lam,order] = sort(diag(L),'descend'); V = V(:,order);

%% plot spectrum
figure(1); clf; plot(lam,'s','markerfacecolor','b'); axis tight
xlabel('$k$'); ylabel('$\lambda_k$');

%% plot eigenvector
figure(2), clf; plotev(t,p,pb,V(:,2),1); axis([0 1 0 1]); colorbar
xlabel('$x$'); ylabel('$y$');

%% compute partition
n1 = 200; x1 = linspace(0,1,n1); [X1,Y1] = meshgrid(x1,x1); 
V1 = eval_p1(p,V(:,1:3),[X1(:) Y1(:)]);     % evaluate eigenvectors on grid
idx = kmeans(V1, size(V1,2));               % kmeans clustering

%% plot partition
figure(3); clf; surf(X1,Y1,reshape(idx,n1,n1)); view(2); shading flat
axis equal; axis tight; xlabel('$x$'); ylabel('$y$'); colorbar
colormap(70*[ 1 1 1; 2 2 2; 3 3 3]/255);


