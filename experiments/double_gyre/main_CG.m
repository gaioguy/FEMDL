addpath('../../src/2d'); clear all; colormap jet; 
at_tf = @(A) squeeze(A(:,:,end,:));

%% rotating double gyre map 
t0 = 0; tf = 1; 
vf = @double_gyre;
DT  = @(x) at_tf(Dflow_map(vf,x,[t0 tf]));  % space derivative of flow map
DL = @(DT) 0.5*(eye(2) + inv(DT)*inv(DT)'); % dynamic Laplace
DLx = @(x) fapply1(DL, DT(x));              % evaluate DL at each row of x

%% triangulation
n = 25; x = linspace(0,1,n);
[X,Y] = meshgrid(x,x); p = [X(:) Y(:)];     % nodes
% p = rand(n^2,2);                          % for Fig. 3
pb = [1:n^2; 1:n^2]';                       % nonperiodic boundary
t = delaunay(p);                            % t is m by 3 matrix of triangles

%% assembly 
deg = 2;                                    % degree of quadrature
tic; A = triquad(p,t,DLx,deg); toc          % integrate DL on triangles
tic; [D,M] = assemble2(p,t,pb,A); toc       % assemble stiffness and mass matrices

%% solve eigenproblem
tic; [V,L] = eigs(D,M,6,'SM'); toc
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord);

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
tic; V1 = eval_p1(p,V(:,1:4),[X1(:) Y1(:)]);  toc   % evaluate eigenvectors on grid
tic; idx = kmeans(V1, size(V1,2),'Replicates',10);  toc         % kmeans clustering

%% plot partition
figure(3); clf; surf(X1,Y1,reshape(idx,n1,n1)); view(2); shading flat
axis equal; axis tight; xlabel('$x$'); ylabel('$y$'); colorbar
colormap(50*[ 1 1 1; 2 2 2; 3 3 3; 4 4 4]/255);

