addpath('../../src'); clear all; colormap jet

%% flow map
t0 = 0; tf = 1; nt = 6; tspan = linspace(t0,tf,nt); 
T = @(x) flow_map(@double_gyre,x,tspan);

%% data points
n = 25; x = linspace(0,1,n);                            
[X,Y] = meshgrid(x,x); p{1} = [X(:) Y(:)];
% p{1} = rand(n^2,2); 
% load p_Meiss
% p1 = rand(20000,2)*diag([1, 0.5])+ones(20000,2)*diag([0 0.5]);
% p2 = rand(200,2)*diag([1, 0.5]);
% p{1} = [p1; p2]; n = size(p{1},1);
pb = [1:n^2; 1:n^2]';

%% time integration
tic; P = T(p{1}); toc
for k = 1:nt, p{k} = [P(:,k) P(:,k+nt)]; end;

%% assembly (for missing data case)
tic; pm = 1;                                % percentage of nodes to remove
% pm = 0.4;                                   % for missing data case
D = sparse(n^2,n^2); M = sparse(n^2,n^2); 
for k = 1:nt
    r = randperm(n^2,floor(pm*n^2))'; 
    tr = delaunay(p{k}(r,:)); 
    t = [r(tr(:,1)), r(tr(:,2)), r(tr(:,3))];
    A = kron([1 0 1],ones(size(t,1),1));      % 2 x 2 identity matrix
    [Dt{k},Mt{k}] = assemble(p{k},t,pb,A); 
    D = D + Dt{k};  M = M + Mt{k}; 
end;
toc

%% remove all zero rows and columns
S = sum(abs(D));
I = find(abs(S)>eps);
D = D(I,I); M = M(I,I); pI = p{1}(I,:); pbI = [1:size(pI,1); 1:size(pI,1)]';

%% solve eigenproblem
I = speye(size(D));
tic; [V,L] = eigs(D,M,10,'SM'); toc
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord); lam

%% plot spectrum
figure(1); plot(lam,'*'); axis tight, axis square
xlabel('$k$'); ylabel('$\lambda_k$');

%% plot eigenvector
figure(2), clf; tI = delaunayn(pI); 
w = V(:,3);
plotf(pI,tI,pbI,w/norm(w),1); axis([0 1 0 1]); colorbar
xlabel('$x$'); ylabel('$y$');

%% compute partition
n1 = 200; x1 = linspace(0,1,n1); [X1,Y1] = meshgrid(x1,x1); 
W = eval_p1(p{1},V(:,1:3),[X1(:) Y1(:)]);       % evaluate eigenvectors on grid
idx = kmeans(W, size(W,2));                     % kmeans clustering

%% plot partition
figure(3); clf; surf(X1,Y1,reshape(idx,n1,n1)); view(2); shading flat
axis equal; axis tight; xlabel('$x$'); ylabel('$y$'); colorbar
colormap(70*[ 1 1 1; 2 2 2; 3 3 3]/255);

