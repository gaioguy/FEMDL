addpath('../../src'); clear all

%% flow map parameters
vf = @(t,x) rotating_double_gyre_vf(t,x);
t0 = 0; tf = 1; nt = 2; ts = linspace(t0,tf,nt); 
T = @(x) flow_map(vf,x,ts);
G = @(x) ones(size(x,1),3);

%% data points
n = 25*25; x = linspace(0,1,n);                            
[X,Y] = meshgrid(x,x); p0 = [X(:) Y(:)];
% p0 = rand(n,2); 
load p_Meiss
% p1 = rand(20000,2)*diag([1, 0.5])+ones(20000,2)*diag([0 0.5]);
% p2 = rand(200,2)*diag([1, 0.5]);
% p0 = [p1; p2]; n = size(p0,1);
pb = [1:n; 1:n]';

%% time integration
P = T(p{1});
for k = 1:nt, p{k} = [P(:,k) P(:,k+nt)]; end;

%% assembly
pm = 1;                                     % percentage of nodes to remove
D = sparse(n,n); M = sparse(n,n); 
for k = 1:nt
    r = randperm(n,floor(pm*n))'; 
    tr = delaunay(p{k}(r,:)); 
    t = [r(tr(:,1)), r(tr(:,2)), r(tr(:,3))];
    [Dt,Mt] = assemble_new(p{k},t,pb,ones(size(t,1),3));
    D = D + Dt;  M = M + Mt; 
end;

%% remove all zero rows and columns
S = sum(abs(D));
I = find(abs(S)>eps);
D = D(I,I); M = M(I,I); pI = p{1}(I,:); pbI = [1:size(pI,1); 1:size(pI,1)]';

%% solve eigenproblem
[V,L] = eigs(D,M,10,'SM'); 
[lam,order] = sort(diag(L),'descend'); V = V(:,order);

%% plot spectrum
figure(1); clf; plot(lam,'s','markerfacecolor','b'); axis tight
xlabel('$k$'); ylabel('$\lambda_k$');

%% plot eigenvector
figure(2), clf; tI = delaunayn(pI); 
plotev(tI,pI,pbI,V(:,3),1); axis([0 1 0 1]); colorbar
xlabel('$x$'); ylabel('$y$');

%% compute partition
n1 = 200; x1 = linspace(0,1,n1); [X1,Y1] = meshgrid(x1,x1); 
W = eval_p1(p{1},V(:,1:3),[X1(:) Y1(:)]);       % evaluate eigenvectors on grid
idx = kmeans(W, size(W,2));                     % kmeans clustering

%% plot partition
figure(3); clf; surf(X1,Y1,reshape(idx,n1,n1)); view(2); shading flat
axis equal; axis tight; xlabel('$x$'); ylabel('$y$'); colorbar
colormap(70*[ 1 1 1; 2 2 2; 3 3 3]/255);
