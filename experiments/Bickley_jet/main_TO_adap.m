addpath('../../src'); clear all; clc; colormap jet

%% flow map 
t0 = 0; days = 60*60*24; tf = 40*days; nt = 11;  % number of time steps for Laplacian
T = @(x) flow_map(@bickleyjet,x,linspace(t0,tf,nt));

%% triangulation
dim = 2; % dimension of state space
nx = 100;  ny = nx/20*6;  n = nx*ny;
[xi,yi] = meshgrid(linspace(0,20.01,nx),linspace(-3,3,ny));
p0 = [xi(:) yi(:)];  

%% time integration
tic; P = T(p0); toc
for k = 1:nt, p{k} = [mod(P(:,k),20) P(:,k+nt)]; end;

%% assembly
tic; 
pm = 1;                               % percentage of nodes to keep
pm = 0.2;                               
D = sparse(n,n); M = sparse(n,n);
for k = 1:nt
    [p1,~,pb] = delaunay_C2(p{k},[0,20]);
    n1 = size(p1,1);
    % r1 = randperm(n1,floor(pm*n1))';
    r1 = randsample(n1-n:n,floor(pm*n))';
    tr = delaunay(p1(r1,:));
    t = r1(tr);
    I = kron([1 0 1],ones(size(t,1),1));
    [Dt,Mt] = assemble(p1,t,pb,I);
    D = D + Dt; M = M + Mt; 
end; toc

%% remove all zero rows and columns
S = sum(abs(D)); I = find(abs(S)>eps); 
D = D(I,I); M = M(I,I); 
[pI,tI,pbI] = delaunay_C2(p{1}(I,:),[0,20]);

%% eigenproblem
tic; [V,L] = eigs(D,M,15,'SM');  toc
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord); 

%% plot spectrum
figure(1); clf; plot(lam,'*'); axis tight, axis square
xlabel('$k$'); ylabel('$\lambda_k$')

%% plot eigenvector
figure(2); ev_no = 4;
plotf(pI,tI,pbI,normed(V(:,ev_no)),0); caxis([-1,1]); colorbar
xlabel('$x$'); ylabel('$y$'); 

%% compute partition
nx1 = 400; ny1 = nx1/20*6; x1 = linspace(0,20,nx1); y1 = linspace(-3,3,ny1);
[X1,Y1] = meshgrid(x1,y1); 
nc = 8; W(pbI(:,1),:) = V(pbI(:,2),:); 
V1 = eval_p1(pI,W(:,1:nc),[X1(:) Y1(:)]);  % evaluate eigenvectors on grid
idx = kmeans(V1,nc,'Replicates',10);       % kmeans clustering

%% plot partition
figure(3); subplot(313);  clf; 
surf(X1,Y1,reshape(idx,ny1,nx1)); view(2); shading flat; colorbar
axis equal; axis tight; ylabel('$y$'); xlabel('$x$');



