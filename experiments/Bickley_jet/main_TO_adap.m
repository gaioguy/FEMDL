addpath('../../src','../../src/2d'); clear all; clc; colormap jet
xmod20 = @(x) [mod(x(:,1,:),20) x(:,2,:)];

%% flow map 
t0 = 0; days = 60*60*24; tf = 40*days; nt = 11;  % number of time steps for Laplacian
T = @(x) xmod20(flowmap(@bickleyjet,x,linspace(t0,tf,nt)));
T1 = @(x) flow_map(@bickleyjet,x,linspace(t0,tf,nt));

%% triangulation
dim = 2; % dimension of state space
nx = 120;  ny = nx/20*6;  n = nx*ny;
[xi,yi] = meshgrid(linspace(0,20.01,nx),linspace(-3,3,ny));
p(:,:,1) = [xi(:) yi(:)]; 

%% time integration
tic; p = T(p(:,:,1)); toc

%% assembly
tic; pm = 1;  % pm = 0.2                         % percentage of nodes to remove                             
D = sparse(n,n); M = sparse(n,n);
for k = 1:nt
    r = randperm(n,floor(pm*n))';                % draw random sample of nodes in p{k}
    [pr,tr,pbr] = delaunay_C2(p(r,:,k),20);
    I = kron([1 0 1],ones(size(tr,1),1));
    [Dr, Mr] = assemble(pr,tr,pbr,I);
    [Ir,Jr,Sr] = find(Dr);  D = D + sparse(r(Ir),r(Jr),Sr,n,n);
    [Ir,Jr,Sr] = find(Mr);  M = M + sparse(r(Ir),r(Jr),Sr,n,n);
end; toc;

%% remove all zero rows and columns
S = sum(abs(D)); I = find(abs(S)>eps); 
D = D(I,I); M = M(I,I); 

%% eigenproblem
tic; [V,L] = eigs(D,M,15,'SM');  toc
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord); 

%% plot spectrum
figure(1); clf; plot(lam,'*'); axis tight, axis square
xlabel('$k$'); ylabel('$\lambda_k$')

%% plot eigenvector
figure(2); ev_no = 9;
pI = p(I,:,1); tI = delaunay(pI); pbI = [1:length(I);1:length(I)]';
plotf(pI,tI,pbI,normed(V(:,ev_no)),0); caxis([-1,1]); colorbar
xlabel('$x$'); ylabel('$y$'); 

%% compute partition
nx1 = 400; ny1 = nx1/20*6; 
[X1,Y1] = meshgrid(linspace(0,20,nx1),linspace(-3,3,ny1)); 
nc = 9;  
V1 = eval_p1(pI,V(:,1:nc),[X1(:) Y1(:)]);  % evaluate eigenvectors on grid
idx = kmeans(V1,nc,'Replicates',20);       % kmeans clustering

%% plot partition
figure(3); subplot(313);  clf; 
surf(X1,Y1,reshape(idx,ny1,nx1)); view(2); shading flat; colorbar
axis equal; axis tight; ylabel('$y$'); xlabel('$x$');



