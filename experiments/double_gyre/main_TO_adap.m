addpath('../../src/2d'); init

%% flow map
t0 = 0; tf = 1; nt = 6; tspan = linspace(t0,tf,nt); 
v = @double_gyre;
T = @(x) flowmap(v,x,tspan);

%% data points
nx = 50; ny = nx; n = nx*ny;
p(:,:,1) = grid2(nx,ny); pb = [1:n; 1:n]';

%% assembly (for missing data case)
p = T(p(:,:,1));                                % time integration
pm = 1;                                         % percentage of nodes to remove
%pm = 0.4;                                       % for missing data case
K = sparse(n,n); M = sparse(n,n); 
for k = 1:nt
    r = randperm(n,floor(pm*n))'; 
    [~,tr] = trimesh(p(r,:,k)); 
    I = repmat(eye(2),[1 1 size(tr,1)]);        % identity matrix
    [Kt,Mt] = assemble2(p(:,:,k),r(tr),pb,I); 
    K = K + Kt;  M = M + Mt; 
end;
% remove all zero rows and columns
S = sum(abs(K)); I = find(abs(S)>eps); K = K(I,I); M = M(I,I); 

%% solve eigenproblem
[V,L] = eigs(K,M,10,'SM'); 
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord); 

figure(1); plot(lam,'*');
[pI,tI,pbI] = trimesh(p(I,:,1));
figure(2), clf; plotf(pI,tI,pbI,normed(V(:,2)),0); colorbar

%% coherent partition
nc = 3;                                             % number of clusters
W = kmeans(V(:,1:nc),nc,'Replicates',20);           % kmeans clustering
figure(3); clf; scatter(p(:,1,1),p(:,2,1),30,W,'filled');
axis equal; axis tight; colormap(jet(nc));

%% advected partition
figure(4); clf; scatter(p(:,1,nt),p(:,2,nt),30,W,'filled'); 
axis equal; colormap(jet(nc));

