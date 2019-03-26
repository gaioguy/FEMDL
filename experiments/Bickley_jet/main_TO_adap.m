addpath('../../src/2d'); init
xmod20 = @(x) [mod(x(:,1,:),20) x(:,2,:)];

%% flow map 
t0 = 0; days = 60*60*24; tf = 40*days; nt = 11;  % number of time steps for Laplacian
v = @bickleyjet;
T = @(x) xmod20(flowmap(v,x,linspace(t0,tf,nt)));

%% triangulation
dom = [0 -3; 20 3]; dx = diff(dom);             % domain 
nx = 100; ny = dx(2)/dx(1)*nx;  n = nx*ny;      % number of grid points
p(:,:,1) = grid2(nx,ny)*diag(dx) + dom(1,:);    % grid 

%% assembly
p = T(p(:,:,1));                                % time integration
pm = 1;  % pm = 0.2                             % uncomplete data case: percentage of nodes to remove                             
K = sparse(n,n); M = sparse(n,n);
for k = 1:nt
    r = randperm(n,floor(pm*n))';               % draw random sample of nodes in p{k}
    [pr,tr,pbr] = delaunay_C2(p(r,:,k),20);
    I = repmat(eye(2),[1 1 size(tr,1)]);        % id tensor
    [Kr,Mr] = assemble2(pr,tr,pbr,I);           % assembly of stiffness and mass
    [Ir,Jr,Sr] = find(Kr);  K = K + sparse(r(Ir),r(Jr),Sr,n,n);
    [Ir,Jr,Sr] = find(Mr);  M = M + sparse(r(Ir),r(Jr),Sr,n,n);
end; 
% remove all zero rows and columns
S = sum(abs(K)); I = find(abs(S)>eps); K = K(I,I); M = M(I,I); 

%% eigenproblem
[V,L] = eigs(K,M,15,'SM'); 
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord); 

figure(1); plot(lam,'*'); 
[pI,tI,pbI] = trimesh(p(I,:,1));
figure(2); plotf(pI,tI,pbI,normed(V(:,2)),0); colorbar

%% coherent partition
nc = 8;                                         % number of clusters
W = kmeans(V(pbI(:,2),1:nc),nc,'Replicates',20); % kmeans clustering
figure(3); clf; scatter(p(:,1,1),p(:,2,1),20,W,'filled');
axis equal; axis tight; colormap(jet(nc));

%% advected partition
figure(4); clf; scatter(p(:,1,nt),p(:,2,nt),20,W,'filled'); 
axis equal; axis(dom(:)); colormap(jet(nc));

