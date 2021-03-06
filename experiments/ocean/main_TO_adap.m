addpath('../../src/2d'); init

%% ocean vector field 
% data by SSALTO/DUACS, distributed by AVISO, http://www.aviso.oceanobs.com/duacs
load('Ocean_geostrophic_velocity.mat','lon','lat','UT','VT','time');
U = griddedInterpolant({lon,lat,time},permute(UT,[2,1,3]),'cubic','none');
V = griddedInterpolant({lon,lat,time},permute(VT,[2,1,3]),'cubic','none');
v = @(t,x) ocean(t,x,U,V);

%% flow 
t0 = time(1); tf = t0+90; 
nt = 11; tspan = linspace(t0,tf,nt);
T = @(x) flowmap(v,x,tspan);                    % flow map

%% (initial) triangulation
dom = [-4 -34; 6 -28]; dx = diff(dom);          % domain 
nx = 200; ny = dx(2)/dx(1)*nx; n = nx*ny;       % number of grid points
g = grid2(nx,ny)*diag(dx) + dom(1,:);           % grid 
m = trimesh(g);

%% assembly (for missing data case)
g = T(g);                                       % advected grid
pm = 1;                                       % percentage of nodes to remove
K = sparse(n,n); M = sparse(n,n); 
for k = 1:nt
    r = randperm(n,floor(pm*n))';               % pm*n random integers in 1:n
    mr = trimesh(g(r,:,k));                     % triangulation on corr. nodes
    I = repmat(eye(2),[1 1 size(mr.t,1)]);      % id tensor
    mt.p = g(:,:,k); mt.t = r(mr.t); mt.pb = [1:n; 1:n]'; 
    [Kt,Mt] = assemble2(mt,I);                  % assembly of stiffness and mass
    K = K + Kt/nt; M = M + Mt/nt; 
end

%% Dirichlet boundary condition
b = m.b;
K(b,:) = 0; K(:,b) = 0; M(b,:) = 0; M(:,b) = 0; 
K(b,b) = speye(length(b),length(b));

%% remove all zero rows and columns (sparse data case)
S = sum(abs(K)); I = find(abs(S)>eps); K = K(I,I); M = M(I,I); 

%%  eigenproblem
[V,L] = eigs(K,M,10,'SM');
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord);

figure(1); plot(lam,'*');                       % plot spectrum
mI = trimesh(g(I,:,1));
figure(2), plotf(mI,V(:,5),0);           % plot eigenvector

%% coherent partition
nc = 4;
W = kmeans(V(:,1:nc), nc, 'Replicates', 10);          % kmeans clustering
figure(3), clf, scatter(g(I,1,1),g(I,2,1),10,W,'filled');
axis equal, axis tight, colormap(jet(nc));

%% advected partition
figure(4); clf; scatter(g(I,1,end),g(I,2,end),10,W,'filled'); 
axis equal; axis([-12 4 -34 -25]); colormap(jet(nc));
