addpath('../../src/2d'); init
xmod20 = @(x) [mod(x(:,1),20) x(:,2)];

%% flow map 
t0 = 0; days = 60*60*24; tf = 40*days; 
v = @bickleyjet;
T = @(x) at_tf(flowmap(v, x, [t0 tf]));     % flow map
DL = @(DT) 0.5*(eye(2) + inv(DT)*inv(DT)'); % dynamic Laplace
DLx = @(x) fapply1(DL, D(T,x));             % evaluate DL at each row of x

%% triangulation
dom = [0 -3; 20 3]; dx = diff(dom);         % domain 
nx = 100; ny = dx(2)/dx(1)*nx;              % number of grid points
p = grid2(nx,ny)*diag(dx) + dom(1,:);       % grid 
[p,t,pb] = delaunay_C2(p, dx(1));

%% assembly
deg = 1;                                    % degree of quadrature
A = triquad(p,t,DLx,deg);                   % integrate DL on triangles
[K,M] = assemble2(p,t,pb,A);                % assemble stiffness and mass matrices

%% eigenproblem
[V,L] = eigs(K,M,15,'SM');
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord);
figure(1); clf; plot(lam,'*')
figure(2), plotf(p,t,pb,normed(V(:,3)),0); colorbar

%% coherent partition
nc = 8;                                         % number of clusters
W = kmeans(V(pb(:,2),1:nc),nc,'Replicates',20); % kmeans clustering
figure(3); clf; scatter(p(:,1),p(:,2),10,W,'filled');
axis equal; axis tight; colormap(jet(nc));

%% advected partition
Tp = xmod20(T(p));                                      % advect grid
figure(4); clf; scatter(Tp(:,1),Tp(:,2),10,W,'filled'); 
axis equal; axis(dom(:)); colormap(jet(nc));


