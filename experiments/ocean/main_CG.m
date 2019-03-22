addpath('../../src/2d'); init
h = 1e-4;
df = @(f,x,d) (f(x+h*d)-f(x-h*d)).'/(2*h);      % derivative of f at x in direction d
D = @(f,x) permute(cat(3,df(f,x,[1,0]),df(f,x,[0,1])),[1 3 2]); 

%% ocean vector field
% data by SSALTO/DUACS, distributed by AVISO, http://www.aviso.oceanobs.com/duacs
load('Ocean_geostrophic_velocity.mat','lon','lat','UT','VT','time');
UI = griddedInterpolant({lon,lat,time},permute(UT,[2,1,3]),'cubic','none');
VI = griddedInterpolant({lon,lat,time},permute(VT,[2,1,3]),'cubic','none');
v = @(t,x) ocean(t,x,UI,VI);

%% flow 
t0 = time(1); tf = t0 + 90; 
T = @(x) at_tf(flowmap(v,x,[t0 tf],1e-6));      % flow map at final time
DT = @(x) D(T,x);                               % space derivative of flow map

DL = @(DT) 0.5*(eye(2) + inv(DT)*inv(DT)');     % dynamic Laplace
DLx = @(x) fapply1(DL, D(T,x));                 % evaluate DL at each row of x

%% triangulation
dom = [-4 -34; 6 -28]; dx = diff(dom);          % domain 
n = 200; m = dx(2)/dx(1)*n;                     % number of grid points
g = grid2(n,m)*diag(dx) + dom(1,:);             % grid 
[p,t,pb,b] = trimesh(g);                        % triangular mesh

%% assembly 
deg = 1;                                        % degree of quadrature
A = triquad(p,t,DLx,deg);                       % integrate DL on triangles
[K,M] = assemble2(p,t,pb,A);                    % assemble stiffness and mass matrices
K(b,:) = 0; K(:,b) = 0; M(b,:) = 0; M(:,b) = 0; % Dirichlet bc
K(b,b) = speye(length(b),length(b));           

%% eigenproblem
[V,L] = eigs(K,M,10,'SM'); 
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord); 

figure(1); plot(lam,'*'); 
figure(2); plotf(p,t,pb,normed(V(:,3)),0); colorbar

%% coherent partition
nc = 4;                                         % number of clusters
W = kmeans(V(:,1:nc),nc,'Replicates',20);       % kmeans clustering
figure(3); scatter(g(:,1),g(:,2),10,W,'filled');
axis equal; axis tight; colormap(jet(nc));

%% advected partition
Tg = T(g);                                      % advect grid
figure(4); scatter(Tg(:,1),Tg(:,2),10,W,'filled'); 
axis equal; axis([-12 4 -34 -25]); colormap(jet(nc));


