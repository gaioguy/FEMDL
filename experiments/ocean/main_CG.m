addpath('../../src'); clear all

%% ocean vector field
% The altimeter data are produced by SSALTO/DUACS and distributed by 
% AVISO, (http://www.aviso.oceanobs.com/duacs).
load('Ocean_geostrophic_velocity.mat','lon','lat','UT','VT','time');
U = griddedInterpolant({lon,lat,time},permute(UT,[2,1,3]),'cubic','none');
V = griddedInterpolant({lon,lat,time},permute(VT,[2,1,3]),'cubic','none');
vf = @(t,x) ocean_vf(t,x,U,V);
t0 = time(1); tf = t0 + 90; nt = 2;
CG = @(x) cg_tensor(vf,x,linspace(t0,tf,nt));

%% triangulation
xmin = -4; xmax = 6; ymin = -34; ymax = -28; 
nx = 250; ny = 0.6*nx; n = nx*ny; 
x1 = linspace(xmin,xmax,nx); y1 = linspace(ymin,ymax,ny);
[X,Y] = meshgrid(x1,y1); p = [X(:) Y(:)];           % nodes
pb = [1:n; 1:n]';                                   % no periodic boundary
tri = delaunayTriangulation(p); t = tri.ConnectivityList; 
b = unique(freeBoundary(tri));                      % b = boundary nodes

%% assembly 
deg = 1;                                            % degree of quadrature
G = compute_G(p,t,CG,deg); 
[D,M] = assemble(p,t,pb,G);

%% Dirichlet boundary condition
D(b,:) = 0; D(:,b) = 0; M(b,:) = 0; M(:,b) = 0;
D(b,b) = speye(length(b),length(b));

%% solve eigenproblem
[V,L] = eigs(D,M,20,'SM');
[lam,order] = sort(diag(L),'descend'); V = V(:,order);

%% plot spectrum
figure(1); clf; plot(lam,'s','markerfacecolor','b');
xlabel('$k$'); ylabel('$\lambda_k$')

%% plot eigenvectors
figure(2), clf; plotev(t,p,pb,V(:,7),0); colorbar
xlabel('lon [$^\circ$]'); ylabel('lat [$^\circ$]'); 

%% compute partition
nx1 = 200; ny1 = 0.6*nx1; x1 = linspace(xmin,xmax,nx1); y1 = linspace(ymin,ymax,ny1);
[X1,Y1] = meshgrid(x1,y1); 
V1 = eval_p1(p,V(:,1:4),[X1(:) Y1(:)]);         % evaluate eigenvectors on grid
idx = kmeans(V1, size(V1,2));                   % kmeans clustering

%% plot partition
figure(3); clf; surf(X1,Y1,reshape(idx,ny1,nx1)); view(2); shading flat
axis equal; axis tight; xlabel('$x$'); ylabel('$y$'); colorbar



