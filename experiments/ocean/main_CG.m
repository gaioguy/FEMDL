addpath('../../src'); clear all; clc; colormap jet
nrmd = @(v) v/norm(v,inf);

%% ocean vector field
% The altimeter data are produced by SSALTO/DUACS and distributed by 
% AVISO, (http://www.aviso.oceanobs.com/duacs).
load('Ocean_geostrophic_velocity.mat','lon','lat','UT','VT','time');
U = griddedInterpolant({lon,lat,time},permute(UT,[2,1,3]),'cubic','none');
V = griddedInterpolant({lon,lat,time},permute(VT,[2,1,3]),'cubic','none');

%% flow
vf = @(t,x) ocean(t,x,U,V);
t0 = time(1); tf = t0 + 90; 
nt = 2; tspan = linspace(t0,tf,nt);
CG = @(x) inv_CG(vf,x,tspan);

%% triangulation
xmin = -4; xmax = 6; ymin = -34; ymax = -28; 
nx = 100; ny = 0.6*nx; n = nx*ny;
x1 = linspace(xmin,xmax,nx); y1 = linspace(ymin,ymax,ny);
[X,Y] = meshgrid(x1,y1); p = [X(:) Y(:)];           % nodes
pb = [1:n; 1:n]';                                   % no periodic boundary

%% triangulation
tri = delaunayTriangulation(p); 
t = tri.ConnectivityList; 
b = unique(freeBoundary(tri));                      % b = boundary nodes

%% assembly 
deg = 5;                                            % degree of quadrature
tic; G = inv_CG_quad(p,t,CG,deg); toc
tic; [D,M] = assemble(p,t,pb,G); toc

%% Dirichlet boundary condition
D(b,:) = 0; D(:,b) = 0; M(b,:) = 0; M(:,b) = 0;
D(b,b) = speye(length(b),length(b));

%% eigenproblem
tic; [V,L] = eigs(D,M,10,'SM'); toc
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord); 

%% plot spectrum
figure(1); clf; plot(lam,'*'); axis square, axis tight
xlabel('$k$'); ylabel('$\lambda_k$')

%% plot eigenvector
figure(2), clf; plotf(p,t,pb,nrmd(V(:,4)),0); colorbar
xlabel('lon [$^\circ$]'); ylabel('lat [$^\circ$]'); 

%% compute partition
nx1 = 200; ny1 = 0.6*nx1; 
x1 = linspace(xmin,xmax,nx1); y1 = linspace(ymin,ymax,ny1);
[X1,Y1] = meshgrid(x1,y1); 
nc = 3;
tic; V1 = eval_p1(p,V(:,1:nc),[X1(:) Y1(:)]); toc    % evaluate eigenvectors on grid
tic; idx = kmeans(V1,nc+1,'Replicates',50); toc       % kmeans clustering

%% plot partition
figure(3); clf; surf(X1,Y1,reshape(idx,ny1,nx1)); view(2); shading flat
axis equal; axis tight; xlabel('lon [$^\circ$]'); ylabel('lat [$^\circ$]'); 
load cmap7; colormap(cmap); colorbar

%% advect abd plot LCS
figure(4); clf; hold on; colormap(cmap); caxis([1 nc+1])
T = @(x) flow_map(vf,x,tspan);
for l = 2:nc+1
    I = find(idx==l); 
    S = [X1(I) Y1(I)];
    TS = T(S); TS = TS(:,[2,4]);
    scatter3(TS(:,1),TS(:,2),zeros(length(I),1),10,cmap(l,:),'filled'); 
end
view(2); axis([-8 2 -33 -27])
xlabel('lon [$^\circ$]'); ylabel('lat [$^\circ$]');


