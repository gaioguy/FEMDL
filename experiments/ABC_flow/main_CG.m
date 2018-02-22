addpath('../../src'); clear all; clc

%% ABC flow 
t0 = 0; tf = 1; tspan = [t0,tf];
% inverse Cauchy-Green tensor
invCG = @(x) inv_CG_3d(@ABC,x,tspan); 

%% nodes
n = 25; x = (0:n-1)/(n-1);
[p1,p2,p3] = meshgrid(x,x,x);
p = [p1(:),p2(:),p3(:)]; 
% periodic boundary conditions
np = reshape(1:(n-1)^3,n-1,n-1,n-1);
np(n,:,:) = np(1,:,:); 
np(:,n,:) = np(:,1,:); 
np(:,:,n) = np(:,:,1);
pb = [(1:n^3)' np(:)];                   

%% triangulation
tri = delaunayTriangulation(p); 
t = tri.ConnectivityList;

%% assembly
deg = 1;  % quad degree
tic; G = inv_CG_quad_3d(p,t,invCG,deg); toc
tic; [D,M] = assemble_3d(p,t,pb,G); toc

%% eigenproblem
tic; [V,L] = eigs(D,M,3,'SM'); toc
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord); lam

%% plot spectrum
figure(1); clf; plot(lam,'*'); 
xlabel('k'); ylabel('$\lambda_k$')

%% plot eigenvector
figure(2); clf; clear w; v = V(:,3);  theta = 0.2;
w(pb(:,1)) = v(pb(:,2)); w = w/norm(w,inf); ww = reshape(w, n, n, n);
h1 = patch(isosurface(p1,p2,p3,ww,theta),'FaceAlpha',0.5);
isonormals(p1,p2,p3,ww,h1); 
h1.FaceColor = 'red'; h1.EdgeColor = 'none';
patch(isocaps(p1,p2,p3,ww,theta,'above'),'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5)

h2 = patch(isosurface(p1,p2,p3,ww,-theta),'FaceAlpha',0.5);
isonormals(p1,p2,p3,ww,h2);
h2.FaceColor = 'blue'; h2.EdgeColor = 'none'; 
patch(isocaps(p1,p2,p3,ww,-theta,'below'),'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5)
patch(isocaps(p1,p2,p3,ww,theta,'below'),'FaceColor','interp','EdgeColor','none','FaceAlpha',0.25)

axis equal; axis([0 1 0 1 0 1]); view(-145,25); caxis([-1,1])
camlight right; camlight left; lighting gouraud; colorbar
xlabel('$x$'); ylabel('$y$'); zlabel('$z$');
set(gcf, 'PaperPosition', [0 0 16 10],'PaperSize', [16 10]);

