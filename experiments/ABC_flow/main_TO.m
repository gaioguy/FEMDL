addpath('../../src'); clear all; clc

%% ABC inverse flow map
t0 = 0; tf = 1; tspan = [tf,t0]; 
T1 = @(x) flow_map_3d(@ABC,x,tspan);

%% nodes
n = 25; x = (0:n-1)/(n-1);
[p1,p2,p3] = meshgrid(x,x,x);
p = [p1(:),p2(:),p3(:)]; 
% periodic boundary conditions
np = reshape(1:(n-1)^3,n-1,n-1,n-1);
np(n,:,:) = np(1,:,:); 
np(:,n,:) = np(:,1,:); 
np(:,:,n) = np(:,:,1);
pb = [(1:n*n*n)' np(:)];        % nodes pb(:,1) are identified with pb(:,2)

%% triangulation
tri = delaunayTriangulation(p); 
t = tri.ConnectivityList;

%% integration and assembly
tic; P = T1(p); T1p = P(:,[2,4,6]); toc
A = kron([1 0 0 1 0 1],ones(size(t,1),1));  % 3 x 3 idenitiy matrix
tic; [D,M] = assemble_3d(p,t,pb,A); toc

%% compute transfer operator and dynamic Laplacian
Phi = get_Phi_3d(p,t,pb);
tic; Alpha = get_alpha_3d(tri,pb,T1p,Phi); toc
DL = 0.5*(D + Alpha'*D*Alpha);

%% eigenproblem
tic; [V,L] = eigs(DL,M,3,'SM'); toc
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord);

%% plot spectrum
figure(1); clf; plot(lam,'*'); 
xlabel('$k$'); ylabel('$\lambda_k$')

%% plot eigenvector
figure(2); clf; clear w; v = V(:,3);  theta = 0.18;
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


