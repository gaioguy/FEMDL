clear all

%% ABC flow map
tspan = [0,1]; nt = 2;
T = @(x) flow_map_3d(@ABC,x,tspan);

%% data points
n = 25; x = (0:n-2)/(n-1);
[px,py,pz] = meshgrid(x,x,x);
p0 = [px(:),py(:),pz(:)]; 

%% time integration
tic; P = T(p0); toc
p1 = [P(:,2) P(:,4) P(:,6)]; 

%% triangulations
[p0,t0,pb0] = delaunay_T3(p0, [0 1 0 1 0 1]);
[p1,t1,pb1] = delaunay_T3(p1, [0 1 0 1 0 1]);

%% assembly
A = kron([1 0 0 1 0 1],ones(size(t0,1),1));     % 3 x 3 identity matrix
tic; [D0,M] = assemble_3d(p0,t0,pb0,A); toc
A = kron([1 0 0 1 0 1],ones(size(t1,1),1));     % 3 x 3 identity matrix
tic; [D1,M1] = assemble_3d(p1,t1,pb1,A); toc
D = (D0 + D1)/2; 

%% eigenproblem
tic; [V,L] = eigs(D,M,4,'SM'); toc
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord);

%% plot spectrum
figure(1); clf; plot(lam,'*'); 
xlabel('$k$'); ylabel('$\lambda_k$')

%% plot eigenvector
figure(2); clf; clear w; v = -V(:,3);  theta = 0.4;
w(pb0(:,1)) = v(pb0(:,2)); 
v = v/norm(v,inf); ww = reshape(v, n-1, n-1, n-1);
h1 = patch(isosurface(px,py,pz,ww,theta),'FaceAlpha',0.5);
isonormals(px,py,pz,ww,h1); 
h1.FaceColor = 'red'; h1.EdgeColor = 'none';
patch(isocaps(px,py,pz,ww,theta,'above'),'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5)

h2 = patch(isosurface(px,py,pz,ww,-theta),'FaceAlpha',0.5);
isonormals(px,py,pz,ww,h2);
h2.FaceColor = 'blue'; h2.EdgeColor = 'none'; 
patch(isocaps(px,py,pz,ww,-theta,'below'),'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5)
patch(isocaps(px,py,pz,ww,theta,'below'),'FaceColor','interp','EdgeColor','none','FaceAlpha',0.25)

axis equal; axis([0 1 0 1 0 1]); view(-145,25); caxis([-1,1])
camlight right; camlight left; lighting gouraud; colorbar
xlabel('$x$'); ylabel('$y$'); zlabel('$z$');


