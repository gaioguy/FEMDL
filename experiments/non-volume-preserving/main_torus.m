addpath('../../src'); clear all

%% periodic triangulation
nx = 40;  ny = nx;  n = (nx-1)*(ny-1); dx = 2*pi/(nx-1); dy = dx;
[xi,yi] = meshgrid(linspace(0,2*pi-dx,nx-1),linspace(0,2*pi-dy,ny-1));
p0{1} = [xi(:) yi(:)];  
tic; [p{1},t{1},pb{1}] = delaunay_T2(p0{1},[0,2*pi,0,2*pi]); toc

%% flow map and initial density
t0 = 0; tf = 1; nt = 2; 
T3 = @(x) [x(:,1) + 0.3*cos(2*x(:,1)), x(:,2)]; 
T4 = @(x) [x(:,1) + x(:,2), x(:,2) + 8*sin(x(:,1) + x(:,2))];
T = @(x) mod(T4(mod(T3(x),2*pi)),2*pi);
h = @(x) 1/(8*pi^2)*(sin(x(:,2)-pi/2)+2);

%% time integration and mesh construction
tic; for k = 1:nt-1
    p0{k+1} = T(p0{k}); 
    [p{k+1},t{k+1},pb{k+1}] = delaunay_T2(p0{k+1},[0,2*pi,0,2*pi]);
end; toc

%% assembly
tic; D = sparse(n,n); M = sparse(n,n); 
a = quad_basis(p{1},t{1},pb{1},h);
for k = 1:nt
    A = kron([1 0 0 1],ones(size(t{k},1),1));      % 2 x 2 identity matrix
    [Dt{k},Mt{k},lt{k}] = assemble_weighted(p{k},t{k},pb{k},A,a);
    D = D + Dt{k}; M = M + Mt{k}; 
end, toc

%% eigenproblem
tic; [V,L] = eigs(D,M,10,'sm'); toc
[lam,order] = sort(diag(L),'descend'); lam, V = V(:,order); 

%% plot spectrum
figure(1); clf; plot(lam,'*'); axis tight, axis square
xlabel('$k$'); ylabel('$\lambda_k$')

%% plot eigenvector and its image
figure(1), clf;  clear w; cl = [-0.8, -0.4, 0, 0.4, 0.8];
w(pb{1}(:,1)) = V(pb{1}(:,2),2); w = w/norm(w,inf); 
tricontour(t{1},p{1}(:,1),p{1}(:,2),-w,cl); axis equal;axis([0 2*pi 0 2*pi]);  box on
%% and its image
figure(2), clf;  clear w; cl = [-0.8, -0.4, 0, 0.4, 0.8];
w(pb{2}(:,1)) = V(pb{2}(:,2),2); w = w/norm(w,inf); 
tricontour(t{2},p{2}(:,1),p{2}(:,2),-w,cl);axis equal; axis([0 2*pi 0 2*pi]);  box on

%% plot initial density
figure(3); clf; colormap summer; cmap = colormap; cmap = cmap(end:-1:1,:); colormap(cmap);
plotf(p{1},t{1},pb{1},a./lt{1},0); axis equal; axis([0 2*pi 0 2*pi]); box on; 
caxis([0.005 0.09]); colorbar

%% ... and its image
figure(4); clf; colormap summer; cmap = colormap; cmap = cmap(end:-1:1,:); colormap(cmap);
plotf(p{2},t{2},pb{2},a./lt{2},0); axis equal; axis([0 2*pi 0 2*pi]); box on; colorbar
caxis([0.005 0.09]); colorbar

%% compute partition
nx1 = 200; ny1 = nx1; 
x1 = linspace(0,2*pi,nx1); y1 = linspace(0,2*pi,ny1);
[X1,Y1] = meshgrid(x1,y1); 
V1 = eval_p1(p0{1},V(:,1:3),[X1(:) Y1(:)]);       % evaluate eigenvectors on grid
idx = kmeans(V1, 3,'Replicates',20);       % kmeans clustering

%% plot partition
figure(3); clf; surf(X1,Y1,reshape(idx,ny1,nx1)); view(2); shading flat
axis equal; axis tight; xlabel('$x$'); ylabel('$y$'); 


