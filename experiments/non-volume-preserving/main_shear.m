addpath('../../src'); clear all

%% flow map T1
t0 = 0; tf = 1; nt = 2;  
T = @(x) [mod(x(:,1) + (cosh(2*x(:,2))-1)/2,4), x(:,2)];
h = @(x) 1/8*(sin(pi*x(:,1))+2);

%% flow map T2
t0 = 0; tf = 1; nt = 2;  
T = @(x) [mod(x(:,1) + x(:,2),4), x(:,2) + 0.1*sin(2*pi*x(:,2))];
h = @(x) ones(size(x,1),1)/4;

%% periodic triangulation
nx = 60;  ny = nx/4;  n = (nx-1)*ny; dx = 4/(nx-1); 
[xi,yi] = meshgrid(linspace(0,4-dx,nx-1),linspace(0,1,ny));
p0{1} = [xi(:) yi(:)];  
[p{1},t{1},pb{1}] = delaunay_C2(p0{1},[0,4]); 

%% time integration and mesh construction
tic; for k = 1:nt-1, 
    p0{k+1} = T(p0{k}); 
    [p{k+1},t{k+1},pb{k+1}] = delaunay_C2(p0{k+1},[0,4]);
end; toc

%% assembly
tic; D = sparse(n,n); M = sparse(n,n); 
a = quad_basis(p{1},t{1},pb{1},h);
for k = 1:nt
    CG = kron([1 0 0 1],ones(size(t{k},1),1));      % 2 x 2 identity matrix
    [Dt,Mt,lt{k}] = assemble_weighted(p{k},t{k},pb{k},CG,a);
    D = D + Dt; M = M + Mt; 
end, toc

%% eigenproblem
tic; [V,L] = eigs(D,M,5,'SM'); toc
[lam,order] = sort(diag(L),'descend'); V = V(:,order); 

%% plot eigenvector
figure(1), clf;  clear w
w(pb{1}(:,1)) = V(pb{1}(:,2),2); w = w/norm(w,inf); 
subplot(211); tricontour(t{1},p{1}(:,1),p{1}(:,2),-w,7);
axis equal, axis([0 4 0 1]), box on
w(pb{2}(:,1)) = V(pb{2}(:,2),2); w = w/norm(w,inf); 
subplot(212); tricontour(t{2},p{2}(:,1),p{2}(:,2),w,7); 
axis equal, axis([0 4 0 1]), box on

%% image density
figure(2); clf; colormap summer; cmap = colormap; cmap = cmap(end:-1:1,:); colormap(cmap);
subplot(211); plotf(p{1},t{1},pb{1},a./lt{1},0); box on
subplot(212); plotf(p{2},t{2},pb{2},a./lt{2},0); box on



