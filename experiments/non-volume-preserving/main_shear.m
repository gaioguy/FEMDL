addpath('../../src/2d'); init

%% flow map T1
t0 = 0; tf = 1; nt = 2;  
T = @(x) [mod(x(:,1) + (cosh(2*x(:,2))-1)/2,4), x(:,2)];
h = @(x) 1/8*(sin(pi*x(:,1))+2);

% %% flow map T2
% t0 = 0; tf = 1; nt = 2;  
% T = @(x) [mod(x(:,1) + x(:,2),4), x(:,2) + 0.1*sin(2*pi*x(:,2))];
% h = @(x) ones(size(x,1),1)/4;

%% periodic triangulation
nx = 60;  ny = nx/4;  n = (nx-1)*ny; dx = 4/(nx-1); 
[xi,yi] = meshgrid(linspace(0,4-dx,nx-1),linspace(0,1,ny));
p0(:,:,1) = [xi(:) yi(:)];  
m{1} = delaunay_C2(p0(:,:,1),4); 
p(:,:,1) = m{1}.p; t{1} = m{1}.t; pb{1} = m{1}.pb;

%% time integration and mesh construction
tic; for k = 1:nt-1, 
    p0(:,:,k+1) = T(p0(:,:,k)); 
    m{k+1} = delaunay_C2(p0(:,:,k+1),4);
    p(:,:,k+1) = m{k+1}.p; t{k+1} = m{k+1}.t; pb{k+1} = m{k+1}.pb;
end; toc

%% assembly
tic; D = sparse(n,n); M = sparse(n,n); 
a = quad_basis(p(:,:,1),t{1},pb{1},h);
for k = 1:nt
    CG = kron([1 0 0 1],ones(size(t{k},1),1));      % 2 x 2 identity matrix
    [Dt{k},Mt{k},lt{k}] = assemble_weighted(p(:,:,k),t{k},pb{k},CG,a);
    D = D + Dt{k}; M = M + Mt{k}; 
end, toc

%% eigenproblem
tic; [V,L] = eigs(D,Mt{1},5,'SM'); toc
[lam,order] = sort(diag(L),'descend'); lam, V = V(:,order); 

%% plot eigenvector
figure(1), clf;  clear w
w(pb{1}(:,1)) = V(pb{1}(:,2),2); 
subplot(211); tricontour(p(:,:,1),t{1},-normed(w'),7);
axis equal, axis([0 4 0 1]), box on
w(pb{2}(:,1)) = V(pb{2}(:,2),2); 
subplot(212); tricontour(p(:,:,2),t{2},normed(w'),7); 
axis equal, axis([0 4 0 1]), box on

%% image density
figure(2); clf; colormap summer; cmap = colormap; cmap = cmap(end:-1:1,:); colormap(cmap);
subplot(211); plotf(p(:,:,1),t{1},pb{1},a./lt{1},0); box on
subplot(212); plotf(p(:,:,2),t{2},pb{2},a./lt{2},0); box on



