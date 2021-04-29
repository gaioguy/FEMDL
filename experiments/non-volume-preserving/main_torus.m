addpath('../../src','../../src/2d'); clear all; init

%% periodic triangulation
nx = 100;  ny = nx;  n = (nx-1)*(ny-1); dx = 2*pi/(nx-1); dy = dx;
% [xi,yi] = meshgrid(linspace(0,2*pi-dx,nx-1),linspace(0,2*pi-dy,ny-1));
% p0{1} = [xi(:) yi(:)];  
p0(:,:,1) = rand((nx-1)^2,2)*diag([2*pi-dx,2*pi-dy]);
m{1} = delaunay_T2(p0(:,:,1), 2*pi, 2*pi);

%% flow map and initial density
t0 = 0; tf = 1; nt = 2; 
T3 = @(x) [x(:,1) + 0.3*cos(2*x(:,1)), x(:,2)]; 
T4 = @(x) [x(:,1) + x(:,2), x(:,2) + 8*sin(x(:,1) + x(:,2))];
T = @(x) mod(T4(mod(T3(x),2*pi)),2*pi);
h = @(x) 1/(8*pi^2)*(sin(x(:,2)-pi/2)+2);

%% time integration and mesh construction
tic; for k = 1:nt-1
    p0(:,:,k+1) = T(p0(:,:,k)); 
    m{k+1} = delaunay_T2(p0(:,:,k+1), 2*pi, 2*pi);
end; toc

%% assembly
tic; D = sparse(n,n); M = sparse(n,n); 
a = quad_basis(m{1}.p,m{1}.t,m{1}.pb,h);
for k = 1:nt
    A = kron([1 0 0 1],ones(size(m{k}.t,1),1));      % 2 x 2 identity matrix
    [Dt{k},Mt{k},lt{k}] = assemble_weighted2(m{k},A,a);
    D = D + Dt{k}; M = M + Mt{k}; 
end, toc

%% eigenproblem
tic; [V,L] = eigs(D,M,10,'sm'); toc
[lam,order] = sort(diag(L),'descend'); V = V(:,order); 

%% plot spectrum
figure(1); clf; plot(lam,'*'); axis tight, axis square
xlabel('$k$'); ylabel('$\lambda_k$')

%% plot eigenvector ...
figure(1), clf;  clear w; cl = [-0.8, -0.4, 0, 0.4, 0.8];
w(m{1}.pb(:,1)) = V(m{1}.pb(:,2),2); w = w/norm(w,inf); 
tricontour(m{1}.p,m{1}.t,normed(w'),cl); axis equal; axis([0 2*pi 0 2*pi]);  box on

%% ... and its image
figure(2), clf;  clear w; cl = [-0.8, -0.4, 0, 0.4, 0.8];
w(m{2}.pb(:,1)) = V(m{2}.pb(:,2),2); 
tricontour(m{2}.p,m{2}.t,normed(w'),cl); axis equal; axis([0 2*pi 0 2*pi]);  box on

%% plot initial density
figure(3); clf; colormap summer; cmap = colormap; cmap = cmap(end:-1:1,:); colormap(cmap);
plotf(m{1},a./lt{1},0); axis equal; axis([0 2*pi 0 2*pi]); box on; 
caxis([0.005 0.09]); colorbar

%% ... and its image
figure(4); clf; colormap(cmap);
plotf(m{2},a./lt{2},0); axis equal; axis([0 2*pi 0 2*pi]); box on; colorbar
caxis([0.005 0.09]); colorbar

