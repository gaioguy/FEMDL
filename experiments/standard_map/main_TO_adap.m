addpath('../../src/2d'); clear all; init

%% standard map
a = 0.971635;
%a = 8;
T = @(x) mod([x(:,1) + x(:,2) + a*sin(x(:,1)), ...
                       x(:,2) + a*sin(x(:,1))], 2*pi);

%% regular triangulation
nx = 50; ny = nx; 
dom = [0 0; 2*pi 2*pi]; dx = diff(dom);
p0 = grid2(nx,ny)*diag(dx*(nx-2)/(nx-1)) + dom(1,:);
mesh0 = delaunay_T2(p0, dx(1), dx(2));

%% assembly
I = eye(2);
[K,M] = assemble2(mesh0,I);     
p1 = T(p0);
mesh1 = delaunay_T2(p1, dx(1), dx(2));
[K1,~] = assemble2(mesh1,I);
K = 0.5*(K+K1);

%% solve eigenproblem
[V,lam] = eigs(K,M,6,'SM'); 
[lam,ord] = sort(diag(lam),'descend'); V = V(:,ord);

figure(1); clf; plot(lam,'*'); 
figure(2); clf; plotf(mesh0,normed(V(:,2)),0); colorbar; 
