addpath('../../src/2d'); clear all; init

%% standard map
a = 0.971635;
%a = 8;
T = @(x) mod([x(:,1) + x(:,2) + a*sin(x(:,1)), x(:,2) + a*sin(x(:,1))],2*pi);
DT = @(x) rowapply(@(x) [1 + a*cos(x(1)) 1; a*cos(x(1)) 1], x);

dL = @(DT) 0.5*(eye(2) + inv(DT)*inv(DT)');
dLx = @(x) fapply1(dL, DT(x));

%% regular triangulation
nx = 50; ny = nx; 
dom = [0 0; 2*pi 2*pi]; dx = diff(dom);
p0 = grid2(nx,ny)*diag(dx*(nx-2)/(nx-1)) + dom(1,:);
mesh = delaunay_T2(p0, dx(1), dx(2));
deg = 2;

%% assembly
A = triquad(mesh,@(x) dLx(x),deg);      % integrate inverse Cauchy Green tensor
[K,M] = assemble2(mesh,A);              % assemble matrices

%% solve eigenproblem               
[V,lam] = eigs(K,M,6,'SM'); 
[lam,ord] = sort(diag(lam),'descend'); V = V(:,ord);

figure(1); clf; plot(lam,'*'); 
figure(2); clf; plotf(mesh,normed(V(:,2)),0); colorbar; 


