addpath('../../src/2d'); clear all; init

%% standard map
a = 0.98;
T = @(a,x) mod([x(:,1) + x(:,2) + a*sin(x(:,1)), ...
                x(:,2) + a*sin(x(:,1))],2*pi);
DT = @(a,x) rowapply(@(x) [1 + a*cos(x(1)) 1;   
                               a*cos(x(1)) 1], x);
Tp = @(a,x) [sin(x(:,1)), sin(x(:,1))];

%% regular triangulation
nx = 100; ny = nx; 
dom = [0 0; 2*pi 2*pi]; dx = diff(dom);
p0 = grid2(nx,ny)*diag(dx*(nx-2)/(nx-1)) + dom(1,:);
mesh0 = delaunay_T2(p0, dx(1), dx(2));
deg = 2;

%% compute u0 
p1 = T(a,p0);
mesh1 = delaunay_T2(p1,2*pi,2*pi);
[V0,lam0,K,M] = solve_TO(mesh0,mesh1);
u0 = V0(:,2); 
u0 = u0/sqrt(u0'*M*u0);         % L2 normalization
figure(1); clf; plotf2(mesh0,u0); 

%% compute linear response
W = Tp(a,p0);
L = assemble_lr(mesh1,W); 
L = 0.5*(L+L');
lamdot0 = (u0'*L*u0)/(u0'*M*u0)
tmp = [(K-lam0(2)*M) -M*u0; (M*u0)' 0]\[-L*u0; 0]; 
udot0 = tmp(1:end-1); 
figure(2); clf; plotf2(mesh0,udot0); 

%% comparison of different methods to obtain coherent set prediction
% exact second eigenvector of perturbed map
da = 0.5; 
p1 = T(a+da,p0);
mesh1 = delaunay_T2(p1,2*pi,2*pi);
[Ve,lame] = solve_TO(mesh0,mesh1);
ue = Ve(:,2);           
figure(3); clf; plotf2(mesh0,ue);
% approximation by linear response
ulr = u0 + da*udot0;     
figure(4); clf; plotf2(mesh0,ulr);  

