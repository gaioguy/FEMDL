addpath('../../src/2d'); clear all; init

%% flow map etc.
t0 = 0; tf = 0.6; 
vf  = @double_gyre; v = @(t,x) rm(vf(t,mr(x)));          % vector field
T   = @(a,x) at_tf(flowmap(vf,x,[t0 a],1e-7));           % flow map
DT  = @(a,x) D(@(x) T(a,x),x);                           % space derivative of flow map
Tp  = @(a,x) v(a,T(a,x));                                % parameter derivative of T

%% triangulation
nx = 100; ny = nx; 
p0 = grid2(nx,ny);
mesh0 = trimesh(p0); mesh0.b = [];

%% compute u0 
p1 = T(tf,p0);
mesh1 = trimesh(p1); mesh1.b = [];
[V0,lam0,K,M] = solve_TO(mesh0,mesh1);
lam0(2)
u0 = V0(:,2); 
u0 = u0/sqrt(u0'*M*u0);                   % L2 normalization
figure(1); clf; plotf2(mesh0,u0);

%% compute linear response
W = Tp(tf,p0);
L = assemble_lr(mesh1,W); 
L = 0.5*(L+L');
tmp = [(K-lam0(2)*M) -M*u0; (M*u0)' 0]\[-L*u0; 0]; 
udot0 = tmp(1:end-1);  
lamdot0 = tmp(end) 
figure(2); clf; plotf2(mesh0,udot0); 

%% prediction of u_eps by linear response
dt = 0.2; 
u0lr = u0 + dt*udot0;
figure(3); clf; plotf2(mesh0,u0lr); colorbar; 

%% exact u_eps, errors of FD approximation
p1 = T(tf+dt,p0);
mesh1 = trimesh(p1); mesh1.b = [];
[Ve,lame] = solve_TO(mesh0,mesh1);
ue = Ve(:,2);
figure(4); clf; plotf2(mesh0,ue); colorbar; 
L2_error = sqrt((u0lr-ue)'*M*(u0lr-ue))
eigerror = abs(lame(2)-(lam0(2)+dt*lamdot0))/abs(lam0(2))


