addpath('../../src/2d'); init

%% flow map etc.
t0 = 0; tf = 0.6; 
vf  = @double_gyre; v = @(t,x) rm(vf(t,mr(x)));          % vector field

T   = @(a,x) at_tf(flowmap(vf,x,[t0 a],1e-5));           % flow map
DT  = @(a,x) D(@(x) T(a,x),x);                           % space derivative of flow map
Tp  = @(a,x) v(a,T(a,x));                                % parameter derivative of T
DTp = @(a,x) D(@(x) Tp(a,x),x);                          % space derivative of T'

dL   = @(DT) 0.5*(eye(2) + inv(DT)*inv(DT)');            % dynamic Laplacian
dLx  = @(a,x) fapply1(dL, DT(a,x));                      % evaluated at each row of x
dLp  = @(DT,DTp) -sym(inv(DT)*DTp*inv(DT)*inv(DT)');     % parameter derivative of dL
dLpx = @(a,x) fapply2(dLp, DT(a,x), DTp(a,x));           % evaluated at each row of x

%% triangulation
nx = 100; ny = nx; 
p = grid2(nx,ny);
m0 = trimesh(p); m0.b = [];

%% compute second eigenvector u0
deg = 5;
[V0,lam0,K,M] = solve_CG(m0,@(x) dLx(tf,x),deg); 
lam0(2)
u0 = V0(:,2); 
figure(1); clf; plotf2(m0,u0);  colorbar

%% compute udot via finite differencing
e = 1e-4; 
[Ve,lame] = solve_CG(m0,@(x) dLx(tf+e,x),deg); 
lamdot_fd = (lame(2)-lam0(2))/e
ue = Ve(:,2);
ue = ue/(ue'*M*u0);             % determine scaling of u_e such that in V_0
udot0_fd = (ue-u0)/e;
figure(2); plotf2(m0,udot0_fd); colorbar; 

%% compute udot0 by CG approach
Aq = triquad(m0,@(x) dLpx(tf,x),deg); 
L = assemble2(m0,Aq);        
tmp = [(K-lam0(2)*M) -M*u0; (M*u0)' 0]\[-L*u0; 0]; 
udot0 = tmp(1:end-1); 
lamdot0 = tmp(end) 
figure(4); clf; plotf2(m0,udot0); colorbar;

%% error of FD approximation 
Id = @(x) eye(2);
[~,~,KId] = solve_CG(m0, Id, deg); 
L2_error = sqrt((udot0_fd-udot0)'*M*(udot0_fd-udot0))
H1_error = sqrt((udot0_fd-udot0)'*(M-KId)*(udot0_fd-udot0))
eig_error = abs(lamdot_fd-lamdot0)

%% prediction of u_eps by linear response
dt = 0.2; 
u0lr = u0 + dt*udot0;
u0lr = u0lr/sqrt(u0lr'*M*u0lr);
figure(4); clf; plotf2(m0,u0lr); colorbar; 

%% exact u_eps
[Ve, lame] = solve_CG(m0, @(x) dLx(tf+dt,x), deg); 
lame(2)
ue = Ve(:,2);
ue = ue/sqrt(ue'*M*ue);
figure(5); clf; plotf2(m0,ue); colorbar;
L2_error = sqrt((u0lr-ue)'*M*(u0lr-ue))
eigerror = abs(lame(2)-(lam0(2)+dt*lamdot0))/abs(lam0(2))

%% determine level set of optimal dynamic isoperimetric ratio
T1 = @(x) T(tf+dt,x);
m1 = trimesh(T1(p));
[c,minipr,minc,iprs] = curve_of_optimal_ipr(m0,T1,u0,0,max(u0));

%% velocity field for level set evolution 
x = unique(p(:,1)); y = unique(p(:,2));
[X,Y] = meshgrid(x,y);
U0 = reshape(u0,nx,ny);
Ue = reshape(ue,nx,ny);
U0P = reshape(udot0,nx,ny);
[gx,gy] = gradient(U0,1/nx,1/ny);
I = 1:6:nx;
normg = sqrt(gx(I,I).^2+gy(I,I).^2);
figure(6); clf
quiver(X(I,I),Y(I,I),-gx(I,I).*U0P(I,I)./normg,-gy(I,I).*U0P(I,I)./normg);
axis equal; axis tight; xlabel('$x$'); ylabel('$y$');
c = minc;
hold on; contour(X,Y,U0,c*[1 1],'k','linewidth',2);
contour(X,Y,U0+dt*U0P,c*[1 1],'g','linewidth',2);
contour(X,Y,Ue,c*[1 1],'r--','linewidth',2);

%% prediction of u_eps by Taylor approximation of the flow
dtvf = @(t,x) df_fd(@(t) vf(t,x), t, 1)';
vft  = @(t,x) (t<tf).*vf(t,x) + (t>=tf).*(vf(tf,x) + (t-tf).*dtvf(tf,x));
Tt   = @(a,x) at_tf(flowmap(@(t,x) vft(t,x), x, [t0 a]));
DTt  = @(a,x) D(@(x) Tt(a,x),x);                           
dLtx  = @(a,x) fapply1(dL, DTt(a,x));

V0t = solve_CG(m0, @(x) dLtx(tf+dt,x), deg);
u0t = V0t(:,2);
figure(7); clf; plotf2(m0,u0t); colorbar; 

