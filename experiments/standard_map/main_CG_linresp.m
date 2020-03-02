addpath('../../src/2d'); clear all; init

%% standard map
a = 0.98;
%a = 8;
T = @(a,x) mod([x(:,1) + x(:,2) + a*sin(x(:,1)), x(:,2) + a*sin(x(:,1))],2*pi);
DT = @(a,x) rowapply(@(x) [1 + a*cos(x(1)) 1; a*cos(x(1)) 1], x);
Tp = @(a,x) [sin(x(:,1)), sin(x(:,1))];
DTp = @(a,x) rowapply(@(x) [ cos(x(1)), 0; cos(x(1)), 0], x);

dL = @(DT) 0.5*(eye(2) + inv(DT)*inv(DT)');
dLx = @(a,x) fapply1(dL, DT(a,x));
dLp = @(DT,DTp) -sym(inv(DT)*DTp*inv(DT)*inv(DT)');
dLpx = @(a,x) fapply2(dLp, DT(a,x), DTp(a,x));

%% regular triangulation
nx = 100; ny = nx; 
dom = [0 0; 2*pi 2*pi]; dx = diff(dom);
p0 = grid2(nx,ny)*diag(dx*(nx-2)/(nx-1)) + dom(1,:);
mesh0 = delaunay_T2(p0, dx(1), dx(2));
deg = 2;

%% compute u0
[V0,lam0,K,M] = solve_CG(mesh0, @(x) dLx(a,x), deg); 
u0 = V0(:,2);
u0 = u0/sqrt(u0'*M*u0);         % L2 normalization
figure(1); clf; plotf2(mesh0,u0); xlabel('$x$'); ylabel('$y$');

%% compute u0' via finite differencing in the eigenvectors
e = 1e-6;
[Ve,lame] = solve_CG(mesh0, @(x) dLx(a+e,x), deg); 
lam0p_fd = (lame(2)-lam0(2))/e
ue = Ve(:,2); 
ue = ue/(ue'*M*u0);           % determine scaling of u_e such that in V_0
u0p_fd = (ue-u0)/e;
figure(2); clf; plotf2(mesh0,u0p_fd); xlabel('$x$'); ylabel('$y$');

%% compute u0' via linear response (CG approach)
Gp = triquad(mesh0,@(x) dLpx(a,x),deg); 
L = assemble2(mesh0,Gp);        
lam0p = (u0'*L*u0)/(u0'*M*u0)
udot1 = [(K-lam0(2)*M) u0; u0' 0]\[((lam0p*M-L)*u0)+u0; 0]; 
udot = udot1(1:end-1);
figure(3); clf; plotf2(mesh0,udot); xlabel('$x$'); ylabel('$y$');

%% H1 error of finite difference approximation of linear response
Id = @(a,x) eye(2);
[~,~,KId] = solve_CG(mesh0, @(x) Id(a,x), deg); 
H1_error = sqrt((u0p_fd-udot)'*(M-KId)*(u0p_fd-udot))
eig_error = abs(lam0p_fd-lam0p)

%% comparison of different methods to obtain coherent set prediction
% exact second eigenvector of perturbed map
da = 0.5; 
Ve = solve_CG(mesh0, @(x) dLx(a+da,x), deg);
ue = Ve(:,2);           
figure(4); clf; plotf2(mesh0,ue); xlabel('$x$'); 
% approximation by linear response
ulr = u0 + da*udot;     
figure(5); clf; plotf2(mesh0,ulr); xlabel('$x$'); 

%% velocity field for level set evolution, level set evolution
figure(6); clf
x = unique(p0(:,1)); y = unique(p0(:,2));
[X,Y] = meshgrid(x,y);
U0 = reshape(u0,nx,ny); Udot = reshape(udot,nx,ny); Ue = reshape(ue,nx,ny);
[gx,gy] = gradient(U0,2*pi/nx,2*pi/ny);
I = 1:nx/20:nx;
normg = sqrt(gx(I,I).^2+gy(I,I).^2);
quiver(X(I,I),Y(I,I),-gx(I,I).*Udot(I,I)./normg,-gy(I,I).*Udot(I,I)./normg);
view(0,90); axis equal; axis tight; xlabel('$x$'); ylabel('$y$');
level = -0.2;
hold on; contour(X,Y,U0,[level level],'k','linewidth',1);
contour(X,Y,Ue,[level level],'k--','linewidth',1);

%% estimate of the linear response via Taylor approximation of the flow
da = 0.5; 
Tt   = @(a1,x) T(a,x) + (a1-a)*Tp(a,x);
DTt  = @(a,x) D_fd(@(x) Tt(a,x),x);                           
dLtx  = @(a,x) fapply1(dL, DTt(a,x));
Vt = solve_CG(mesh0, @(x) dLtx(a+da,x), deg);
ut = Vt(:,2);         
figure(7); clf; plotf2(mesh0,ut); xlabel('$x$'); 

      