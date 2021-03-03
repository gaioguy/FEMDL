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
deg = 5;                   % degree of Gauss quadrature for space integrals
k = 2;                     % number of eigenfunction to consider

%% compute u0
[V0,lam0,K,M] = solve_CG(mesh0, @(x) dLx(a,x), deg); 
lambda0 = lam0(k)
u0 = V0(:,k);
u0 = u0/sqrt(u0'*M*u0);         % L2 normalization
figure(1); clf; plotf2(mesh0,u0); xlabel('$x$'); ylabel('$y$');

%% compute u0' via finite differencing in the eigenvectors
e = 1e-6;
[Ve,lame] = solve_CG(mesh0, @(x) dLx(a+e,x), deg); 
lamdot0_fd = (lame(k)-lam0(k))/e
ue = Ve(:,k); 
ue = ue/(ue'*M*u0);           % determine scaling of u_e such that in V_0
udot0_fd = (ue-u0)/e;
figure(2); clf; plotf2(mesh0,udot0_fd); xlabel('$x$'); ylabel('$y$');

%% compute udot0 by CG approach
Gp = triquad(mesh0,@(x) dLpx(a,x),deg); 
L = assemble2(mesh0,Gp);        
tmp = [(K-lam0(k)*M) -M*u0; (M*u0)' 0]\[-L*u0; 0]; 
udot0 = tmp(1:end-1); 
lamdot0 = tmp(end); 
figure(4); clf; plotf2(mesh0,udot0); xlabel('$x$'); %ylabel('$y$');

%% H1 error of finite difference approximation of linear response
Id = @(a,x) eye(2);
[~,~,KId] = solve_CG(mesh0, @(x) Id(a,x), deg); 
H1_error = sqrt((udot0_fd-udot0)'*(M-KId)*(udot0_fd-udot0))
eig_error = abs(lamdot0_fd-lamdot0)

%% comparison of different methods to obtain coherent set prediction
% exact second eigenvector of perturbed map
da = 0.5; 
[Ve, lame] = solve_CG(mesh0, @(x) dLx(a+da,x), deg); 
lambdae = lame(k)
ue = Ve(:,k);      
ue = ue/sqrt(ue'*M*ue);         % L2 normalization
figure(4); clf; plotf2(mesh0,ue); xlabel('$x$'); 
% approximation by linear response
ulr = u0 + da*udot0; 
ulr = ulr/sqrt(ulr'*M*ulr);         % L2 normalization
figure(5); clf; plotf2(mesh0,ulr); xlabel('$x$'); 
L2_error = sqrt((ue-ulr)'*M*(ue-ulr))

%% compare dynamic isoperimetric ratio of some level set of u0 and ue
level = 0.2;
%dipr0 = dipr_of_levelset(mesh0, @(x) T(a,x), u0, level)
% dipre = dipr_of_levelset(mesh0, @(x) T(a+da,x), ue, level);
[c,minipr,minc,iprs] = curve_of_optimal_ipr(mesh0,@(x)T(a,x),u0,0,max(u0));
level = minc;

%% velocity field for level set evolution, level set evolution
figure(6); clf
x = unique(p0(:,1)); y = unique(p0(:,2));
[X,Y] = meshgrid(x,y);
U0 = reshape(u0,nx,ny); Udot = reshape(udot0,nx,ny); Ue = reshape(ue,nx,ny);
[gx,gy] = gradient(U0,2*pi/nx,2*pi/ny);
I = 1:nx/20:nx;
normg = sqrt(gx(I,I).^2+gy(I,I).^2);
quiver(X(I,I),Y(I,I),-gx(I,I).*Udot(I,I)./normg,-gy(I,I).*Udot(I,I)./normg);
view(0,90); axis equal; axis tight; xlabel('$x$'); ylabel('$y$');
c = 0.2;
hold on; 
contour(X,Y,U0,[c c],'k','linewidth',2);
contour(X,Y,U0+da*Udot,[c c],'g','linewidth',2);
contour(X,Y,Ue,[c c],'r--','linewidth',2);

%% estimate of the linear response via Taylor approximation of the flow
da = 0.5; 
Tt   = @(a1,x) T(a,x) + (a1-a)*Tp(a,x);
DTt  = @(a,x) D_fd(@(x) Tt(a,x),x);                           
dLtx  = @(a,x) fapply1(dL, DTt(a,x));
Vt = solve_CG(mesh0, @(x) dLtx(a+da,x), deg);
ut = Vt(:,k);         
figure(7); clf; plotf2(mesh0,ut); xlabel('$x$'); 


      