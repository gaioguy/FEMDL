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
u0 = V0(:,2); u0 = u0/sqrt(u0'*M*u0);                   % L2 normalization
figure(1); clf; plotf2(m0,u0);  colorbar

%% compute udot via finite differencing
e = 1e-4; 
[Ve,lame] = solve_CG(m0,@(x) dLx(tf+e,x),deg); 
lamdot_fd = (lame(2)-lam0(2))/e
ue = Ve(:,2);
ue = ue/(ue'*M*u0);             % determine scaling of u_e such that in V_0
udot_fd = (ue-u0)/e;
figure(2); plotf2(m0,udot_fd); colorbar; 

%% compute udot by CG approach
Aq = triquad(m0,@(x) dLpx(tf,x),deg); 
L = assemble2(m0,Aq);  
lamdot = (u0'*L*u0)/(u0'*M*u0)
udot1 = [(K-lam0(2)*M) u0; u0' 0]\[((lamdot*M-L)*u0)+u0; 0]; 
udot = udot1(1:end-1);
figure(3); plotf2(m0,udot); colorbar; 

%% error of FD approximation of linear response
Id = @(x) eye(2);
[~,~,KId] = solve_CG(m0, Id, deg); 
H1_error = sqrt((udot_fd-udot)'*(M-KId)*(udot_fd-udot))
eig_error = abs(lamdot_fd-lamdot)

%% prediction of u_eps by linear response
dt = 0.2; 
u0lr = u0 + dt*udot;
figure(4); clf; plotf2(m0,u0lr); colorbar; 

%% exact u_eps
Ve = solve_CG(m0, @(x) dLx(tf+dt,x), deg);
ue = Ve(:,2);
figure(5); clf; plotf2(m0,ue); colorbar; 

%% velocity field for level set evolution 
x = unique(p(:,1)); y = unique(p(:,2));
[X,Y] = meshgrid(x,y);
U0 = reshape(u0,nx,ny);
U0P = reshape(udot,nx,ny);
[gx,gy] = gradient(U0,1/nx,1/ny);
I = 1:6:nx;
normg = sqrt(gx(I,I).^2+gy(I,I).^2);
figure(6); clf
quiver(X(I,I),Y(I,I),-gx(I,I).*U0P(I,I)./normg,-gy(I,I).*U0P(I,I)./normg);
axis equal; axis tight; xlabel('$x$'); ylabel('$y$');
level = 0.5;
hold on; contour(X,Y,U0,level*[1 1],'k','linewidth',1);
contour(X,Y,U0+dt*U0P,level*[1 1],'k--','linewidth',1);

%% prediction of u_eps by Taylor approximation of the flow
dtvf = @(t,x) df_fd(@(t) vf(t,x), t, 1)';
vft  = @(t,x) (t<tf).*vf(t,x) + (t>=tf).*(vf(tf,x) + (t-tf).*dtvf(tf,x));
Tt   = @(a,x) at_tf(flowmap(@(t,x) vft(t,x), x, [t0 a]));
DTt  = @(a,x) D(@(x) Tt(a,x),x);                           
dLtx  = @(a,x) fapply1(dL, DTt(a,x));

V0t = solve_CG(m0, @(x) dLtx(tf+dt,x), deg);
u0t = V0t(:,2);
figure(7); clf; plotf2(m0,u0t); colorbar; 

