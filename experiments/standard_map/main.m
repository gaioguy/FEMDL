addpath('../../src'); clear all; colormap jet;
dott = @(A) A*A';
symmetric = @(A) 0.5*(A+A');
mymod = @(x,m) mod(real(x),m) + i*mod(imag(x),m);

%% modified standard map
% a = 0.971635;
% %a = 20;
% T = @(a,x) mod([x(:,1) + x(:,2) + a*sin(x(:,1)+x(:,2)), x(:,2) + a*sin(x(:,1))],2*pi);
% DT = @(a,x) [1 + a*cos(x(1)+x(2)) 1 + a*cos(x(1)+x(2));
%              a*cos(x(1)) 1];
% w = @(a,x) [sin(x(:,1)+x(:,2)), sin(x(:,1))];
% Dw = @(a,x) [ cos(x(1)+x(2)), cos(x(1)+x(2));
%               cos(x(1)), 0];

%% standard map
a = 0.971635;
T = @(a,x) mymod([x(:,1) + x(:,2) + a*sin(x(:,1)), x(:,2) + a*sin(x(:,1))],2*pi);
DT = @(a,x) [1 + a*cos(x(1)) 1;
                 a*cos(x(1)) 1];
w = @(a,x) [sin(x(:,1)), sin(x(:,1))];
Dw = @(a,x) [ cos(x(1)), 0;
              cos(x(1)), 0];

%% regular triangulation
nx = 200;  ny = nx;  n = (nx-1)*(ny-1); dx = 2*pi/(nx-1); dy = dx;
[xi,yi] = meshgrid(linspace(0,2*pi-dx,nx-1),linspace(0,2*pi-dy,ny-1));
p0 = [xi(:) yi(:)];
[p,t,pb] = delaunay_T2(p0,[0,2*pi,0,2*pi]);

%% compute u0' via finite differencing in the eigenvectors
deg = 2;
e = 1e-3;
A0 = @(a,x) 0.5*(eye(2) + dott(inv(DT(a,x))));
[~,u0,lambda0,K0,M] = get2nd_ev(p,t,pb,A0,a,deg);
u0 = u0/sqrt(u0'*M*u0);         % L2 normalization
[~,ue,lambdae,~,~] = get2nd_ev(p,t,pb,A0,a+e,deg);
lambda0p_fd = (lambdae-lambda0)/e
ue = ue/(ue'*(M)*u0);           % determine scaling of u_e such that in V_0
u0p_fd = (ue-u0)/e;
figure(1); clf; plotf(p,t,pb,u0p_fd,0); colormap jet; title('finite differencing');

%% compute u0' via linear response (CG approach)
Ap = @(a,x) -symmetric(inv(DT(a,x))*Dw(a,x)*dott(inv(DT(a,x))));
Gp = inv_CG_quad(p,t,@(x) rowfun(@(x) Ap(a,x),x), deg);        
[L,~] = assemble(p,t,pb,Gp);        
lambda0p = (u0'*L*u0)/(u0'*M*u0)
v = [(K0-lambda0*M) u0; u0' 0]\[((lambda0p*M-L)*u0)+u0; 0]; 
u0p = v(1:end-1);
figure(2); clf; plotf(p,t,pb,u0p,0); colormap jet; title('CG');

%% compute u0' via linear response (TO approach)
% assembly on initial domain
I = kron([1 0 1],ones(size(t,1),1));
[K,M] = assemble(p,t,pb,I);        
% assembly on image domain
Tp0 = T(a,p0);
[p1,t1,pb1] = delaunay_T2(Tp0,[0,2*pi,0,2*pi]);
I = kron([1 0 1],ones(size(t1,1),1));
[K1,~] = assemble(p1,t1,pb1,I);
% solve eigenproblem for dynamic Laplacian
[V,lambda] = eigs(0.5*(K+K1),M,6,'SM'); 
[lambda,ord] = sort(diag(lambda),'descend'); 
lambda0 = lambda(2);
u0 = V(:,2); 
% compute linear response
W = w(a,p0);
L = assemble_lr(p1,t1,pb1,W); 
L = 0.5*(L+L');
lambda0p = (u0'*L*u0)/(u0'*M*u0)
v = [(K0-lambda0*M) u0; u0' 0]\[((lambda0p*M-L)*u0)+u0; 0]; 
u0p_TO = v(1:end-1);
figure(3); clf; plotf(p,t,pb,u0p_TO,0); colormap jet; title('TO');

%% H1 error of finite difference approximation of linear response
Id = @(a,x) eye(2);
[~,~,~,K,~] = get2nd_ev(p,t,pb,Id,a,deg); 
H1_error = sqrt((u0p_fd-u0p)'*(M-K)*(u0p_fd-u0p))
eig_error = abs(lambda0p_fd-lambda0p)


